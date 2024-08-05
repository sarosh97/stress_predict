clc
clear all
close all

%% 1. Computing feature matrix and Label matrix(vector) for LOSO classification

algo_type = 'RF';    %KNN, RF, DT, LDA, AB
numclasses = 2;      %2, 3, 4

features = load('features_60.mat');

subj = 1:15;          %all subjects

mod = {'acc', 'ecg', 'emg60', 'resp', 'temp', 'acc', 'bvp', 'temp', 'eda'}; 

%initializing minimum number of rows to infinite value
%Done to find minimum number of rows of each class, of each subject
min_row0(subj) = inf;
min_row1(subj) = inf;
min_row2(subj) = inf;
min_row3(subj) = inf;

%Loop for saving features in cells
for i = 1:length(subj)    %loop through subjects (subj_all or subj_m or subj_f - also change the value at line 90)
    for j = 1:length(mod)   %loop through modalities
        subj_name = ['S', num2str(i)];
        
        %6 feature modalties from chest and 4 from wrist
        if j<=5
            X = features.feat.(subj_name).chest.(mod{j}).feat;
            Y = features.feat.(subj_name).chest.(mod{j}).labels;   
        else
            X = features.feat.(subj_name).wrist.(mod{j}).feat;
            Y = features.feat.(subj_name).wrist.(mod{j}).labels;
        end    
        
        idx0 = find(Y==0);  %indices where the class is "0"
        idx1 = find(Y==1);  %indices where the class is "1"
        idx2 = find(Y==2);  %indices where the class is "2"
        idx3 = find(Y==3);  %indices where the class is "3"
        
        class0_feat{i,j} = X(idx0,:);   %i=subjects, j=modalities
        class1_feat{i,j} = X(idx1,:);   %i=subjects, j=modalities
        class2_feat{i,j} = X(idx2,:);   %i=subjects, j=modalities
        class3_feat{i,j} = X(idx3,:);   %i=subjects, j=modalities
        
        class0_label{i,j} = Y(idx0,1);   %i=subjects, j=modalities
        class1_label{i,j} = Y(idx1,1);   %i=subjects, j=modalities
        class2_label{i,j} = Y(idx2,1);   %i=subjects, j=modalities
        class3_label{i,j} = Y(idx3,1);   %i=subjects, j=modalities
        
        %to see, for each class and each subject, which is smallest number of rows so that
        %all features of a class, of a subject, have same number of rows
        current_row0 = length(idx0);
        current_row1 = length(idx1);
        current_row2 = length(idx2);
        current_row3 = length(idx3);
        
        if current_row0 < min_row0(i)
            min_row0(i) = current_row0;
        end
        if current_row1 < min_row1(i)
            min_row1(i) = current_row1;
        end
        if current_row2 < min_row2(i)
            min_row2(i) = current_row2;
        end
        if current_row3 < min_row3(i)
            min_row3(i) = current_row3;
        end 
    end
end


%loop to compute feature matrix and label matrix(column vector) where columns are features of all
%modalities and rows are computed features
%Label matrix, and corresponding feature matrix, computed like [0,1,2,3,0,1,2,3,0,1,2,3....]
%Features and labels of classes are saved for each subject separately

k = [15, 18, 5, 12, 5, 15, 18, 5, 28];  %Number of features from each modality
n = 1;      %for feature/label rows
for i = 1:length(subj)      %loop through subject 
    subj_name = ['S', num2str(i)];
    m = 1;  %for feature/label columns
    for j = 1:length(mod)   %loop through modalities
        feat0 = class0_feat{i,j}(1:min_row0(i),:);
        feat1 = class1_feat{i,j}(1:min_row1(i),:);
        feat2 = class2_feat{i,j}(1:min_row2(i),:);
        feat3 = class3_feat{i,j}(1:min_row3(i),:);
        
        feat_col = m : m+k(j)-1;  %done to concatenate columns of each modality
        feat_mat.(subj_name)(:,feat_col) = [feat0; feat1; feat2; feat3];
        m = m + k(j);             %update column number
    end
    
    label0 = 0*ones(min_row0(i),1);
    label1 = 1*ones(min_row1(i),1);
    label2 = 2*ones(min_row2(i),1);
    label3 = 3*ones(min_row3(i),1);
    
    label_mat.(subj_name) = [label0; label1; label2; label3];

end

%% 2. LOSO classification

subj = 15;
for i = 1:subj
    X_train = [];
    Y_train = [];
    for j = 1:subj
        if j ~= i
            subj_name = ['S', num2str(j)];
            X_train = [X_train; feat_mat.(subj_name)];
            Y_train = [Y_train; label_mat.(subj_name)];
        end
    end
    
    X_test = feat_mat.(['S' num2str(i)]); %current subject index as test matrix
    Y_test = label_mat.(['S' num2str(i)]);
    
    %to remove any imaginary numbers by taking abs value
    [row_train,col_train] = find(imag(X_train));
    [row_test,col_test] = find(imag(X_test));
    X_train(row_train, col_train) = abs(X_train(row_train, col_train));
    X_test(row_test, col_test) = abs(X_test(row_test, col_test));
    

    [Y_pred, Y_test_new] = classification_function(X_train, Y_train, X_test, Y_test, algo_type, numclasses);

    accuracies(i) = sum(Y_pred == Y_test_new) / length(Y_test_new);

    %f1 score computation only for "stress" class
    confmat = confusionmat(Y_test_new, Y_pred);
    if length(unique(Y_test_new)) == 2
        TP = confmat(2,2);
        TN = confmat(1,1);
        FP = confmat(1,2);
        FN = confmat(2,1);
    elseif length(unique(Y_test_new)) == 3
        TP = confmat(3,3);
        TN = sum(confmat(1:2,1:2), 'all');
        FP = sum(confmat(1:2,end), 'all');
        FN = sum(confmat(end,1:2), 'all');
    elseif length(unique(Y_test_new)) == 4
        TP = confmat(4,4);
        TN = sum(confmat(1:3,1:3),'all');
        FP = sum(confmat(1:3,end), 'all');
        FN = sum(confmat(end,1:3), 'all');
    end

        prec = TP/(TP+FP);
        rec = TP/(TP+FN);
        f1(i) = 2*(prec*rec)/(prec+rec);

    i
end

accuracy = mean(accuracies)*100
acc_range = std(accuracies)*100

%in case TP=0, f1 score cannot be computed so we simply put f1=0 for that
%part
f1(isnan(f1)) = 0;

f1score = mean(f1)*100
f1range = std(f1)*100



%Custom function for selection of classifier and number of classes
function [Y_pred, Y_test_new] = classification_function(X_train, Y_train, X_test, Y_test, algo_type, numclasses)

if numclasses == 2
    Y_train(Y_train<3) = 0;
    Y_train(Y_train==3) = 1;
    Y_test(Y_test<3) = 0;
    Y_test(Y_test==3) = 1;
    ABtype = 'AdaBoostM1';
elseif numclasses == 3
    X_train(Y_train==2,:) = [];
    Y_train(Y_train==2,:) = [];
    X_test(Y_test==2,:) = [];
    Y_test(Y_test==2,:) = [];
    Y_train(Y_train==3) = 2;
    Y_test(Y_test==3) = 2;
    ABtype = 'AdaBoostM2';
elseif numclasses == 4
    X_train = X_train;
    Y_train = Y_train;
    ABtype = 'AdaBoostM2';
end

Y_test_new = Y_test;

if strcmp(algo_type, 'KNN')
    k = 9;
    mdl_knn = fitcknn(X_train, Y_train, 'NumNeighbors', k);
    Y_pred = predict(mdl_knn, X_test);
elseif strcmp(algo_type, 'RF')
    numTrees = 100;
    mdl_rf = TreeBagger(numTrees, X_train, Y_train, 'Method', 'classification'); 
    Y_pred_str = predict(mdl_rf, X_test);
    Y_pred = cellfun(@str2double, Y_pred_str);
elseif strcmp(algo_type, 'DT')
    mdl_dt = fitctree(X_train, Y_train);
    Y_pred = predict(mdl_dt, X_test);
elseif strcmp(algo_type, 'LDA')
    mdl_lda = fitcdiscr(X_train, Y_train, 'discrimType','pseudoLinear');
    Y_pred = predict(mdl_lda, X_test);
elseif strcmp(algo_type, 'AB')
    numTrees = 100;
    tr = templateTree('MaxNumSplits',20);
    mdl_ab = fitcensemble(X_train, Y_train, 'Method', ABtype,...
                          'NumLearningCycles', numTrees, 'Learners', tr);
    Y_pred = predict(mdl_ab, X_test);
else
    disp('Wrong selection/error')
end

end