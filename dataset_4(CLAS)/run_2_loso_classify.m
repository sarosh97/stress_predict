clc
clear all
close all

%% 1. LOSO classification

algo_type = 'DT';    %KNN, RF, DT, LDA, AB

%No need to define numclasses since this dataset is only used for binary
%classification

%"features_30.mat" refers to features computed using 30s window shift. Same
%for 5s shift. For 0.25s shift we use "features.mat"
features = load('features_5.mat');

subj = 1:60;
subj(4) = [];

for i = 1:length(subj)
    X_train = [];
    Y_train = [];
    for j = 1:length(subj)
        if j ~= i
            subj_name = ['S', num2str(subj(j))];
            
            %only EDA features
            X_train = [X_train; features.features.(subj_name).feat];
            Y_train = [Y_train; features.features.(subj_name).labels];
        end
    end
    
    subj_name_test = ['S', num2str(subj(i))];
    
    X_test = features.features.(subj_name_test).feat;
    
    Y_test = features.features.(subj_name_test).labels;
    
    
    Y_pred = classification_function(X_train, Y_train, X_test, Y_test, algo_type); 
    
    accuracies(i) = sum(Y_pred == Y_test) / numel(Y_test);

    %only for binary classification
    confmat = confusionmat(Y_test, Y_pred);
    TP = confmat(2,2);
    TN = confmat(1,1);
    FP = confmat(1,2);
    FN = confmat(2,1);

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
function Y_pred = classification_function(X_train, Y_train, X_test, Y_test, algo_type)

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
    mdl_ab = fitcensemble(X_train, Y_train, 'Method', 'AdaBoostM1',...
                          'NumLearningCycles', numTrees, 'Learners', tr);
    Y_pred = predict(mdl_ab, X_test);
else
    disp('Wrong selection/error')
end

end
