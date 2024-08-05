%The data was available in .pkl format which was first converted to .mat
%using a python script

%Labels are numbered as baseline=0, amusement(fun)=1, meditation=2, 
%stress(tsst)=3 

function [data, labels_num, start_time, end_time, subjective_labels] = load_wesad_data(filename)

%Loading the Data and labels
subj = 2:17;    
subj(11) = [];      %subjects' filenames from 2 to 17 excluding 12
data = cell(1,15);          %Data to be stored in cell arrays

labels_num = nan(15,5);       %labels data, as numbers, stored as matrix
start_time = nan(15,5);   %Rows=subjects, columns=start time of each condition
end_time = nan(15,5);     %Rows=subjects, columns=end time of each condition

for i = 1:length(subj)
    temp1 = load([filename '\S' num2str(subj(i)) '\S' num2str(subj(i)) '.mat']);
    data{i} = temp1.pickle_data.signal;
    
    label_data = importdata([filename '\S' num2str(subj(i)) '\S' num2str(subj(i)) '_quest.csv']);
    
    labels_txt = label_data.textdata(2,2:6);     %Baseline, fun, stress, meditation (currently in string format, need to assign a number to each condition)
    
    %temp2 = ~isnan(label_data.data(1,1:5));    %temporary variable to remove NaN
    start_time(i,1:length(labels_txt)) = label_data.data(1,1:5)*60;   %start time of each test/condition, in seconds
    
    %temp2 = ~isnan(label_data.data(2,1:5));    %temporary variable to remove NaN
    end_time(i,1:length(labels_txt)) = label_data.data(2,1:5)*60;     %end time of each test/condition, in seconds
    
    %mapping text label to numbers
    labels_num(i,1:length(labels_txt)) = txt_label_to_num(labels_txt);
    
    %subjective labels (mean of SSSQscores divided by 5 to normalize between 0 and 1)
    subjective_labels(i,1) = mean(label_data.data(end,1:6))/5;
end
end


function labels_num = txt_label_to_num(labels_txt)

%Labels are numbered as baseline=0, amusement(fun)=1, meditation=2, stress(tsst)=3 

for i = 1:length(labels_txt)
    if strcmp(labels_txt{i}, 'Base')
        labels_num(i) = 0;
    elseif strcmp(labels_txt{i}, 'Fun')
        labels_num(i) = 1;
    elseif strcmp(labels_txt{i}, 'Medi 1') || strcmp(labels_txt{i}, 'Medi 2')
        labels_num(i) = 2;
    elseif strcmp(labels_txt{i}, 'TSST')
        labels_num(i) = 3;
    else
        labels_num(i) = nan;
    end
end
end