clc
clear all

fs_e = 4;   %Sampling frequency of EDA

%Defining window length for feature_extraction portion 
win_len_e = 60*fs_e;

%Defining window shift
win_shift_e = 5*fs_e;

%since location of tag is different in each file, I am manually providing
%it here. The order is subject-wise.
%Not taking some subjects (6, 17, 18, 29)
tags_loc = [4, 3, 3, 3, 3, nan, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, nan, nan, ...
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3];
tags_loc([6,17,18, 29,30]) = [];

subj = strsplit(sprintf('%02d ',1:35));
subj([6,17,18,29,30,36]) = [];

for i = 1:length(subj)
    
    %time log
    tags = importdata(['D:\Research Papers and Data\Wearable sensing\dataset_2\Raw_data\S' subj{i} '\tags_S' subj{i} '.csv']);
    tags = round(tags);

    %Indices for the duration of stress test
    n_e = [((tags(2)-tags(1))*fs_e:(tags(3)-tags(1))*fs_e-1)]';

    data_raw_e = importdata(['D:\Research Papers and Data\Wearable sensing\dataset_2\Raw_data\S' subj{i} '\EDA.csv']);

    %first value of each file is timestamp (not same as first timestamp of 'logs')
    data_log_e = data_raw_e(1);

    %Selecting the data that starts from the timestamp given in 'logs' file
    data_e = data_raw_e(round((tags(1)-data_log_e))*fs_e:end,1);

    %Selecting only the stress data and detrending it (except TEMP and EDA) 
    data_stress_e = (data_e(n_e,:));

    %defining time in seconds
    t_e = n_e/fs_e;
 
    %EDA feature extraction
    current_idx = 1;
    end_idx = n_e(end) - n_e(1) - win_len_e;
    M = 1;
    while (current_idx < end_idx )        %sliding windows and computing features for each condition
        current_win = [current_idx:current_idx + win_len_e - 1];  %indices for the current window
        time_win = t_e(current_win,1);
        signal = data_stress_e(current_win,1);
        feat_eda(M,:) = feat_ext_eda(signal, time_win, fs_e);
    
        current_idx = current_idx + win_shift_e;
        M = M+1;
    end

    %features saved as struct
    subj_name = ['S', subj{i}];
    features.(subj_name) = feat_eda;
    
    disp(['Features of S' subj{i} ' computed'])
end

if win_shift_e == 0.25*fs_e
    save('features_relax.mat', 'features')

elseif win_shift_e == 5*fs_e
    save('features_relax_5.mat', 'features')

elseif win_shift_e == 30*fs_e
    save('features_relax_30.mat', 'features')
end
