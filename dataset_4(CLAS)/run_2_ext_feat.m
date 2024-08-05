clc
clear all
close all

%extracting only EDA features
filename = 'D:\Research Papers and Data\Wearable sensing\CLAS (dataset_4)';

fs_e = 256;   %Sampling frequency of EDA
new_fs = 16;

win_len_e = 60*new_fs;
win_shift_e = 5*new_fs;
        
subj = 1:60;
subj(4) = [];

for i = 1:length(subj)
    
    test_name = readtable([filename '\Block_details\Part' num2str(subj(i)) '_Block_Details.csv']);
    idx = contains(test_name.BlockType,'Neutral');
    
    temp = 1:length(idx);
    neutral_files = temp(idx);
    
    M = 1;
    data_e_relax_big = [];
    %relax phase feature extraction
    %First concatenating the "neutral" in the entire dataset of a single
    %subject and then extracting features
    for j = neutral_files    
        
        data_raw = readtable([filename '\Participants\Part' num2str(subj(i)) '\by_block\' test_name.EDA_PPGFile{j}]);
        data = data_raw.gsr;
        data_e_relax_big = [data_e_relax_big; data]; 
    end
    
    data_resampled = resample(data_e_relax_big, new_fs, fs_e); %resampling for quicker processing
    data_e_relax = data_resampled(new_fs:end-new_fs); %chopping off the ends (about 1 fs_new) to remove transients
    t_e_relax = (0:length(data_e_relax)-1)'/new_fs;
    
    current_idx = 1;
    end_idx = length(t_e_relax) - win_len_e - 1;
    while (current_idx < end_idx )        %sliding windows and computing features for each condition
        current_win = current_idx:current_idx + win_len_e - 1;  %indices for the current window
        
        time_win = t_e_relax(current_win,1);
        signal_e = data_e_relax(current_win,:);
        feat_mat(M,:) = feat_ext_eda(signal_e, time_win, new_fs);
        
        label_mat(M,1) = 0;
    
        current_idx = current_idx + win_shift_e;
        M = M+1;
    end
            
    
    %Stress phase feature extraction
    idx1 = contains(test_name.BlockType,'Math Test');
    idx2 = contains(test_name.BlockType,'Stroop Test');
    idx3 = contains(test_name.BlockType,'IQ Test');
    
    temp = 1:length(idx1); %or idx2 or idx3, same number of elements
    stress_files = temp(idx1 | idx2 | idx3);
  
    for j = stress_files    
        if strcmp(test_name.BlockType{j}, 'Math Test')
            data_raw = readtable([filename '\Participants\Part' num2str(subj(i)) '\by_block\' test_name.EDA_PPGFile{j}(1:end-4) 'mathtest.csv']);
            data = data_raw.gsr;
            data_resampled = resample(data, new_fs, fs_e); %resampling for quicker processing
            data_e_stress = data_resampled(new_fs:end-new_fs); %chopping off the ends (about 1 fs_new) to remove transients
            t_e_stress = (0:length(data_e_stress)-1)'/new_fs;

            %Stress phase feature extraction
            current_idx = 1;
            end_idx = length(t_e_stress) - win_len_e - 1;
            while (current_idx < end_idx )        %sliding windows and computing features for each condition
                current_win = current_idx:current_idx + win_len_e - 1;  %indices for the current window
        
                time_win = t_e_stress(current_win,1);
                signal_e = data_e_stress(current_win,:);
                feat_mat(M,:) = feat_ext_eda(signal_e, time_win, new_fs);
        
                label_mat(M,1) = 1;
    
                current_idx = current_idx + win_shift_e;
                M = M+1;
            end
            
        elseif strcmp(test_name.BlockType{j}, 'Stroop Test')
            data_raw = readtable([filename '\Participants\Part' num2str(subj(i)) '\by_block\' test_name.EDA_PPGFile{j}(1:end-4) 'strooptest.csv']);
            data = data_raw.gsr;
            data_resampled = resample(data, new_fs, fs_e); %resampling for quicker processing
            data_e_stress = data_resampled(new_fs:end-new_fs); %chopping off the ends (about 1 fs_new) to remove transients
            t_e_stress = (0:length(data_e_stress)-1)'/new_fs;

            %Stress phase feature extraction
            current_idx = 1;
            end_idx = length(t_e_stress) - win_len_e - 1;
            while (current_idx < end_idx )        %sliding windows and computing features for each condition
                current_win = current_idx:current_idx + win_len_e - 1;  %indices for the current window
        
                time_win = t_e_stress(current_win,1);
                signal_e = data_e_stress(current_win,:);
                feat_mat(M,:) = feat_ext_eda(signal_e, time_win, new_fs);
        
                label_mat(M,1) = 1;
    
                current_idx = current_idx + win_shift_e;
                M = M+1;
            end

        elseif strcmp(test_name.BlockType{j}, 'IQ Test')
            data_raw = readtable([filename '\Participants\Part' num2str(subj(i)) '\by_block\' test_name.EDA_PPGFile{j}(1:end-4) 'IQtest.csv']);
            data = data_raw.gsr;
            data_resampled = resample(data, new_fs, fs_e); %resampling for quicker processing
            data_e_stress = data_resampled(new_fs:end-new_fs); %chopping off the ends (about 1 fs_new) to remove transients
            t_e_stress = (0:length(data_e_stress)-1)'/new_fs;

            %Stress phase feature extraction
            current_idx = 1;
            end_idx = length(t_e_stress) - win_len_e - 1;
            while (current_idx < end_idx )        %sliding windows and computing features for each condition
                current_win = current_idx:current_idx + win_len_e - 1;  %indices for the current window
        
                time_win = t_e_stress(current_win,1);
                signal_e = data_e_stress(current_win,:);
                feat_mat(M,:) = feat_ext_eda(signal_e, time_win, new_fs);
        
                label_mat(M,1) = 1;
    
                current_idx = current_idx + win_shift_e;
                M = M+1;
            end
        end
    end
    

    %features saved as struct
    subj_name = ['S', num2str(subj(i))];
    features.(subj_name).feat = feat_mat;
    features.(subj_name).labels = label_mat;
    
    disp(['Features of S' num2str(subj(i)) ' computed'])
end

if win_shift_e == 0.25*new_fs
    save('features.mat', 'features')

elseif win_shift_e == 5*new_fs
    save('features_5.mat', 'features')

elseif win_shift_e == 30*new_fs
    save('features_30.mat', 'features')
end

