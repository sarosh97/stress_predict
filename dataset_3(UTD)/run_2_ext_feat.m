clc
clear all
close all

%extracting only EDA features
filename = 'D:\Research Papers and Data\Wearable sensing\UTD (dataset_3)';

fs_e = 8;   %Sampling frequency of EDA

%Defining window length for feature_extraction portion
win_len_e = 60*fs_e;

%Defining window shift
win_shift_e = 0.25*fs_e;

subj = 1:20;

for i = 1:length(subj)
    
    data_1 = readtable([filename '\Subject' num2str(subj(i)) '\Subject' num2str(subj(i)) 'AccTempEDA.csv']);
    
    %Data_1 from Affectiva (Only emotional or cognitive stress)
    relax_idx = contains(data_1.Label,'Relax');
    stress_idx = [contains(data_1.Label,'CognitiveStress') | contains(data_1.Label,'PhysicalStress')];
    
    data_e = data_1.EDA;
    data_e_relax = data_e(relax_idx,:);
    data_e_stress = data_e(stress_idx,:);
    
    %defining time in seconds 
    t_e_relax = (0:length(data_e_relax)-1)'/fs_e;
    t_e_stress = (0:length(data_e_stress)-1)'/fs_e;

    %Relax phase feature extraction
    current_idx = 1;
    end_idx = length(t_e_relax) - win_len_e - 1;
    M = 1;
    while (current_idx < end_idx )        %sliding windows and computing features for each condition
        current_win = [current_idx:current_idx + win_len_e - 1];  %indices for the current window
        
        time_win = t_e_relax(current_win,1);
        signal_e = data_e_relax(current_win,:);
        feat_eda = feat_ext_eda(signal_e, time_win, fs_e);
        
        feat_mat(M,:) = feat_eda;
        label_mat(M,1) = 0;
    
        current_idx = current_idx + win_shift_e;
        M = M+1;
    end

    %Stress phase feature extraction
    current_idx = 1;
    end_idx = length(t_e_stress) - win_len_e - 1;
    while (current_idx < end_idx )        %sliding windows and computing features for each condition
        current_win = [current_idx:current_idx + win_len_e - 1];  %indices for the current window
        
        time_win = t_e_stress(current_win,1);
        signal_e = data_e_stress(current_win,:);
        feat_eda = feat_ext_eda(signal_e, time_win, fs_e);
        
        feat_mat(M,:) = feat_eda;
        label_mat(M,1) = 1;
    
        current_idx = current_idx + win_shift_e;
        M = M+1;
    end
    
    %features saved as struct
    idx_val = ~isnan(feat_mat(:,1));   %to remove NaNs
    feat_mat = feat_mat(idx_val,:);
    label_mat = label_mat(idx_val,1);
    subj_name = ['S', num2str(subj(i))];
    features.(subj_name).feat = feat_mat;
    features.(subj_name).labels = label_mat;
    
    disp(['Features of S' num2str(subj(i)) ' computed'])
end



if win_shift_e == 0.25*fs_e
    save('features.mat', 'features')

elseif win_shift_e == 5*fs_e
    save('features_5.mat', 'features')

elseif win_shift_e == 30*fs_e
    save('features_30.mat', 'features')
end
