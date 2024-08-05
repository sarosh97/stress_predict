
clc
clear all
close all

%% Loading Data using custom function 

%Data is a 1xj cell array, where j = #subjects
%labels_num, start_time, end_time are 15x5 matrices
%Labels are numbered as baseline=0, amusement(fun)=1, meditation=2, 
%stress(tsst)=3 

%Filename should be the location of the folder containing data
filename = 'D:\Research Papers and Data\Wearable sensing\WESAD (dataset_1)';
[data, labels_num, start_time, end_time] = load_wesad_data(filename);

%Window Shift size
shift_size = 30;

%% Separating data modality-wise and converting raw data to SI units

%Each variable is 1xj cell where j = #subjects
%Each cell has two arrays, one for concerned data and other for time in sec

%Defining samping Frequencies
c_fs = 700;   %sampling frequency of respiban (attached on chest)
w_a_fs = 32;  %sampling rate of wrist accelerometer
w_b_fs = 64;  %sampling rate of wrist PPG sensor
w_e_fs = 4;   %sampling rate of wrist EDA sensor
w_t_fs = 4;   %sampling rate of wrist temperature sensor

%Defining Low-Pass filter for Chest Respiration Data
fc1 = [0.35];     %Low-pass at 0.35 Hz
fc2 = [1];        %0.65 Hz of transition band
Rp = 3;
Rs = 40;
Wp = fc1/(c_fs/2);
Ws = fc2/(c_fs/2);
[N, Wp] = ellipord(Wp, Ws, Rp, Rs);
[b_resp, a_resp] = ellip(N, Rp, Rs, Wp);

%Defining Low-Pass filter for Chest EMG Data
fc1 = [20];     %Low-pass filter at 50 Hz
fc2 = [25];     %5 Hz of transition band
Rp = 3;
Rs = 40;
Wp = fc1/(c_fs/2);
Ws = fc2/(c_fs/2);
[N, Wp] = ellipord(Wp, Ws, Rp, Rs);
[b_emg, a_emg] = ellip(N, Rp, Rs, Wp);

vcc = 3;
chan_bit = 2^16;

for i = 1:length(data)
    
    %time vector for chest based modalities
    %can use any modality to define time vector since chest-modalities have same fs
    tc = [0:length(data{i}.chest.ACC)-1]'/c_fs; %(column vector)
    
    %chest accelerometer n*3 matrix, unit = g
    cmin = 28000;
    cmax = 38000;
    chest_accX{i} = [detrend((data{i}.chest.ACC(:,1) - cmin)./(cmax - cmin)*2-1) , tc];
    chest_accY{i} = [detrend((data{i}.chest.ACC(:,2) - cmin)./(cmax - cmin)*2-1) , tc];
    chest_accZ{i} = [detrend((data{i}.chest.ACC(:,3) - cmin)./(cmax - cmin)*2-1) , tc];
    
    %chest ecg, unit = mv
    chest_ecg{i} = [detrend(((data{i}.chest.ECG /chan_bit-0.5)*vcc)) , tc];
    
    %chest emg, unit = mv
    data_filt_emg = filtfilt(b_emg,a_emg,((data{i}.chest.EMG /chan_bit-0.5)*vcc));
    chest_emg_60{i} = [detrend(data_filt_emg) , tc];   %to be processed as 60s window segments
    
    %chest temp, unit = C
    ntcv = (data{i}.chest.Temp*vcc)/(chan_bit-1);
    ntcr = ((1e-4)*ntcv)./(vcc - ntcv); 
    a0 = 1.12764514e-3;
    a1 = 2.34282709e-4;
    a2 = 8.77303013e-8;
    tmpK = 1./(a0 + a1*log(ntcr) + a2*[log(ntcr)].^3);
    chest_temp{i} = [(tmpK - 273.15) , tc];
    
    %chest respiration, unit = %
    filt_resp = filtfilt(b_resp,a_resp, (data{i}.chest.Resp/chan_bit - 0.5)*100);
    chest_resp{i} = [detrend(filt_resp) , tc];
    
    
    %time vector for wrist based modalities (column vector)
    ta = [0:length(data{i}.wrist.ACC)-1]'/w_a_fs;
    tb = [0:length(data{i}.wrist.BVP)-1]'/w_b_fs;
    te = [0:length(data{i}.wrist.EDA)-1]'/w_e_fs;
    tt = [0:length(data{i}.wrist.TEMP)-1]'/w_t_fs;
    
    %wrist accelerometer n*3 matrix, unit = 1/64 g
    wrist_accX{i} = [detrend(data{i}.wrist.ACC(:,1)) , ta];
    wrist_accY{i} = [detrend(data{i}.wrist.ACC(:,2)) , ta];
    wrist_accZ{i} = [detrend(data{i}.wrist.ACC(:,3)) , ta];
    
    %wrist BVP (from PPG), unit = ?
    wrist_bvp{i} = [detrend(data{i}.wrist.BVP) , tb];

    %wrist EDA , unit = uS
    wrist_eda{i} = [(data{i}.wrist.EDA) , te];
    
    %wrist Temp , unit = C
    wrist_temp{i} = [data{i}.wrist.TEMP , tt];
    
    i    %keep track of progress
end

clear data    %to save memory

%% Feature extraction 1: Chest accelerometer

%extracting chest accelerometer features
win_shift = shift_size*c_fs;    %window shift of 0.25 seconds
win_len_acc = 5*c_fs;    %window length of 5 seconds for accelerometer signals

subj = 15;
for i = 1:subj                  %loop through subjects
    k = 1;
    tc = chest_accX{i}(:,2);    %time vector for each subject. Same for within subject variables.
    
    for j = 1:size(start_time,2)        %loop through testing conditions
        subj_name = ['S', num2str(i)];  %to create structure name based on subject number
        
        [~, st_idx] = min(abs(tc-start_time(i,j)));
        [~, end_idx] = min(abs(tc-end_time(i,j)));
        current_idx = st_idx;
        
        while (current_idx < end_idx - win_len_acc + 1)        %sliding windows and computing features for each condition
            current_win = [current_idx:current_idx + win_len_acc - 1];  %indices for the current window
            signal_acc = [chest_accX{i}(current_win,1), chest_accY{i}(current_win,1), chest_accZ{i}(current_win,1)];
            
            features_acc = feat_ext_acc(signal_acc , c_fs); %accelerometer features
            
            %appending each value row-wise to each variable
            feat.(subj_name).chest.acc.feat(k,:) = features_acc;
            
            %storing labels alongside features
            feat.(subj_name).chest.acc.labels(k,1) = labels_num(i,j);
            
            current_idx = current_idx + win_shift;
            k = k+1;
        end
    end
    disp(['Chest ACC features of S', num2str(i), ' computed'])
end


%% Feature extraction 2: Wrist accelerometer

clear features  %the varibale to be used again

%extracting wrist accelerometer features
win_shift = shift_size*w_a_fs;    %window shift of 0.25 seconds
win_len_acc = 5*w_a_fs;    %window length of 5 seconds for accelerometer signals

subj = 15;
for i = 1:subj                  %loop through subjects
    k = 1;
    ta = wrist_accX{i}(:,2);    %time vector of each subject. Same for within subject variables.
    
    for j = 1:size(start_time,2)        %loop through testing conditions
        subj_name = ['S', num2str(i)];  %to create structure name based on subject number
        
        [~, st_idx] = min(abs(ta-start_time(i,j))); %to find the index of start_time in the time vector
        [~, end_idx] = min(abs(ta-end_time(i,j)));  %to find the index of end_time in the time vector
        current_idx = st_idx;
        
        while (current_idx < end_idx - win_len_acc + 1)        %sliding windows and computing features for each condition
            current_win = [current_idx:current_idx + win_len_acc - 1];  %indices for the current window
            signal = [wrist_accX{i}(current_win,1), wrist_accY{i}(current_win,1), wrist_accZ{i}(current_win,1)];
            features = feat_ext_acc(signal , w_a_fs);
            
            %appending each value row-wise to each variable
            feat.(subj_name).wrist.acc.feat(k,:) = features;
            
            %storing labels alongside features
            feat.(subj_name).wrist.acc.labels(k,1) = labels_num(i,j);
            
            current_idx = current_idx + win_shift;
            k = k+1;
        end
    end
    disp(['Wrist ACC features of S', num2str(i), ' computed'])  %keep track of progress
end

%% Feature extraction 3: Chest ECG-EMG60-Resp-Temp and Wrist BVP-TEMP-EDA 

clear features  %the varibale to be used again

%extracting chest-based modalities' features
win_shift_c = shift_size*c_fs;    %window shift of 0.25 seconds (chest)
win_len_phy_c = 60*c_fs;    %window length of 60 seconds for physiological signals (chest)

win_shift_wb = shift_size*w_b_fs;    %window shift of 0.25 seconds (wrist_bvp)
win_len_phy_wb = 60*w_b_fs;    %window length of 60 seconds for physiological signals (wrist_bvp)

win_shift_wt = shift_size*w_t_fs;    %window shift of 0.25 seconds (wrist_temp)
win_len_phy_wt = 60*w_t_fs;    %window length of 60 seconds for physiological signals (wrist_temp)

win_shift_we = shift_size*w_e_fs;    %window shift of 0.25 seconds (wrist_eda)
win_len_phy_we = 60*w_e_fs;    %window length of 60 seconds for physiological signals (wrist_eda)

subj = 15;
for i = 1:subj                  %loop through subjects
    K = 1;  %inner-most loop variable, wrsit_bvp
    L = 1;  %inner-most loop variable, wrist_temp
    M = 1;  %inner-most loop variable, wrist_eda
    N = 1;  %inner-most loop variable, chest ECG, EMG60, Temp, Resp
    
    tw_b = wrist_bvp{i}(:,2);  %time vector of each subject for wrist_bvp
    tw_t = wrist_temp{i}(:,2); %time vector of each subject for wrist_temp
    tw_e = wrist_eda{i}(:,2); %time vector of each subject for wrist_eda
    tc = chest_ecg{i}(:,2);    %time vector of each subject. Same for within subject variables.
    
    for j = 1:size(start_time,2)        %loop through testing conditions
        subj_name = ['S', num2str(i)];  %to create structure name based on subject number
        
        %Features from Wrist_BVP
        [~, st_idx] = min(abs(tw_b - start_time(i,j))); %to find the index of start_time in the time vector
        [~, end_idx] = min(abs(tw_b - end_time(i,j)));  %to find the index of end_time in the time vector
        current_idx = st_idx;
        
        while (current_idx < end_idx - win_len_phy_wb + 1)        %sliding windows and computing features for each condition
            current_win = [current_idx:current_idx + win_len_phy_wb - 1];  %indices for the current window
            signal = wrist_bvp{i}(current_win,1);
            features_bvp = feat_ext_bvp_ecg(signal , w_b_fs, 'bvp');
            
            %appending each value row-wise to each variable
            feat.(subj_name).wrist.bvp.feat(K,:) = features_bvp;     
            
            %storing labels alongside features
            feat.(subj_name).wrist.bvp.labels(K,1) = labels_num(i,j);
            
            current_idx = current_idx + win_shift_wb;
            K = K+1;
        end
        
        
        %Features from Wrist_Temp
        [~, st_idx] = min(abs(tw_t - start_time(i,j))); %to find the index of start_time in the time vector
        [~, end_idx] = min(abs(tw_t - end_time(i,j)));  %to find the index of end_time in the time vector
        current_idx = st_idx;
        
        while (current_idx < end_idx - win_len_phy_wt + 1)        %sliding windows and computing features for each condition
            current_win = [current_idx:current_idx + win_len_phy_wt - 1];  %indices for the current window
            signal = wrist_temp{i}(current_win,1);
            features_tempw = feat_ext_temp(signal, w_t_fs);
            
            %appending each value row-wise to each variable
            feat.(subj_name).wrist.temp.feat(L,:) = features_tempw;     
            
            %storing labels alongside features
            feat.(subj_name).wrist.temp.labels(L,1) = labels_num(i,j);
            
            current_idx = current_idx + win_shift_wt;
            L = L+1;
        end
        
        
        %Features from Wrist_EDA
        [~, st_idx] = min(abs(tw_e - start_time(i,j))); %to find the index of start_time in the time vector
        [~, end_idx] = min(abs(tw_e - end_time(i,j)));  %to find the index of end_time in the time vector
        current_idx = st_idx;
        
        while (current_idx < end_idx - win_len_phy_we + 1)        %sliding windows and computing features for each condition
            current_win = [current_idx:current_idx + win_len_phy_we - 1];  %indices for the current window
            time_win_eda = tw_e(current_win,1);
            signal = wrist_eda{i}(current_win,1);
            features_eda = feat_ext_eda(signal, time_win_eda, w_e_fs);
            
            %appending each value row-wise to each variable
            feat.(subj_name).wrist.eda.feat(M,:) = features_eda;     
            
            %storing labels alongside features
            feat.(subj_name).wrist.eda.labels(M,1) = labels_num(i,j);
            
            current_idx = current_idx + win_shift_we;
            M = M+1;
        end
        
        
        %Features from Chest ECG, Respiration, Temperature
        [~, st_idx] = min(abs(tc-start_time(i,j))); %to find the index of start_time in the time vector
        [~, end_idx] = min(abs(tc-end_time(i,j)));  %to find the index of end_time in the time vector
        current_idx = st_idx;
        
        while (current_idx < end_idx - win_len_phy_c + 1)        %sliding windows and computing features for each condition
            current_win = [current_idx:current_idx + win_len_phy_c - 1];  %indices for the current window
            signal_e = chest_ecg{i}(current_win,1);       %ECG
            signal_emg_60 = chest_emg_60{i}(current_win,1);  %EMG
            signal_r = chest_resp{i}(current_win,1);      %Respiration
            signal_t = chest_temp{i}(current_win,1);      %Temperature
            
            time_win_resp = tc(current_win,1);
            
            feat_ecg = feat_ext_bvp_ecg(signal_e , c_fs, 'ecg');
            feat_emg60 = feat_ext_emg(signal_emg_60 , c_fs);
            feat_resp = feat_ext_resp(signal_r , time_win_resp, c_fs);
            feat_temp = feat_ext_temp(signal_t, c_fs);
            
            %appending each value row-wise to each variable
            feat.(subj_name).chest.ecg.feat(N,:) = feat_ecg;
            feat.(subj_name).chest.emg60.feat(N,:) = feat_emg60;
            feat.(subj_name).chest.resp.feat(N,:) = feat_resp;
            feat.(subj_name).chest.temp.feat(N,:) = feat_temp;
            
            %storing labels alongside features
            feat.(subj_name).chest.ecg.labels(N,1) = labels_num(i,j);
            feat.(subj_name).chest.emg60.labels(N,1) = labels_num(i,j);
            feat.(subj_name).chest.resp.labels(N,1) = labels_num(i,j);
            feat.(subj_name).chest.temp.labels(N,1) = labels_num(i,j);
            
            current_idx = current_idx + win_shift_c;
            N = N+1;
        end
        
    end
    disp(['Wrist BVP-Temp-EDA and Chest ECG-EMG60-Resp-Temp features of S', num2str(i), ' computed'])   %keep track of progress
end


if shift_size == 0.25
    save('features.mat', 'feat')

elseif shift_size == 5
    save('features_5.mat', 'feat')

elseif shift_size == 30
    save('features_30.mat', 'feat')
end

