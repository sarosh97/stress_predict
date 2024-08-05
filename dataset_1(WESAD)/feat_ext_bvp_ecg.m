%This function extracts the features from ECG and BVP data
%Input should be an nx1 array 
%The output is a 1xm array where m = #different features
%The argument cond should be a string to deinfe which the input modality
%cond='ecg' or cond='bvp'

%The extracted features are, and in this order,
                           %mean of heart rate (uHR)
                           %standard deviation of heart rate (sHR)
                           %rms of heart rate (rHR)
                           %standard deviation of NN intervals (sdnn)
                           %rms of successive differences between beats (rmssd)
                           %NN intervals differing more than 50ms (nn50)
                           %percentage of NN intervals differing more than 50ms (pnn50)
                           %energy in different frequency bands
                                                %ultra-low frequency(e_ulf)
                                                %low-frequency (e_lf)
                                                %high-frequency (e_hf)
                                                %ultra-high frequency (e_uhf)
                           %ratio of lf and hf = lf/hf (lf_hf)
                           %relative power of each frequency band
                                                %r_ulf
                                                %r_lf
                                                %r_hf
                                                %r_uhf
                           %norm of LF (n_lf)
                           %norm of HF (n_HF)
                           
                       
function feat = feat_ext_bvp_ecg(signal, fs, cond)

if cond == 'bvp'
    %Heart Rate related features from wrist BVP
    [pks, locs] = findpeaks(-signal,'MinPeakDistance', 0.5*fs);    %systolic feet and their location
    th = prctile(signal(locs), 90);     %threshold at 90th percentile to remove spurious peaks
    locs(pks<th) = [];
    
    ibi_w = diff(locs)/fs;   %Inter-Beat Interval

    feat(1,1) = mean(60./ibi_w);  %mean heart rate bpm
    feat(1,2) = std(60./ibi_w);   %standard deviation of heart rate
    feat(1,3) = rms(60./ibi_w);   %rms of heart rate
    feat(1,4) = std(ibi_w)*1000; %standard deviation of IBI (NN) intervals in ms
    feat(1,5) = sqrt(mean(diff(ibi_w).^2))*1000; %in ms
    feat(1,6) = sum(abs(diff(ibi_w)*1000) > 50);   %NN intervals differing more than 50ms
    feat(1,7) = 100*(feat(1,6))/length(ibi_w);   %percentage of NN intervals differing more than 50ms
    
    %Energy in different frequency bands
    N = length(signal);     
    freq = (0:round(N/2))*fs/N;
    FT_w = abs(fft(signal)/N);    %normalized FFT
    FT_w = FT_w(1:round(N/2));    %taking only one side of fft
    
    ulf = [0.01 0.04];   %ultra-low frequency band
    lf = [0.04 0.15];    %low frequency band
    hf = [0.15 0.4];     %high frequency band
    uhf = [0.4 1];       %ultra-high frequency band

    ulf_idx = find(freq >= ulf(1) & freq <= ulf(2));
    lf_idx = find(freq >= lf(1) & freq <= lf(2));
    hf_idx = find(freq >= hf(1) & freq <= hf(2));
    uhf_idx = find(freq >= uhf(1) & freq <= uhf(2));

    feat(1,8) = sum(FT_w(ulf_idx).^2);        %energy in ULF
    feat(1,9) = sum(FT_w(lf_idx).^2);        %energy in LF
    feat(1,10) = sum(FT_w(hf_idx).^2);        %energy in HF
    feat(1,11) = sum(FT_w(uhf_idx).^2);       %energy in UHF
    
    feat(1,12) = (feat(1,10))/(feat(1,11));   %LF/HF
    
    total_energy = sum(FT_w.^2);
    feat(1,13) = (feat(1,9))/total_energy;    %relative power ulf
    feat(1,14) = (feat(1,10))/total_energy;   %relative power lf
    feat(1,15) = (feat(1,11))/total_energy;   %relative power hf
    feat(1,16) = (feat(1,12))/total_energy;   %relative power uhf
    
    feat(1,17) = norm(FT_w(lf_idx));          %norm of LF
    feat(1,18) = norm(FT_w(hf_idx));          %norm of HF
    
    
    
elseif cond == 'ecg'
    %Heart Rate related features from chest ECG
    [pks, locs] = findpeaks(signal,'MinPeakDistance', 0.5*fs);    %R peaks and their location
    rr_int = diff(locs)/fs;   %R-R Interval
    
    feat(1,1) = mean(60./rr_int);  %mean heart rate bpm
    feat(1,2) = std(60./rr_int);   %standard deviation of heart rate
    feat(1,3) = rms(60./rr_int);   %rms of heart rate
    feat(1,4) = std(rr_int)*1000; %standard deviation of IBI intervals in ms
    feat(1,5) = sqrt(mean(diff(rr_int).^2))*1000; %in ms
    feat(1,6) = sum(abs(diff(rr_int)*1000) > 50);   %NN intervals differing more than 50ms
    feat(1,7) = 100*(feat(1,6))/length(rr_int);   %percentage of NN intervals differing more than 50ms
    
    %Energy in different frequency bands
    N = length(signal);     
    freq = (0:round(N/2))*fs/N;
    FT_c = abs(fft(signal)/N);    %normalized FFT
    FT_c = FT_c(1:round(N/2));    %taking only one side of fft
    
    ulf = [0.01 0.04];   %ultra-low frequency band
    lf = [0.04 0.15];    %low frequency band
    hf = [0.15 0.4];     %high frequency band
    uhf = [0.4 1];       %ultra-high frequency band

    ulf_idx = find(freq >= ulf(1) & freq <= ulf(2));
    lf_idx = find(freq >= lf(1) & freq <= lf(2));
    hf_idx = find(freq >= hf(1) & freq <= hf(2));
    uhf_idx = find(freq >= uhf(1) & freq <= uhf(2));

    feat(1,8) = sum(FT_c(ulf_idx).^2);       %energy in ULF
    feat(1,9) = sum(FT_c(lf_idx).^2);       %energy in LF
    feat(1,10) = sum(FT_c(hf_idx).^2);       %energy in HF
    feat(1,11) = sum(FT_c(uhf_idx).^2);      %energy in UHF
    
    feat(1,12) = (feat(1,10))/(feat(1,11));  %LF/HF
    
    total_energy = sum(FT_c.^2);
    feat(1,13) = (feat(1,9))/total_energy;   %relative power ulf
    feat(1,14) = (feat(1,10))/total_energy;  %relative power lf
    feat(1,15) = (feat(1,11))/total_energy;  %relative power hf
    feat(1,16) = (feat(1,12))/total_energy;  %relative power uhf
    
    feat(1,17) = norm(FT_c(lf_idx));          %norm of LF
    feat(1,18) = norm(FT_c(hf_idx));          %norm of HF
    
end
end