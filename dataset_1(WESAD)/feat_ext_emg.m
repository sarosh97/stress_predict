%This function extracts the features from EMG data
%Input should be an nx1 array 
%The output is a 1xm array where m = #different features


%the extracted features are, and in this order,
                           %number of peaks
                           %mean of peaks' amplitudes
                           %standard deviation of peaks' amplitudes
                           %sum of peaks' amplitudes
                           %normalized sum of peaks' amplitudes
                           
                           
function feat = feat_ext_emg(signal, fs)
    
    [pks, locs] = findpeaks(signal,'MinPeakDistance', 0.02*fs);    %Peaks should be atleast 0.02s apart
    feat(1,1) = length(pks);           %number of peaks
    feat(1,2) = mean(pks);             %mean of peaks' amplitudes
    feat(1,3) = std(pks);              %standard deviation of peaks
    feat(1,4) = sum(pks);              %sum of peaks
    feat(1,5) = sum(pks)/length(pks);  %normalized sum of peaks
end