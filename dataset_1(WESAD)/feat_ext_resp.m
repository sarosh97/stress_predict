%This function extracts the features from respiration data
%Input should be an nx1 array

%The extracted features are, and in this order,
                           %mean of inhalation duration
                           %standard deviation of inhalation duration
                           %mean of exhalation duration
                           %standard deviation of exhalation duration
                           %mean of I/E ratio
                           %standard deviation of I/E ratio
                           %mean of stretch
                           %std of stretch
                           %minute volume
                           %breath rate
                           %mean of respiration duration in one minute
                           %std of respiration duration in one minute
                           
function feat = feat_ext_resp(signal, t, fs)

[pks_up, locs_up] = findpeaks(signal,'MinPeakDistance', 1.5*fs);    %threshold from Kurt Plare et al 2011

th = prctile(signal(locs_up), 25);  %to remove spurious peaks
locs_up(pks_up<th) = [];
pks_up(pks_up<th) = [];

%finding exhalation valley which is the minimum value between consecutive inhalation peaks
for i = 1:length(locs_up)
    if i == length(locs_up)  %if last value
        [pks_down(i), temp] = min(signal(locs_up(i):end));
        locs_down(i,1) = locs_up(i) + temp - 1;
    else
        [pks_down(i), temp] = min(signal(locs_up(i):locs_up(i+1)));
        locs_down(i,1) = locs_up(i) + temp - 1;
    end
end

exh_dur = (locs_down - locs_up)/fs;                 %exhalation duration
inh_dur = (locs_up - [1; locs_down(1:end-1)])/fs;   %inhalation duration

feat(1,1) = mean(inh_dur);              %mean of inhalation duration
feat(1,2) = std(inh_dur);               %std of inhalation duration    
feat(1,3) = mean(exh_dur);              %mean of exhalation duration
feat(1,4) = std(exh_dur);               %std of exhalation duration
feat(1,5) = mean(inh_dur./exh_dur);     %mean of I/E ratio
feat(1,6) = std(inh_dur./exh_dur);      %standard deviation of I/E ratio
feat(1,7) = mean(signal(locs_up) - signal(locs_down));  %mean of stretch
feat(1,8) = std(signal(locs_up) - signal(locs_down));   %std of stretch

%computing inhalation volume
i = 1;  
for j = 1:length(inh_dur)
    inh_vol(j) = trapz(t(i:locs_up(j)), signal(i:locs_up(j)));  %(i-1) valley to (i) peak -> testing can be done using "plot(t, signal, t(locs_down), signal(locs_down), 'k*', t(locs_up), signal(locs_up), 'r*')"
    i = locs_down(j);
end

feat(1,9) = sum(inh_vol);       %minute volume (inhalation volume in one minute)
feat(1,10) = length(locs_up)/60; %breath rate
feat(1,11) = mean(inh_dur + exh_dur);    %mean of respiration duration in one minute
feat(1,12) = std(inh_dur + exh_dur);    %std of respiration duration in one minute

end