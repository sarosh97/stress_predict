%Note that this function uses a function "cvxEDA". It can be downloaded
%from https://www.mathworks.com/matlabcentral/fileexchange/53326-cvxeda

%This function extracts the features from wrist EDA data
%Input should be an nx1 array 
%The output is a 1xm array where m = #different features

%The extracted features are, and in this order,
                           %mean of EDA
                           %standard deviation of EDA
                           %minimum amplitude of EDA
                           %maximum amplitude of EDA
                           %skewness of EDA 
                           %kurtosis of EDA
                           %dyanmic range of EDA
                           %mean of phasic component
                           %standard deviation of phasic component
                           %skewness of phasic component
                           %kurtois of phasic component
                           %minimum of phasic component
                           %maximum of phasic component
                           %range of phasic component
                           %mean of 1st derivative of phasic component
                           %STD of 1st derivative of phasic component
                           %mean of 2nd derivative of phasic component
                           %STD of 2nd derivative of phasic component
                           %mean of tonic component
                           %standard deviation of tonic component
                           %Number of SCR segments
                           %Correlation of SCR and time
                           %Sum of startle magnitudes
                           %mean of startle magnitudes
                           %standard deviation of startle magnitudes
                           %sum of response durations
                           %mean of response durations
                           %standard deviation of response durations
                           %sum of AUC of startle segments
                           %Slope of EDA (mean value)
                           %Max frequency of EDA
                           %max(pks)/mean(SCL)
                           %sum, mean, std, skew, kurt, max, min, range of PSD of SCR
                           
function feat = feat_ext_eda(signal, t, fs)

[SCR, ~, SCL] = cvxEDA(signal, 1/fs);
%SCR startle magnitudes
[pks, locs] = findpeaks(SCR, 'MinPeakDistance', 1.5*fs);

if length(pks) >= 2
   
    feat(1,1) = mean(SCR);  %mean of phasic component
    feat(1,2) = std(SCR);   %standard deviation of phasic component
    feat(1,3) = skewness(SCR);  %skewness of phasic component
    feat(1,4) = kurtosis(SCR);  %kurtosis of phasic component
    feat(1,5) = min(SCR);  %minimum of phasic component
    feat(1,6) = max(SCR);  %maximum of phasic component
    feat(1,7) = max(SCR) - min(SCR);  %range of phasic component
    feat(1,8) = mean(diff(SCR)*fs);   %mean of 1st derivative of SCR
    feat(1,9) = std(diff(SCR)*fs);    %std of 1st derivative of SCR
    feat(1,10) = mean(diff(diff(SCR)*fs)*fs);    %mean of 2nd derivative of SCR
    feat(1,11) = std(diff(diff(SCR)*fs)*fs);    %STD of 2ND derivative of SCR

    feat(1,11) = mean(SCL);  %mean of tonic component
    feat(1,12) = std(SCL);  %standard deviation of tonic component

    %peaks are found but the local minima between consecutive peaks need to be 
    %computed to calculate response durations
    d = [0;diff(SCR)*fs];
    for i = 1:length(locs)
        if i == length(locs)
            [pk, temp] = min(d(locs(i):end));
            idx2(i) = locs(i) + temp - 1;
        else
            [pk2, temp2] = min(SCR(locs(i):locs(i)+1.5*fs));
            post_idx = locs(i+1);
            idx2(i) = locs(i) + temp2 - 1;
        end
        
        if i==1 && (locs(i)-1.5*fs)<1 
            pk1 = SCR(1);
            idx1(1) = 1;
        else
            [pk1, temp1] = min(SCR(locs(i)-1.5*fs:locs(i)));
            pre_idx = locs(i);
            idx1(i) = pre_idx + temp1 - 1.5*fs - 1;
        end
        
        pre_duration(i) = t(locs(i)) - t(idx1(i));
        post_duration(i) = t(idx2(i)) - t(locs(i)+1);
        total_duration(i) = pre_duration(i) + post_duration(i);
    
        startle_mag(i) = pks(i) - pk1;
        startle_AUC(i) = trapz(t(idx1(i):locs(i)) , SCR(idx1(i):locs(i)));
    
    end

    temp = corrcoef(SCL, t);
    feat(1,13) = temp(2,1);   %correlation coefficient between SCL and t
    feat(1,14) = length(pks);        %Number of SCR segments
    feat(1,15) = sum(startle_mag);   %Sum of startle magnitudes
    feat(1,16) = mean(startle_mag);  %mean of startle magnitudes
    feat(1,17) = std(startle_mag);   %standard deviation of startle magnitudes
    feat(1,18) = sum(pre_duration);  %sum of response durations
    feat(1,19) = mean(pre_duration); %mean of response durations
    feat(1,20) = std(pre_duration);  %standard deviation of response durations
    feat(1,21) = sum(startle_AUC);   %sum of AUC of startle segments
    feat(1,22) = mean(diff(signal)*fs);

    %frequency based features
    N = length(signal);     
    freq = (0:round(N/2))*fs/N;

    FT = abs(fft(signal)/N);    %normalized FFT 
    FT = FT(1:round(N/2));         %taking only one side of fft 

    [~, idX] = max(FT);

    feat(1,23) = freq(idX);
    feat(1,24) = max(pks)/mean(SCL);

    [Pxx,F] = pwelch(SCR,[],[],[],fs);
    feat(1,25) = trapz(F,Pxx);
    feat(1,26) = skewness(Pxx);
    feat(1,27) = kurtosis(Pxx);
    feat(1,28) = max(Pxx);
    
else
    feat(1,1:28) = nan;
end

end