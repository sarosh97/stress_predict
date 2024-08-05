%This function extracts the features from the accelerometer data
%Input should be an nx3 matrix
%The output is 1x15 array
%The extracted features are, and in this order,
                   %mean (uA, where A=X,Y,Z)                    1x3
                   %absolute mean (abs_uA)                      1x1
                   %standard deviation (sA, where A=X,Y,Z)      1x3
                   %absolute standard deviation (abs_sA)        1x1
                   %absolute integral (intgA, where A=X,Y,Z)    1x3
                   %absolute integral overall (abs_intgA)       1x1
                   %peak frequency (pk_freqA, where A=X,Y,Z)    1x3
                       
function feat = feat_ext_acc(signal, fs)

%mean
feat(1,1:3) = mean(signal);             %1x3
feat(1,4) = sqrt(sum(feat(1,1:3).^2));  %1x1

%standard deviation
feat(1,5:7) = std(signal);              %1x3
feat(1,8) = sqrt(sum(feat(1,5:7).^2));  %1x1

%absolute integral
feat(1,9:11) = trapz(abs(signal));      %1x3
feat(1,12) = sqrt(sum(feat(1,9:11).^2));%1x1

%peak frequency
N = length(signal);     
freq = (0:round(N/2))*fs/N;

FT_X = abs(fft(signal(:,1))/N);    %normalized FFT
FT_X = FT_X(1:round(N/2));         %taking only one side of fft 
FT_Y = abs(fft(signal(:,2))/N);    %normalized FFT
FT_Y = FT_Y(1:round(N/2));         %taking only one side of fft 
FT_Z = abs(fft(signal(:,3))/N);    %normalized FFT
FT_Z = FT_Z(1:round(N/2));         %taking only one side of fft 

[~, idx_x] = max(FT_X);
[~, idx_y] = max(FT_Y);
[~, idx_z] = max(FT_Z);

feat(1,13:15) = [freq(idx_x), freq(idx_y), freq(idx_z)]; %1x3
end