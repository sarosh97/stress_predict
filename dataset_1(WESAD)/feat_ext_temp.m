%This function extracts the features from temperature data
%PInput should be an nx1 array 
%The output is a 1xm array where m = #different features

%The extracted features are, and in this order,
                           %mean (uHR)
                           %standard deviation (sHR)
                           %minimum value
                           %maximum value
                           %dynamic range
                           %slope
                           
function feat = feat_ext_temp(signal, fs)

feat(1,1) = mean(signal);
feat(1,2) = std(signal);
feat(1,3) = min(signal);
feat(1,4) = max(signal);
feat(1,5) = max(signal) - min(signal);
%feat(1,6) = diff(signal)*fs;

end