function [data_artrem,data_filt] = filterEPs(subject,session)
% filters artifact-removed data with highpass, lowpass, and bandpass filters
    % saves filtered data to subject file

fileName ="iTBS_"+subject+"_"+session+".mat"; % example: subject = "92a04"; session = "pre";
cd("C:\Users\Katherine\Documents\AAA UW\Herron Lab\Data\" + subject)
load(fileName)

order = 3;
freqs = {.5 200 [53,63] [117,123] [177,183]};
ftype = ["high" "low" "stop" "stop" "stop"];
data_filt = double(data_artrem);
for i=1:length(freqs)
    freq = freqs{i};
    [b,a] = butter(order,freq/(fs/2),ftype(i));
%     freqz(b, a, 0:.1:400, fs);
    data_filt = filtfilt(b,a,data_filt);
end

save(fileName,'data_filt','-append')
end