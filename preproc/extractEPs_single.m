function [fs,responses,avgresponses,selected_chans] = extractEPs_single(subject,fileName,window_time,state)
% returns all channel responses and the average channel response
    % to stimulation at ONE pair of channels (ie selected stim pair)

if isempty(window_time)
    window_time = [-0.5,1];
end

cd("C:\Users\Katherine\Documents\AAA UW\Herron Lab\Data\" + subject)
load(fileName);
if state=="filt"
    usedata = data_filt;
elseif state=="filt2"
    usedata = data_filt2;
elseif state=="artrem"
    usedata = data_artrem;
elseif state=="raw"
    usedata = data;
end

stimpairs = [anode cathode];
N_trials = length(stimpairs);
N_chan = size(usedata,2);
N_pair = length(unique(stimpairs,'rows'))/2;

if N_pair == 1
    selected_chans = [anode(1) cathode(1)];
    N_stims = length(stimpairs);
else
    errID = 'MyFunction:IncompatibleDataset';
    msg = 'Unable to perform this extraction on a dataset with >1 stim pair.';
    incorrectState = MException(errID,msg);
    throw(incorrectState);
end

window_samps = fs*window_time;
window_samps = [floor(window_samps(1)), ceil(window_samps(2))];
window_sz = window_samps(2) - window_samps(1) + 1;

responses = zeros(N_chan,N_stims,window_sz);
avgresponses = zeros(N_chan,window_sz);

for C_r=1:N_chan % response channel
    for stim=1:N_stims % stim number delivered at stim pair
        window = onsets_samps(stim) + window_samps;
        start = window(1);
        stop = window(2);
        responses(C_r,stim,:) = transpose(usedata(start:stop,C_r));
    end
    avgresponses(C_r,:) = mean(responses(C_r,:,:),2);
end