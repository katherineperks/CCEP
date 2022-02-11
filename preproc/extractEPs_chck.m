function [fs,response_bypair,avgresponse_bypair,selected_chans] = extractEPs_chck(subject,fileName,window_time,state)
% returns all channel responses and the average channel response
    % to stimulation at each pair of channels (ie stim pair)
    % uses chck to exclude bad trials (artifacts) from the channel responses

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
    errID = 'MyFunction:IncompatibleDataset';
    msg = 'Unable to exclude trials before performing artifact removal.';
    incorrectState = MException(errID,msg);
    throw(incorrectState);
end

stimpairs = [anode cathode];
N_chan = size(usedata,2);
N_pair = length(unique(stimpairs,'rows'))/2;

onsets_bychan = cell(1,N_chan); % stim onsets for each channel
onsets_bypair = cell(1,N_pair); % stim onsets for each stim pair
art_bychan = cell(1,N_chan); % unremoved artifacts (exclude these trials)
for i=1:N_chan
    onsets_bychan{i} = onsets_samps(anode==i | cathode==i);
    art_bychan{i} = onsets_samps(chck(:,i)==1);
end
if N_pair > 7
    for i=1:7
        onsets_bypair{i} = onsets_samps((anode==i & cathode==i+1) | (anode==i+1 & cathode==i));
    end
    for i=9:N_chan-1
        onsets_bypair{i-1} = onsets_samps((anode==i & cathode==i+1) | (anode==i+1 & cathode==i));
    end
else
    for i=1:N_pair
        onsets_bypair{i} = onsets_samps((anode==i & cathode==i+1) | (anode==i+1 & cathode==i));
    end
end

[~,chan] = max(cellfun(@length,onsets_bypair)); % identify main stim pair
if chan > 7
    selected_chans = [chan+1 chan+2];
else
    selected_chans = [chan chan+1];
end

window_samps = fs*window_time;
window_samps = [floor(window_samps(1)), ceil(window_samps(2))];
window_sz = window_samps(2) - window_samps(1) + 1;

response_bypair = cell(N_pair,N_chan);
avgresponse_bypair = zeros(N_pair,N_chan,window_sz);

for P_s=1:N_pair % stim pair
    N_stims = length(onsets_bypair{P_s});
    for C_r=1:N_chan % response channel
        trial = 0;
        N_trials = N_stims - length(art_bychan{C_r});
        response_bypair{P_s,C_r} = zeros(N_trials,window_sz);
        for stim=1:N_stims % stim number delivered at stim pair
            if ~ismember(onsets_bypair{P_s}(stim),art_bychan{C_r}) % only include trials that passed artifact removal
                trial = trial + 1;
                window = onsets_bypair{P_s}(stim) + window_samps;
                start = window(1);
                stop = window(2);
                response_bypair{P_s,C_r}(trial,:) = transpose(usedata(start:stop,C_r));
            end
        end
        avgresponse_bypair(P_s,C_r,:) = mean(response_bypair{P_s,C_r},1);        
    end
end

% response_bychan = cell(N_chan);
% for C_s=1:N_chan % stim channel
%     N_stims = length(onsets_bychan{C_s});
%     for C_r=1:N_chan % response channel
%         response_bychan{C_s,C_r} = zeros(N_stims,window_sz);
%         for stim=1:N_stims % stim number delivered at stim channel
%             window = onsets_bychan{C_s}(stim) + window_samps;
%             start = window(1);
%             stop = window(2);
%             response_bychan{C_s,C_r}(stim,:) = transpose(usedata(start:stop,C_r));
%         end
%     end
% end