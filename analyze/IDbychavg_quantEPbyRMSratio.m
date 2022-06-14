% this method does NOT perform baseline subtraction,
% identifies EP components using the zero-crossing indices of the channel average, 
% and quantifies each component using the ratio between RMSresp and RMSbase
% as well as the ratio between RMSresp overall and RMSbase
clear
subject = "0d5e8e";
fileName = "EP_Measure-4.mat";
cd("C:\Users\Katherine\Documents\AAA UW\Herron Lab\Data\" + subject)
load(fileName)
saveName = "IDbychavg_quantEPbyRMSratio";
if ~exist(saveName, 'dir')
   mkdir(saveName)
end
cd(saveName)

stimpairs = [anode cathode];
N_trials = length(stimpairs);
N_chan = size(data_artrem,2);
N_pair = length(unique(stimpairs,'rows'))/2;

if N_pair == 1
    selected_chans = [anode(1) cathode(1)]; % stim pair
    N_stims = length(stimpairs);
end

% extract EPs from each channel in filtered data
usedata = data_filt;
window_time = [0,1.5];
window_samps = fs*window_time;
window_samps = [floor(window_samps(1)), ceil(window_samps(2))];
window_sz = window_samps(2) - window_samps(1) + 1;

N_chansub = size(usedata,2); % look at a subset of all channels
responses = zeros(N_chansub,N_stims,window_sz);
avgresponses = zeros(N_chansub,window_sz);
for C_idx=1:N_chansub % response channel index (out of 20)
    C_r = C_idx+17-1; % response channel number (out of 126)
    for stim=1:N_stims % stim number delivered at stim pair
        window = onsets_samps(stim) + window_samps;
        start = window(1);
        stop = window(2);
        responses(C_idx,stim,:) = transpose(usedata(start:stop,C_idx));
    end
    avgresponses(C_idx,:) = mean(responses(C_idx,:,:),2);
end

% flags
showplot = false;
saveplot = false;
if fileName == "EP_Measure-3.mat"
    suffix = "_stim23_24";
elseif fileName == "EP_Measure-4.mat"
    suffix = "_stim27_28";
end
% function to find zero-crossing indices
zci = @(v) find(v.*circshift(v,1) <= 0);
N_components = 4;
% quantification metrics
peak_amps = zeros(N_chansub,N_stims,N_components);
peak_latens = zeros(N_chansub,N_stims,N_components);
RMS = zeros(N_chansub,N_stims,N_components); % ratio between RMSresp of each component and RMSbase
RMSoverall = zeros(N_chansub,N_stims); % ratio between RMSresp overall and RMSbase
RMSbase = zeros(N_chansub,N_stims);
RMSresp = zeros(N_chansub,N_stims,N_components);

for C_idx=1:N_chansub % response channel index (out of 20)
    C_r = C_idx+17-1; % response channel number (out of 126)
    % identify EP components from zci of channel average
    idx = zci(avgresponses(C_idx,:));
    art_idx = idx(idx<0.01*fs); % potential stim artifact
    real_idx = [floor(0.01*fs), idx(idx>=0.015*fs)];
    real_idx = [real_idx(diff(real_idx)>=0.015*fs), real_idx(end)]; % reliable signal
    
    % analyze channels with non-zero responses
    if length(idx) ~= length(responses)
        if showplot
            f = figure();
        else
            f = figure('Visible','off');
        end
        plot(1:window_sz,avgresponses(C_idx,:));
        hold on
        plot(real_idx,zeros(length(real_idx)),'ro')
        xline(0.01*fs,'r') % 10ms
        xline(0.5*fs,'r') % 500ms
        yline(0,'r')
        if saveplot
            figName = "components_ch"+num2str(C_r)+suffix;
            saveas(gcf,figName,'png')
        end
        
        % apply component windows identified in channel average to all trials
        for stim=1:N_stims                 
            % compute RMSbase from 200ms period before stim onset (within channel)
            window = onsets_samps(stim) + window_samps;
            start = window(1);
            stop = window(2);
            baseline = usedata(start-floor(0.2*fs):start-1,C_idx);
            RMSbase(C_idx,stim) = rms(baseline);
            % compute RMSoverall (ratio between RMS over entire response window and RMSbase)
            respwindow = responses(C_idx,stim,floor(0.01*fs):end); % exclude the first 10ms of response from quantification
            RMSoverall(C_idx,stim) = rms(respwindow)/rms(baseline);
            % extract single-trial peak amplitude, peak latency, and RMS ratio of each component
            for i=1:N_components
                if i == length(real_idx)
                    peak_amps(C_idx,stim,i:end) = NaN;
                    peak_latens(C_idx,stim,i:end) = NaN;
                    RMSresp(C_idx,stim,i:end) = NaN;
                    RMS(C_idx,stim,i:end) = NaN;
                    break
                end
                component = squeeze(responses(C_idx,stim,real_idx(i):real_idx(i+1)));
                [peak_amp, peak_idx] = max(abs(component));
                peak_loc = peak_idx + real_idx(i) - 1;
                peak_amps(C_idx,stim,i) = component(peak_idx);
                peak_latens(C_idx,stim,i) = peak_loc/fs;
                RMSresp(C_idx,stim,i) = rms(component);
                RMS(C_idx,stim,i) = rms(component)/rms(baseline);
            end           
        end
        
    end
end
save(saveName+suffix,'peak_amps','peak_latens','RMS','RMSoverall','RMSbase','RMSresp')
%%
% overall distribution of peak amplitudes
if fileName == "EP_Measure-3.mat"
    nonzero_peak_amps = peak_amps([3:6,9:end],:,:);
elseif fileName == "EP_Measure-4.mat"
    nonzero_peak_amps = peak_amps([3:10,13:end],:,:);
end
figure()
width = 5e-5;
h1 = histogram(nonzero_peak_amps(:,:,1),'Normalization','probability','BinWidth',width);
hold on
h2 = histogram(nonzero_peak_amps(:,:,2),'Normalization','probability','BinWidth',width); 
h3 = histogram(nonzero_peak_amps(:,:,3),'Normalization','probability','BinWidth',width); 
h4 = histogram(nonzero_peak_amps(:,:,4),'Normalization','probability','BinWidth',width,'FaceColor','#77AC30'); 
xlabel('Amplitude')
ylabel('Proportion')
legend('1','2','3','4')
saveas(gcf,'overalldist_peakamps'+suffix,'png')
component1_mean_amp = nanmean(nonzero_peak_amps(:,:,1),'all')
component2_mean_amp = nanmean(nonzero_peak_amps(:,:,2),'all')
component3_mean_amp = nanmean(nonzero_peak_amps(:,:,3),'all')
component4_mean_amp = nanmean(nonzero_peak_amps(:,:,4),'all')

% overall distribution of peak latencies
if fileName == "EP_Measure-3.mat"
    nonzero_peak_latens = peak_latens([3:6,9:end],:,:);
elseif fileName == "EP_Measure-4.mat"
    nonzero_peak_latens = peak_latens([3:10,13:end],:,:);
end
figure()
width = .02; % 20ms bins
h1 = histogram(nonzero_peak_latens(:,:,1),'Normalization','probability','BinWidth',width); 
hold on
h2 = histogram(nonzero_peak_latens(:,:,2),'Normalization','probability','BinWidth',width); 
h3 = histogram(nonzero_peak_latens(:,:,3),'Normalization','probability','BinWidth',width); 
h4 = histogram(nonzero_peak_latens(:,:,4),'Normalization','probability','BinWidth',width,'FaceColor','#77AC30'); 
xlabel('Latency (s)')
ylabel('Proportion')
legend('1','2','3','4')
saveas(gcf,'overalldist_peaklatens'+suffix,'png')
component1_mean_latency = nanmean(nonzero_peak_latens(:,:,1),'all')
component2_mean_latency = nanmean(nonzero_peak_latens(:,:,2),'all')
component3_mean_latency = nanmean(nonzero_peak_latens(:,:,3),'all')
component4_mean_latency = nanmean(nonzero_peak_latens(:,:,4),'all')

% overall distribution of RMS
if fileName == "EP_Measure-3.mat"
    nonzero_RMS = RMS([3:6,9:end],:,:);
elseif fileName == "EP_Measure-4.mat"
    nonzero_RMS = RMS([3:10,13:end],:,:);
end
figure()
width = .5;
h1 = histogram(nonzero_RMS(:,:,1),'Normalization','probability','BinWidth',width); 
hold on
h2 = histogram(nonzero_RMS(:,:,2),'Normalization','probability','BinWidth',width); 
h3 = histogram(nonzero_RMS(:,:,3),'Normalization','probability','BinWidth',width); 
h4 = histogram(nonzero_RMS(:,:,4),'Normalization','probability','BinWidth',width,'FaceColor','#77AC30'); 
xlabel('RMS ratio')
ylabel('Proportion')
legend('1','2','3','4')
saveas(gcf,'overalldist_rms'+suffix,'png')
component1_mean_rms = nanmean(nonzero_RMS(:,:,1),'all')
component2_mean_rms = nanmean(nonzero_RMS(:,:,2),'all')
component3_mean_rms = nanmean(nonzero_RMS(:,:,3),'all')
component4_mean_rms = nanmean(nonzero_RMS(:,:,4),'all')

% overall distribution of RMSoverall
if fileName == "EP_Measure-3.mat"
    nonzero_RMSoverall = RMSoverall([3:6,9:end],:);
elseif fileName == "EP_Measure-4.mat"
    nonzero_RMSoverall = RMSoverall([3:10,13:end],:);
end
figure()
width = 1;
h1 = histogram(nonzero_RMSoverall,'Normalization','probability','BinWidth',width);
xlabel('RMS ratio')
ylabel('Proportion')
saveas(gcf,'overalldist_rmsoverall'+suffix,'png')

% % channel distribution of peak amplitudes
% figure()
% C_idx = 20; C_r = C_idx+17-1;
% width = 1e-5;
% h1 = histogram(peak_amps(C_idx,:,1),'Normalization','probability','BinWidth',width);
% hold on
% h2 = histogram(peak_amps(C_idx,:,2),'Normalization','probability','BinWidth',width);
% h3 = histogram(peak_amps(C_idx,:,3),'Normalization','probability','BinWidth',width);
% h4 = histogram(peak_amps(C_idx,:,4),'Normalization','probability','BinWidth',width,'FaceColor','#77AC30');
% xlabel('Amplitude')
% ylabel('Proportion')
% legend('1','2','3','4')
% saveas(gcf,'peakamps_ch'+num2str(C_r)+suffix,'png')

% % channel distribution of RMS
% figure()
% C_idx = 20; C_r = C_idx+17-1;
% width = 5e-6;
% h1 = histogram(RMS(C_idx,:,1),'Normalization','probability','BinWidth',width);
% hold on
% h2 = histogram(RMS(C_idx,:,2),'Normalization','probability','BinWidth',width);
% h3 = histogram(RMS(C_idx,:,3),'Normalization','probability','BinWidth',width);
% h4 = histogram(RMS(C_idx,:,4),'Normalization','probability','BinWidth',width,'FaceColor','#77AC30');
% xlabel('RMS')
% ylabel('Proportion')
% legend('1','2','3','4')
% saveas(gcf,'rms_ch'+num2str(C_r)+suffix,'png')
