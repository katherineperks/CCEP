% this method performs baseline subtraction,
% identifies EP components using the zero-crossing indices of the sum of Gaussians fit on each trial, 
% and quantifies each component using peak amplitude, peak latency and RMS
clear
subject = "0d5e8e";
fileName = "EP_Measure-4.mat";
cd("C:\Users\Katherine\Documents\AAA UW\Herron Lab\Data\" + subject)
load(fileName)
saveName = "IDbyGauss_quantEPbyRMS";
if ~exist(saveName, 'dir')
   mkdir(saveName)
end
cd(saveName)

if fileName == "EP_Measure-3.mat"
    suffix = "_stim23_24";
elseif fileName == "EP_Measure-4.mat"
    suffix = "_stim27_28";
end
load("IDbyGauss"+suffix+".mat")

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
        % subtract mean of 200ms period before stim onset (within channel)
        responses(C_idx,stim,:) = squeeze(responses(C_idx,stim,:) - mean(usedata(start-floor(0.2*fs):start-1,C_idx)));
    end
    avgresponses(C_idx,:) = mean(responses(C_idx,:,:),2);
end
%%
% if fileName == "EP_Measure-3.mat"
%     suffix = "_stim23_24";
% elseif fileName == "EP_Measure-4.mat"
%     suffix = "_stim27_28";
% end
% % store best fit solution, RSS, and R for each channel and each trial
% all_parameters = zeros(N_chansub,N_stims,13); % channel x trial x parameter
% all_residuals = zeros(N_chansub,N_stims); % channel x trial's residual
% all_corrcoefs = zeros(N_chansub,N_stims); % channel x trial's goodness of fit
% opts1 = optimoptions('lsqnonlin','Display','off');
% 
% d = 1:(window_sz-floor(0.01*fs)+1); % exclude the first 10ms of response from quantification
% for C_idx=1:N_chansub % response channel index (out of 20)
%     C_r = C_idx+17-1; % response channel number (out of 126)
%     % analyze channels with non-zero responses
%     if sum(avgresponses(C_idx,:) == zeros(1,window_sz)) == window_sz
%         continue
%     end
%     % fit sum of Gaussians on a single trial basis
%     for stim=1:N_stims
%        temp_residuals = zeros(1,5);
%        temp_parameters = zeros(5,13);
%        y(:) = responses(C_idx,stim,floor(0.01*fs):end); % exclude the first 10ms of response from quantification
%        fun = @(x) ( x(3)*exp(-1/2*((d-x(2))/x(1)).^2) + x(6)*exp(-1/2*((d-x(5))/x(4)).^2) + x(9)*exp(-1/2*((d-x(8))/x(7)).^2) + x(12)*exp(-1/2*((d-x(11))/x(10)).^2) + x(13) ) - (y);
%        for g=1:5
%            % random guesses for each parameter, parameter bounds consistent across Gaussian components:
%            % sigma between 10ms to .25s
%            % mu between 10ms to 1.5s (all mu values below account for exclusion of first 10ms)
%            % amplitude between +/-1e-2 = +/-0.01
%            % baseline parameter included but not really interpretable since already baseline-subtracted
%            x1 = randi([floor(0.01*fs) floor(0.25*fs)]);   x2 = 0.03;  x3 =  1e-3*rand(); % x1=sigma,   x2=mu,  x3=amplitude for component 1
%            x4 = randi([floor(0.01*fs) floor(0.25*fs)]);   x5 = 0.08;  x6 = -1e-4*rand(); % x4=sigma,   x5=mu,  x6=amplitude for component 2
%            x7 = randi([floor(0.01*fs) floor(0.25*fs)]);   x8 = 0.2;   x9 =  5e-4*rand(); % x7=sigma,   x8=mu,  x9=amplitude for component 3
%            x10 = randi([floor(0.01*fs) floor(0.25*fs)]); x11 = 0.4;  x12 = -5e-4*rand(); % x10=sigma, x11=mu, x12=amplitude for component 4
%            x13 = -2e-4+(2e-4+2e-4).*rand(); % x13=baseline shared across components
%            x0 = [x1, round((x2-0.01)*fs), x3, x4, round((x5-0.01)*fs), x6, x7, round((x8-0.01)*fs), x9, x10, round((x11-0.01)*fs), x12, x13];
%            lb = [floor(0.01*fs), 1, -1e-2, ...
%                  floor(0.01*fs), 1, -1e-2, ...
%                  floor(0.01*fs), 1, -1e-2, ...
%                  floor(0.01*fs), 1, -1e-2, ...
%                  -1e-4];
%            ub = [floor(0.25*fs), length(y), 1e-2, ...
%                  floor(0.25*fs), length(y), 1e-2, ...
%                  floor(0.25*fs), length(y), 1e-2, ...
%                  floor(0.25*fs), length(y), 1e-2, ...
%                  1e-4];
%            [x,temp_residuals(g)] = lsqnonlin(fun,x0,lb,ub,opts1);
%            temp_parameters(g,:) = x;
%        end
%        [all_residuals(C_idx,stim),bestfit_index] = min(temp_residuals);
%        all_parameters(C_idx,stim,:) = temp_parameters(bestfit_index,:);
%        fit = fun(all_parameters(C_idx,stim,:))+y;
%        actual = y;
%        R = corrcoef(fit,actual);
%        all_corrcoefs(C_idx,stim) = R(1,2);
%     end
%     save("IDbyGauss"+suffix,'all_parameters','all_residuals','all_corrcoefs')
% end
% save("IDbyGauss"+suffix,'all_parameters','all_residuals','all_corrcoefs')
%%
% flags
showplot = false;
saveplot = true;
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
RMS = zeros(N_chansub,N_stims,N_components);
component_idx = zeros(N_chansub,N_stims,N_components+1); % channel x trial x component's start/end idx

d = 1:(window_sz-floor(0.01*fs)+1); % exclude the first 10ms of response from quantification
fitfun = @(x) ( x(3)*exp(-1/2*((d-x(2))/x(1)).^2) + x(6)*exp(-1/2*((d-x(5))/x(4)).^2) + x(9)*exp(-1/2*((d-x(8))/x(7)).^2) + x(12)*exp(-1/2*((d-x(11))/x(10)).^2) + x(13) );
for C_idx=1:N_chansub % response channel index (out of 20)
    C_r = C_idx+17-1; % response channel number (out of 126)
    if sum(avgresponses(C_idx,:) == zeros(1,window_sz)) == window_sz
        continue
    end
    for stim=1:N_stims
       % identify EP components from zci of sum of Gaussians fit
       y(:) = responses(C_idx,stim,floor(0.01*fs):end); % exclude the first 10ms of response from quantification
       fit = fitfun(all_parameters(C_idx,stim,:));
       idx = zci(fit);
       if isempty(idx)
            peak_amps(C_idx,stim,:) = NaN;
            peak_latens(C_idx,stim,:) = NaN;
            RMS(C_idx,stim,:) = NaN;
            component_idx(C_idx,stim,:) = NaN;
            continue
       end
       real_idx = [0, idx] + floor(0.01*fs);
       real_idx = [real_idx(diff(real_idx)>=5), real_idx(end)]; % component windows must contain at least 6 datapoints
       component_idx(C_idx,stim,1) = floor(0.01*fs); % first component always starts at 10ms after stim
       
        if showplot && rem(stim,50)==0
            f = figure();
        else
            f = figure('Visible','off');
        end
        plot(1:window_sz,squeeze(responses(C_idx,stim,:)));
        hold on
        plot(floor(0.01*fs):window_sz,fit) % sum of Gaussians fit (starts at 10ms)
        plot(real_idx,zeros(length(real_idx)),'ro')
        xline(0.01*fs,'r') % 10ms
        xline(0.5*fs,'r') % 500ms
        yline(0,'r')
        yline(all_parameters(C_idx,stim,13),'k--')
        if saveplot && stim==N_stims/2
            figName = "components_ch"+num2str(C_r)+suffix+"_trial"+num2str(stim);
            saveas(gcf,figName,'png')
        end
                 
        % extract single-trial peak amplitude, peak latency, and RMS of each component
        for i=1:N_components
            if i == length(real_idx)
                peak_amps(C_idx,stim,i:end) = NaN;
                peak_latens(C_idx,stim,i:end) = NaN;
                RMS(C_idx,stim,i:end) = NaN;
                component_idx(C_idx,stim,i+1:end) = NaN;
                break
            end
            component = squeeze(responses(C_idx,stim,real_idx(i):real_idx(i+1)));
            component_idx(C_idx,stim,i+1) = real_idx(i+1);
            [peak_amp, peak_idx] = max(abs(component));
            peak_loc = peak_idx + real_idx(i) - 1;
            peak_amps(C_idx,stim,i) = component(peak_idx);
            peak_latens(C_idx,stim,i) = peak_loc/fs; % convert peak latency from units of samples to seconds
            RMS(C_idx,stim,i) = rms(component);
        end           
    end
    save(saveName+suffix,'peak_amps','peak_latens','RMS','component_idx')
end
save(saveName+suffix,'peak_amps','peak_latens','RMS','component_idx')
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
width = 2e-5;
h1 = histogram(nonzero_RMS(:,:,1),'Normalization','probability','BinWidth',width); 
hold on
h2 = histogram(nonzero_RMS(:,:,2),'Normalization','probability','BinWidth',width); 
h3 = histogram(nonzero_RMS(:,:,3),'Normalization','probability','BinWidth',width); 
h4 = histogram(nonzero_RMS(:,:,4),'Normalization','probability','BinWidth',width,'FaceColor','#77AC30'); 
xlabel('RMS')
ylabel('Proportion')
legend('1','2','3','4')
saveas(gcf,'overalldist_rms'+suffix,'png')
component1_mean_rms = nanmean(nonzero_RMS(:,:,1),'all')
component2_mean_rms = nanmean(nonzero_RMS(:,:,2),'all')
component3_mean_rms = nanmean(nonzero_RMS(:,:,3),'all')
component4_mean_rms = nanmean(nonzero_RMS(:,:,4),'all')

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
