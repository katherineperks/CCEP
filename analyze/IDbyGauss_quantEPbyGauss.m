% this method performs baseline subtraction,
% identifies EP components using the zero-crossing indices of the sum of Gaussians fit on each trial, 
% and quantifies each component using the parameters from a Gaussian fit
clear
subject = "0d5e8e";
fileName = "EP_Measure-4.mat";
cd("C:\Users\Katherine\Documents\AAA UW\Herron Lab\Data\" + subject)
load(fileName)
saveName = "IDbyGauss_quantEPbyGauss";
if ~exist(saveName, 'dir')
   mkdir(saveName)
end
cd(saveName)

if fileName == "EP_Measure-3.mat"
    suffix = "_stim23_24";
elseif fileName == "EP_Measure-4.mat"
    suffix = "_stim27_28";
end
load("..\IDbyGauss_quantEPbyRMS\IDbyGauss"+suffix+".mat")

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
% flags
showplot = false;
saveplot = true; % plots should be identical to IDbyGauss_quantEPbyRMS plots since using same Gaussian fits
if fileName == "EP_Measure-3.mat"
    suffix = "_stim23_24";
elseif fileName == "EP_Measure-4.mat"
    suffix = "_stim27_28";
end
% function to find zero-crossing indices
zci = @(v) find(v.*circshift(v,1) <= 0);
N_components = 4;
% quantification metrics
Gauss_parameters = zeros(N_chansub,N_stims,N_components,4); % channel x trial x component x parameter
Gauss_residuals = zeros(N_chansub,N_stims,N_components); % channel x trial x component's residual
Gauss_corrcoefs = zeros(N_chansub,N_stims,N_components); % channel x trial x component's goodness of fit
component_idx = zeros(N_chansub,N_stims,N_components+1); % channel x trial x component's start/end idx
opts1 = optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-10);

d = 1:(window_sz-floor(0.01*fs)+1); % exclude the first 10ms of response from quantification
Gausfitfun = @(x) ( x(3)*exp(-1/2*((d-x(2))/x(1)).^2) + x(6)*exp(-1/2*((d-x(5))/x(4)).^2) + x(9)*exp(-1/2*((d-x(8))/x(7)).^2) + x(12)*exp(-1/2*((d-x(11))/x(10)).^2) + x(13) );
for C_idx=1:N_chansub % response channel index (out of 20)
    C_r = C_idx+17-1; % response channel number (out of 126)
    if sum(avgresponses(C_idx,:) == zeros(1,window_sz)) == window_sz
        continue
    end
    for stim=1:N_stims
       % identify EP components from zci of sum of Gaussians fit
       y(:) = responses(C_idx,stim,floor(0.01*fs):end); % exclude the first 10ms of response from quantification
       Gausfit = Gausfitfun(all_parameters(C_idx,stim,:));
       idx = zci(Gausfit);
       if isempty(idx)
            Gauss_parameters(C_idx,stim,:,:) = NaN;
            Gauss_residuals(C_idx,stim,:) = NaN;
            Gauss_corrcoefs(C_idx,stim,:) = NaN;
            component_idx(C_idx,stim,:) = NaN;
            continue
       end
       if idx(1) == 1
           idx = idx(2:end);
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
        plot(floor(0.01*fs):window_sz,Gausfit) % sum of Gaussians fit (starts at 10ms)
        plot(real_idx,zeros(length(real_idx)),'ro')
        xline(0.01*fs,'r') % 10ms
        xline(0.5*fs,'r') % 500ms
        yline(0,'r')
        if saveplot && stim==N_stims/2
            figName = "components_ch"+num2str(C_r)+suffix+"_trial"+num2str(stim);
            saveas(gcf,figName,'png')
        end
                      
        % fit Gaussian to each component and save best-fit parameters
        for i=1:N_components
            if i == length(real_idx)
                Gauss_parameters(C_idx,stim,i:end,:) = NaN;
                Gauss_residuals(C_idx,stim,i:end) = NaN;
                Gauss_corrcoefs(C_idx,stim,i:end) = NaN;
                component_idx(C_idx,stim,i+1:end) = NaN;
                break
            end
            temp_residuals = zeros(1,3);
            temp_parameters = zeros(3,4);
            component = squeeze(responses(C_idx,stim,real_idx(i):real_idx(i+1)))';
            component_idx(C_idx,stim,i+1) = real_idx(i+1);
            t = 1:length(component);
            fun = @(x) ( x(1)*exp(-0.5*(t-x(2)).^2./(x(3)^2)) + x(4) ) - (component);
            for g=1:3
                % random guesses for each parameter: x1=amplitude, x2=mu, x3=sigma, x4=baseline
                [peak,x2] = max(abs(component)); x1 = component(x2); x3 = round(length(component)/4); x4 = 0;
                x0 = [x1, x2, x3, x4];
                if x1 >= 0
                    lb = [0, 1, 1, 0];
                    ub = [peak, length(component), length(component), 0];
                else
                    lb = [-peak, 1, 1, 0];
                    ub = [0, length(component), length(component), 0];
                end
                [x,temp_residuals(g)] = lsqnonlin(fun,x0,lb,ub,opts1);
                temp_parameters(g,:) = x;
            end
            [Gauss_residuals(C_idx,stim,i), bestfit_index] = min(temp_residuals);
            Gauss_parameters(C_idx,stim,i,:) = temp_parameters(bestfit_index,:);
            fit = fun(Gauss_parameters(C_idx,stim,i,:))+component;
            R = corrcoef(fit,component);
            Gauss_corrcoefs(C_idx,stim,i) = R(1,2);
            % convert mu, sigma from units of samples to seconds
            Gauss_parameters(C_idx,stim,i,2) = Gauss_parameters(C_idx,stim,i,2)/fs;
            Gauss_parameters(C_idx,stim,i,3) = Gauss_parameters(C_idx,stim,i,3)/fs;
        end           
    end
    save(saveName+suffix,'Gauss_parameters','Gauss_residuals','Gauss_corrcoefs','component_idx')
end
save(saveName+suffix,'Gauss_parameters','Gauss_residuals','Gauss_corrcoefs','component_idx')
%%
% Gauss_parameters(:,:,:,2) = Gauss_parameters(:,:,:,2)/fs;
% Gauss_parameters(:,:,:,3) = Gauss_parameters(:,:,:,3)/fs;
% overall distribution of amplitude
if fileName == "EP_Measure-3.mat"
    nonzero_Gauss_parameters = Gauss_parameters([3:6,9:end],:,:,:);
elseif fileName == "EP_Measure-4.mat"
    nonzero_Gauss_parameters = Gauss_parameters([3:10,13:end],:,:,:);
end
figure()
width = 5e-5;
h1 = histogram(nonzero_Gauss_parameters(:,:,1,1),'Normalization','probability','BinWidth',width);
hold on
h2 = histogram(nonzero_Gauss_parameters(:,:,2,1),'Normalization','probability','BinWidth',width); 
h3 = histogram(nonzero_Gauss_parameters(:,:,3,1),'Normalization','probability','BinWidth',width); 
h4 = histogram(nonzero_Gauss_parameters(:,:,4,1),'Normalization','probability','BinWidth',width,'FaceColor','#77AC30'); 
xlabel('Amplitude')
ylabel('Proportion')
legend('1','2','3','4')
saveas(gcf,'overalldist_amp'+suffix,'png')
component1_mean_amp = nanmean(nonzero_Gauss_parameters(:,:,1,1),'all')
component2_mean_amp = nanmean(nonzero_Gauss_parameters(:,:,2,1),'all')
component3_mean_amp = nanmean(nonzero_Gauss_parameters(:,:,3,1),'all')
component4_mean_amp = nanmean(nonzero_Gauss_parameters(:,:,4,1),'all')

% overall distribution of mu
figure()
width = .02; % 20ms bins
h1 = histogram(nonzero_Gauss_parameters(:,:,1,2),'Normalization','probability','BinWidth',width); 
hold on
h2 = histogram(nonzero_Gauss_parameters(:,:,2,2),'Normalization','probability','BinWidth',width); 
h3 = histogram(nonzero_Gauss_parameters(:,:,3,2),'Normalization','probability','BinWidth',width); 
h4 = histogram(nonzero_Gauss_parameters(:,:,4,2),'Normalization','probability','BinWidth',width,'FaceColor','#77AC30'); 
xlabel('Mu (s)')
ylabel('Proportion')
legend('1','2','3','4')
saveas(gcf,'overalldist_mu'+suffix,'png')
component1_mean_mu = nanmean(nonzero_Gauss_parameters(:,:,1,2),'all')
component2_mean_mu = nanmean(nonzero_Gauss_parameters(:,:,2,2),'all')
component3_mean_mu = nanmean(nonzero_Gauss_parameters(:,:,3,2),'all')
component4_mean_mu = nanmean(nonzero_Gauss_parameters(:,:,4,2),'all')

% overall distribution of sigma
figure()
width = .01; % 10ms bins
h1 = histogram(nonzero_Gauss_parameters(:,:,1,3),'Normalization','probability','BinWidth',width); 
hold on
h2 = histogram(nonzero_Gauss_parameters(:,:,2,3),'Normalization','probability','BinWidth',width); 
h3 = histogram(nonzero_Gauss_parameters(:,:,3,3),'Normalization','probability','BinWidth',width); 
h4 = histogram(nonzero_Gauss_parameters(:,:,4,3),'Normalization','probability','BinWidth',width,'FaceColor','#77AC30'); 
xlabel('Sigma (s)')
ylabel('Proportion')
legend('1','2','3','4')
saveas(gcf,'overalldist_sigma'+suffix,'png')
component1_mean_sigma = nanmean(nonzero_Gauss_parameters(:,:,1,3),'all')
component2_mean_sigma = nanmean(nonzero_Gauss_parameters(:,:,2,3),'all')
component3_mean_sigma = nanmean(nonzero_Gauss_parameters(:,:,3,3),'all')
component4_mean_sigma = nanmean(nonzero_Gauss_parameters(:,:,4,3),'all')
