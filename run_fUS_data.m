

% Author: Aybuke Erol (a.erol@tudelft.nl)


% References

% [1] https://www.tensorlab.net/demos/sobi.html#convolutive-mixtures


clear; clc; rng default; close all;
addpath(genpath('../../../../../Toolboxes/tensorlab4.0beta'))
addpath(genpath('utils'))


%% Load the data


load('fUS_Data/fUS_time_series.mat'); % Load the fUS time-series into a 
                                                   % matrix called as x
lgds = {'SC','LGN','V1'}; % Order of the regions as stored in x


%% Parameters


Fs = 4;                         % sampling rate of the experiment (fixed, 
                                  % should not be changed)
M = size(x,2);                  % number of regions
R = 2;                          % nb sources
L_s = 8;                        % HRF filter length in seconds
L = L_s*Fs;                     % HRF filter length in samples
K = L;                       	% number of time-lags, i.e., frontal slices
Lacc = ceil(R*L/(M-R));         % shifting for convolutive mixtures
numBTDs = 20;                   % number of BTD repetitions
N = length(x);                  % number of time points in the experiment


%% Load the experimental paradigm (ep)


% The visual stimulus was displayed with a frame rate of 25 fps. Thus, the
% experimental paradigm should be interpolated to the time axis of the
% experiment (at Fs = 4)

movie_fps = 25;
load('fUS_Data/ep_original.mat'); % the ep at its original time scale 

screen_frame_times = 0:1/movie_fps:(length(ep_original)-1)/movie_fps;
t_axis = 0:1/Fs:(N-1)/Fs; % time-axis of the experiment

ep = interp1(screen_frame_times,ep_original,t_axis,'nearest');


% Find the true start and end times of the ep, will be used later

true_start_times = find(diff(ep)==1)+1;
stim_on = t_axis(true_start_times);
stim_off = t_axis(find(diff(ep)==-1)+1);
stim_on = stim_on(1:length(stim_off));


%% Initializations


u = 0:1/Fs:L/Fs; % time-axis of HRFs
cmap = lines(M); % colormap for visualization

all_estim_hrfs = zeros(L+1,M,numBTDs); % estimated HRFs at each repetition 
                                                                  % of BTD
costs = zeros(numBTDs,1); % final value of the BTD objective function

x_ext = zeros(N-Lacc,M*Lacc); % shifted output signal

for m = 1:M

    for l = 0:Lacc-1

        x_ext(:,(m-1)*Lacc+l+1) = x(l+1:end-Lacc+l,m);

    end

end

x_n = x_ext'; x_n_cut = circshift(x_n(:,floor(L/2)+1:end),-floor(L/2));


%% Display raw fUS data


fig = figure; ymin = -3; ymax = 10; yaxis = 0:4:8;

for m = 1:M

    subplot(M,1,m);
    p1 = plot(t_axis,x(:,m),'LineWidth',2,'Color',cmap(m,:));
    for i = 1:length(stim_on)
        hold on;
        p4 = fill([stim_on(i) stim_on(i) stim_off(i) stim_off(i)],...
            [ymin ymax ymax ymin],'r');
        set(p4,'facealpha',.1);
        set(p4,'EdgeColor','none');
    end
    xlim([0 t_axis(end)]); ylim([ymin ymax]); xlim([0 t_axis(end)]);
    set(gca,'YTick',yaxis); set(gca,'FontSize',15); 
    title(lgds{m}); ax = gca;

    switch m

        case 1
            ax.Position = [0.1 0.71 0.88 0.245];
            set(gca,'XTick',[]); 
        case 2
            ax.Position = [0.1 0.41 0.88 0.245];
            set(gca,'XTick',[]); 
        case 3
            ax.Position = [0.1 0.11 0.88 0.245];
            xlabel('Time (s)');

    end

end

han=axes(fig,'visible','off'); han.YLabel.Visible='on';
ylbl = ylabel(han,'Power Doppler Amplitude [a.u.]','FontSize',18);
ylbl.Position(2) = 0.5;
ylbl.Position(1) = -0.1;


%% Perform BTD


T = scov(x_ext,0:K-1); % tensor of lagged output autocorrelations

for testno = 1:numBTDs

    [sol,cost] = btd_deconv(T,M,R,Lacc,u);

    for m = 1:M

        estim_hrf = fliplr(sol.factors.H1(Lacc*(m-1)+1,1:L+1));

        all_estim_hrfs(:,m,testno) = estim_hrf;

    end

end

final_hrfs = find_solution(costs,all_estim_hrfs,Fs);


%% Display HRFs


H_T_est = [];
widths = zeros(M,1); % full-width-at-half-maximum for each estimated HRF
peak_lats = zeros(M,1); % peak latency for each estimated HRF

for m = 1:M
    
    curr_hrf = normalize(final_hrfs(:,m));

    H_T_est = cat(1,H_T_est,struct_toeplitz(final_hrfs(:,m),[],...
        [Lacc L+Lacc],zeros(Lacc-1,1),zeros(Lacc-1,1)));
        
    widths(m) = compute_fwhm(curr_hrf)/Fs; % in seconds
    [~,max_idx] = max(curr_hrf);
    peak_lats(m) = max_idx/Fs; % in seconds

    plot(u,curr_hrf - curr_hrf(1),'LineWidth',2); hold on;

end

set(gca,'FontSize',14); xlabel('Time (s)'); legend(lgds,'FontSize',15);
legend boxoff; yticks(0:3); grid on;


%% Source Signal Estimation


svd_thres = 4; % tolerance for pseudo-inverse of H
ep_rec = medfilt1(estimate_source(H_T_est,svd_thres,x_n_cut,N,L,Lacc),25);


p2 = plot(t_axis,ep_rec,'LineWidth',2,'Color','k'); ymin = -2; ymax = 4;

for i = 1:length(stim_on)

    hold on;
    p1 = fill([stim_on(i) stim_on(i) stim_off(i) stim_off(i)],...
        [ymin ymax ymax ymin],'r');
    set(p1,'facealpha',.1);
    set(p1,'EdgeColor','none');

end

xlim([0 t_axis(end)]); ylim([ymin ymax]); xlim([0 t_axis(end)]);
set(gca,'YTick',[]); set(gca,'FontSize',15); xlabel('Time (s)'); hold on;
legend([p1 p2],{'Experimental Paradigm','Estimated Source Signal'},...
    'FontSize',15);
thres = 0.35; hold on; 
plot(t_axis,thres*ones(1,N),'--k','DisplayName','Global Threshold');
hold on; [pks,locs] = findpeaks(ep_rec,t_axis);
plot(locs+1,pks+.2,'rv','MarkerFaceColor','red',...
    'DisplayName','Local Maxima of Estimation')


start_idx_rec = zeros(length(locs),1);
end_idx_rec = zeros(length(locs),1);
binarized_ep = zeros(1,N);

for i = 1:length(locs)

    mid_idx = find(t_axis == locs(i));

    for j = 1:100
        if ep_rec(mid_idx) - ep_rec(mid_idx - j) > ep_rec*.15
            break;
        end
    end

    start_idx_rec(i) = mid_idx - j;

    for j = 1:100
        if ep_rec(mid_idx) - ep_rec(mid_idx + j) > ep_rec*.15
            break;
        end
    end
    
    end_idx_rec(i) = mid_idx + j;

end

for i = 1:length(locs)

    binarized_ep(start_idx_rec(i):end_idx_rec(i)) = 1;

end


one_period = -5:1/Fs:10; % construct the time axis for 5 seconds before the 
                                    % stimulus shown until 10 seconds after

sc_seg = []; lgn_seg = []; v1_seg = [];

for i = true_start_times

    curr_time_seg = i-5*Fs:i+10*Fs;

    % Segment the regional responses for later averaging
    sc_seg = cat(2,sc_seg,x(curr_time_seg,1));
    lgn_seg = cat(2,lgn_seg,x(curr_time_seg,2));
    v1_seg = cat(2,v1_seg,x(curr_time_seg,3));

end

figure; 
hold on; sc = mean(sc_seg,2);
p3 = plot(one_period,sc-sc(1),'LineWidth',2,'Color',cmap(1,:));
hold on; lgn = mean(lgn_seg,2);
p4 = plot(one_period,lgn-lgn(1),'LineWidth',2,'Color',cmap(2,:));
hold on; v1 = mean(v1_seg,2);
p5 = plot(one_period,v1-v1(1),'LineWidth',2,'Color',cmap(3,:));
p1 = fill([0 0 4 4],[-1 2.5 2.5 -1],'r'); hold on;
set(p1,'facealpha',.1); set(p1,'EdgeColor','none');
avg_ep = mean(ep_seg,2); 
avg_ep(avg_ep<=.65) = 0; avg_ep(avg_ep>.65) = 1; 
start = find(diff(avg_ep)==1)/Fs-5; endd = find(diff(avg_ep)==-1)/Fs-5;
p2 = fill([start start endd endd],[-1 2.5 2.5 -1],'k');
set(p2,'facealpha',.1); set(p2,'EdgeColor','none');
legend([p1,p2,p3,p4,p5],{'True EP','Estimated EP','SC','LGN','V1'});
xlim([-5,10]);
legend boxoff; xlabel('Time (s)'); yticks(0:2); set(gca,'FontSize',15);
box on;

