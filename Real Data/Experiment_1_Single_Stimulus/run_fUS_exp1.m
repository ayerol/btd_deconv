

% Author: Aybuke Erol (a.erol@tudelft.nl)


% References

% [1] https://www.tensorlab.net/demos/sobi.html#convolutive-mixtures


clear; clc; rng default; close all;
addpath(genpath('../../../../../../Toolboxes/tensorlab4.0beta'))
addpath(genpath('../../utils'))


%% Load the data


load('Exp_Data/fUS_time_series.mat'); % load the fUS time-series
                                        % to a matrix x
lgds = {'SC','LGN','V1'}; % order of the regions as stored in x

load('Exp_Data/ep.mat'); % the experimental paradigm (ep)
Fs = 4; % sampling rate of the experiment


%% Parameters


% Specify how many task ('t') and artifact ('a') sources -->
source_list = {'t','a'};

M = size(x,2);                  % number of regions
R = length(source_list);        % number of sources
L_s = 8;                        % HRF filter length in seconds
L = round(L_s*Fs);              % HRF filter length in samples
K = L;                       	% number of time-lags, i.e. frontal slices
Lacc = ceil(R*L/max(1,M-R));    % window size for matricizing convolutive 
                                  % mixtures 
numBTDs = 20;                   % number of BTD repetitions
N = length(x);                  % number of time points in the experiment
t_axis = 0:1/Fs:(N-1)/Fs;       % time axis of the experiment
cmap = lines(M);                % colormap for visualization


%% Initializations


u = 0:1/Fs:L/Fs; % time-axis of HRFs

all_estim_filters = cell(1,R); % estimated convolutive mixing filters
                                          % at each repetition of BTD
for r = 1:R

    all_estim_filters{r} = zeros(L+1,M,numBTDs);

end

costs = zeros(numBTDs,1); % final value of the BTD objective function

x_ext = zeros(N-Lacc,M*Lacc); % shifted output signal

for m = 1:M

    for l = 0:Lacc-1

        x_ext(:,(m-1)*Lacc+l+1) = x(l+1:end-Lacc+l,m);

    end

end

x_n = x_ext'; x_n_cut = circshift(x_n(:,floor(L/2)+1:end),-floor(L/2));


%% Display raw fUS data

    
ep_start_times = find(diff(ep)==1)+1;
stim_on = t_axis(ep_start_times);
stim_off = t_axis(find(diff(ep)==-1)+1);
stim_on = stim_on(1:length(stim_off));

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

    disp(['Computing BTD #' num2str(testno) ' out of ' num2str(numBTDs)]);
    
    [sol,cost] = btd_deconv(T,M,Lacc,u,source_list);

    for m = 1:M

        for r = 1:R

            temp = sol.factors.(['H' num2str(r)]);
            estim_filter = fliplr(temp(Lacc*(m-1)+1,1:L+1));
            all_estim_filters{r}(:,m,testno) = estim_filter;

        end

    end

end

final_estim_filters = cell(1,R);

for r = 1:R

    final_sol = find_solution(costs,all_estim_filters{r},Fs);
    final_estim_filters{r} = final_sol;

end


%% Display HRFs


H_est = cell(1,R);

for r = 1:R

    H_est{r} = []; ctr = 1;

    if source_list{r} == 't', figure; end

    for m = 1:M

        curr_hrf = final_estim_filters{r}(:,m);

        H_est{r} = cat(1,H_est{r},struct_toeplitz(curr_hrf,[],...
            [Lacc L+Lacc],zeros(Lacc-1,1),zeros(Lacc-1,1)));

        if source_list{r} == 't'

            curr_hrf = normalize(curr_hrf);
            plot(u,curr_hrf-curr_hrf(1),'Color',cmap(ctr,:),...
                'LineWidth',2); 
            hold on; ctr = ctr + 1;

        end

    end

end

set(gca,'FontSize',14); xlabel('Time (s)'); grid on;
title('Estimated HRFs');
lgd  = legend(lgds,'FontSize',14); lgd.EdgeColor = [1 1 1];


%% Source signal estimation


svd_thres = 4; % tolerance for pseudo-inverse of H
ep_rec = estimate_source(H_est{1},svd_thres,x_n_cut,N,L,Lacc,1);
ep_rec = medfilt1(ep_rec,25);

figure;
p2 = plot(t_axis,ep_rec,'LineWidth',2,'Color','k'); 
ymin = -2; ymax = 4;

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
    'DisplayName','Local Maxima of Estimation');


%% Binarization of the source signal


start_idx_rec = zeros(length(locs),1);
end_idx_rec = zeros(length(locs),1);
binarized_ep = zeros(1,N);

% Search for when a drastic change happens in the estimated source -->
for i = 1:length(locs)

    mid_idx = find(t_axis == locs(i));

    for j = 1:100
        if ep_rec(mid_idx) - ep_rec(mid_idx - j) > ep_rec(mid_idx)*.15
            break;
        end
    end

    start_idx_rec(i) = mid_idx - j;

    for j = 1:100
        if ep_rec(mid_idx) - ep_rec(mid_idx + j) > ep_rec(mid_idx)*.15
            break;
        end
    end
    
    end_idx_rec(i) = mid_idx + j;

end

for i = 1:length(locs)

    binarized_ep(start_idx_rec(i):end_idx_rec(i)) = 1;

end

one_period = -5:1/Fs:10; % Construct the time axis for a segment covering 
         % 5 seconds before the stimulus was shown until 10 seconds after

sc_seg = []; lgn_seg = []; v1_seg = []; % Segmented ROI responses
est_ep_seg = []; % Segmented version of the estimated source

for i = ep_start_times

    curr_time_seg = i-5*Fs:i+10*Fs;

    % Segment the regional responses and the binarized source estimation 
    % for later averaging -->
    sc_seg = cat(2,sc_seg,x(curr_time_seg,1));
    lgn_seg = cat(2,lgn_seg,x(curr_time_seg,2));
    v1_seg = cat(2,v1_seg,x(curr_time_seg,3));
    est_ep_seg = cat(2,est_ep_seg,binarized_ep(curr_time_seg)');

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
avg_ep = mean(est_ep_seg,2); thres = .4;
avg_ep(avg_ep<=thres) = 0; avg_ep(avg_ep>thres) = 1; 
start = (find(diff(avg_ep)==1)+1)/Fs-5; 
endd = (find(diff(avg_ep)==-1)+1)/Fs-5;
p2 = fill([start start endd endd],[-1 2.5 2.5 -1],'k');
set(p2,'facealpha',.1); set(p2,'EdgeColor','none');
legend([p1,p2,p3,p4,p5],{'True EP','Estimated EP','SC','LGN','V1'});
xlim([-5,10]);
legend boxoff; xlabel('Time (s)'); yticks(0:2); set(gca,'FontSize',15);
box on;

