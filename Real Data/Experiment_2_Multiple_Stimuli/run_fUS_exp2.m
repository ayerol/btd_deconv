

% Author: Aybuke Erol (a.erol@tudelft.nl)


% References

% [1] https://www.tensorlab.net/demos/sobi.html#convolutive-mixtures


clear; clc; rng default; close all;
addpath(genpath('../../../../../../Toolboxes/tensorlab4.0beta'))
addpath(genpath('../../utils'))


%% Load the data


load('Exp_Data/fUS_time_series.mat'); % load the fUS time-series
                                        % to a matrix x
roi_lgds= {'Left V2','Left V1','Left SC','Right SC','Right V1','Right V2'}; 
                                      % order of the regions as stored in x

load('Exp_Data/ep.mat'); % the experimental paradigm (ep)
Fs = 3.7202; % sampling rate of the experiment
ep_lgds = {'LM','SL','F','SR','RM'}; % order of the stimulus location 
% timings as stored in ep -->
% Leftmost / Slight-left / Front / Slight-right / Rightmost


%% Parameters


% Specify how many task ('t') and artifact ('a') sources -->
source_list = {'t','t','t','t','t','a'}; 

M = size(x,2);                  % number of regions
R = length(source_list);        % number of sources
L_s = 8;                        % HRF filter length in seconds
L = round(L_s*Fs);              % HRF filter length in samples
K = L;                       	% number of time-lags, i.e. frontal slices
Lacc = 150;                     % window size for matricizing convolutive 
                                                               % mixtures 
numBTDs = 20;                   % number of BTD repetitions
N = length(x);                  % number of time points in the experiment
t_axis = 0:1/Fs:(N-1)/Fs;       % time axis of the experiment


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


%% HRFs


H_est = cell(1,R);

for r = 1:R

    H_est{r} = []; ctr = 1;

    for m = 1:M

        curr_hrf = final_estim_filters{r}(:,m);

        H_est{r} = cat(1,H_est{r},struct_toeplitz(curr_hrf,[],...
            [Lacc L+Lacc],zeros(Lacc-1,1),zeros(Lacc-1,1)));

    end

end

H_all = [];

for r = 1:R

    H_all = cat(2,H_all,H_est{r});

end


%% Source signal estimation


num_task_sources = sum(strcmp(source_list,'t'));

ep_rec = medfilt1(estimate_source(H_all,min(size(H_all)),x_n_cut,N,L,...
    Lacc,R),20);
all_corrs = zeros(num_task_sources,num_task_sources);

for i = 1:num_task_sources
    for j = 1:num_task_sources
        all_corrs(i,j) = corr(ep(:,i),ep_rec(:,j));
    end
end

% % These orderings are determined after analyzing the estimated sources
% % to visualize them with the stimulus condition they match the most -->

estim_source_order = fliplr([3,1,5,4,2]); 
true_source_order = fliplr([1,5,2,4,3]); 
ep_lgds_ordered = ep_lgds(true_source_order);

for r = 1:num_task_sources
    if all_corrs(true_source_order(r),estim_source_order(r)) < 0
        all_corrs(:,estim_source_order(r)) = ...
            -all_corrs(:,estim_source_order(r));
    end
end

% Plot settings -->
ax_pos = ([0.105 0.285 0.465 0.64 0.82]);
seg_len = [.16 .16 .155 .16 .16]; 

figure; cmap = lines(num_task_sources);

for r = 1:num_task_sources

    subplot(num_task_sources,r,1);
    
    curr_ep = ep(:,true_source_order(r));

    if corr(ep_rec(:,estim_source_order(r)),curr_ep) < 0
        ep_rec(:,estim_source_order(r)) = -ep_rec(:,estim_source_order(r));
    end

    p2 = plot(t_axis,ep_rec(:,estim_source_order(r)),...
        'LineWidth',2,'Color','k');
    ymin = min(ep_rec(:,estim_source_order(r))); 
    ymin = ymin - abs(ymin)*0.1; 
    ymax = max(ep_rec(:,estim_source_order(r))); 
    ymax = ymax + ymax*0.1;
    ep_start_times = find(diff(curr_ep)==1)+1;
    stim_on = t_axis(ep_start_times);
    stim_off = t_axis(find(diff(curr_ep)==-1)+1);
    stim_on = stim_on(1:length(stim_off));

    for i = 1:length(stim_on)

        hold on;
        p1 = fill([stim_on(i) stim_on(i) stim_off(i) stim_off(i)],...
            [ymin ymax ymax ymin],cmap(num_task_sources-r+1,:));
        set(p1,'facealpha',.3);
        set(p1,'EdgeColor','none');

    end

    xlim([0 t_axis(end)]); ylim([ymin ymax]); xlim([0 t_axis(end)]);
    set(gca,'YTick',[]); set(gca,'FontSize',15);

    legend([p1 p2],{ep_lgds_ordered{r},...
        ['Estimated Source Signal #' num2str(estim_source_order(r))]},...
        'FontSize',14,'Location','southeast','Orientation','horizontal');

    ax = gca;
    ax.Position = [0.01 ax_pos(r) 0.98 seg_len(r)];

    if r == 1
        xlabel('Time (s)');
    else
        set(gca,'XTick',200:200:1200,'XTickLabel',{});
    end


end


%% Correlation matrix


all_corrs = all_corrs(fliplr(true_source_order),...
    fliplr(estim_source_order)); % flip to place the best results at the 
                                   % beginning of the correlation matrix
figure; imagesc(all_corrs); 
cmap = custom_colormap(-.2,max(all_corrs(:)));
colormap(cmap);
set(gca,'YTick',1:5);
ep_lgds_ordered = fliplr(ep_lgds_ordered);
set(gca,'YTickLabels',ep_lgds_ordered);
set(gca,'XTick',1:5);
set(gca,'XTickLabels',{'#3','#1','#5','#4','#2'});
caxis([-.2 .5]); cb = colorbar;
cb.Ticks = -.2:.2:.4;
ylabel('True Stimulus Times'); xlabel('Estimated Source Signals');
set(gca,'FontSize',16   )
cb.Title.String = 'Pearson Correlation Coefficient';
pos = get(cb,'Position');
cb.Title.Position = [70 170];
cb.Title.Rotation = 90; 
cb.Title.FontSize = 17;

% Compute the false correlations -->
temp = all_corrs; temp(temp<0) = 0; false_corrs = zeros(1,length(temp));

for r = 1:length(temp)

    false_corrs(r) = sum(temp(:,r))-temp(r,r);

end

