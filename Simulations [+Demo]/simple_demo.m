

% Author: Aybuke Erol (a.erol@tudelft.nl)


% References

% [1] https://www.tensorlab.net/demos/sobi.html#convolutive-mixtures


clear; clc; close all; rng default;
addpath(genpath('../../../../../Toolboxes/tensorlab4.0beta'))
addpath(genpath('../utils'))


%% Load example simulated data


load('Example_Data/hrfs.mat'); % load the true HRFs (H is an L x M matrix
% where L is the HRF length and M is the number of regions)

Fs = 2; % sampling rate of the simulated data

L = size(H,1)-1; % corresponds to L/Fs = 8 seconds
M = size(H,2); % there are 3 regions

u = 0:1/Fs:L/Fs; % time-axis of HRFs


cmap = lines(M); % colormap for visualization
figure(1);

for m = 1:M

    subplot(M,1,m);
    plot(u,H(:,m)/max(H(:,m)),'--','LineWidth',2,...
        'Color',cmap(m,:));
    set(gca,'YTick',[]);
    legend(['True HRF, m = ', num2str(m)],'FontSize',15);
    legend boxoff;
    ax = gca; set(gca,'FontSize',15);
    if m == 1
        set(gca,'XTick',[]);
        ax.Position = [0.01 0.705 0.98 0.285];
    elseif m == 2
        set(gca,'XTick',[]);
        ax.Position = [0.01 0.405 0.98 0.285];
    else
        set(gca,'XTick',0:8);
        ax.Position = [0.01 0.105 0.98 0.285];
    end

end
xlabel('Time (s)');


load('Example_Data/source_matrix.mat'); % there are two generated sources:

% the first column represents the event timings (i.e. experimental paradigm

% - ep) that we are interested in (and we will try to estimate it assuming

% it is unknown), whereas the second column stands for the artifact source)


N = length(s); % length of the experiment
t_axis = 0:1/Fs:(N-1)/Fs; % time-axis of the experiment
ep = s(:,1); 
figure(2); stem(t_axis,ep); xlabel('Time (s)'); title('True EP');


%% Parameters


% Specify how many task ('t') and artifact ('a') sources -->
source_list = {'t','a'};

R = size(s,2);          % number of sources
K = L;                  % number of time-lags, i.e., frontal slices
Lacc = ceil(R*L/(M-R)); % shifting for convolutive mixtures
numBTDs = 20;           % number of BTD repetitions
svd_thres = 5;          % tolerance for pseudo-inverse of H
snr_val = 0;            % desired value for the signal-to-noise ratio of 
                          % measurements (in dB)


%% Initializations


costs = zeros(numBTDs,1); % final value of the BTD objective function

all_estim_hrfs = zeros(L+1,M,numBTDs); % all of the estimated HRFs


%% Generate measurement signals


x = zeros(N,M);
sig_pwr = zeros(1,M); noise_pwr = zeros(1,M);

for r = 1:R

    for m = 1:M

        if source_list{r} == 't'

            temp = conv(s(:,r),H(:,m));
            x(:,m) = temp(1:N);
            sig_pwr(m) = var(x(:,m));

        else

            % The added noise power should be according to the specified
            % SNR value

            noise_pwr(m) = sig_pwr(:,m)/var(s(:,r))/db2pow(snr_val);
            x(:,m) = x(:,m) + sqrt(noise_pwr(m))*s(:,r);

        end

    end

end


% Visualization of the measured responses 

% (Note that, the proposed solution assumes only these responses are known,

% and aims to recover back both the sources of interest (Figure 2) and the

% HRFs (Figure 1)) -->

figure;

for m = 1:M

    subplot(M,1,m);
    plot(t_axis,x(:,m),'LineWidth',2,'Color',cmap(m,:));
    set(gca,'YTick',[]);
    legend(['Measured Response, m = ', num2str(m)],'FontSize',15);
    legend boxoff;
    ax = gca; set(gca,'FontSize',15);
    if m == 1
        set(gca,'XTick',[]);
        ax.Position = [0.01 0.705 0.98 0.285];
    elseif m == 2
        set(gca,'XTick',[]);
        ax.Position = [0.01 0.405 0.98 0.285];
    else
        set(gca,'XTick',0:50:t_axis(end));
        ax.Position = [0.01 0.105 0.98 0.285];
    end

end
xlabel('Time (s)');


%% Compute the autocorrelation tensor


x_ext = zeros(N-Lacc,M*Lacc); % shifted output signal

for m = 1:M

    for l = 0:Lacc-1

        x_ext(:,(m-1)*Lacc+l+1) = x(l+1:end-Lacc+l,m);

    end

end

x_n = x_ext'; x_n_cut = circshift(x_n(:,floor(L/2)+1:end),-floor(L/2));
T = scov(x_ext,0:K-1); % tensor of lagged output autocorrelations


%% Perform BTD


for testno = 1:numBTDs

    [sol,cost] = btd_deconv(T,M,Lacc,u,source_list);

    for m = 1:M

        estim_hrf = fliplr(sol.factors.H1(Lacc*(m-1)+1,1:L+1));

        all_estim_hrfs(:,m,testno) = estim_hrf;

    end

    costs(testno) = cost;

end

final_hrfs = find_solution(costs,all_estim_hrfs,Fs);


% Visualization of the estimated HRFs (gives same plot as in the paper) -->

plot_hrfs(final_hrfs,H,u,cmap);


%% Source signal estimation


H_T_est = []; % block column of H that is of interest

for m = 1:M

    H_T_est = cat(1,H_T_est,struct_toeplitz(final_hrfs(:,m),[],...
        [Lacc L+Lacc],zeros(Lacc-1,1),zeros(Lacc-1,1)));

end

ep_rec = medfilt1(estimate_source...
    (H_T_est,svd_thres,x_n_cut,N,L,Lacc,1),12); % estimated ep


% Visualization of the estimated ep (gives same plot as in the paper) -->

stim_on = t_axis(find(diff(ep)==1)+1);
stim_off = t_axis(find(diff(ep)==-1)+1);

plot_source(ep_rec,stim_on,stim_off,t_axis);

