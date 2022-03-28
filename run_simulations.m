

% Author: Aybuke Erol (a.erol@tudelft.nl)


% References

% [1] https://www.tensorlab.net/demos/sobi.html#convolutive-mixtures

% [2] Correa, N., Adali, T., Yi-Ou Li, Calhoun, V.D., 2005. Comparison of 

% blind source separation algorithms for fMRI using a new matlab toolbox: 

% Gift, in: Proceedings. (ICASSP  05). IEEE International Conference on 

% Acoustics, Speech, and Signal Processing, 2005, pp. v/401–v/404


clear; clc; close all; rng default;
addpath(genpath('../../../../../Toolboxes/tensorlab4.0beta'))
addpath(genpath('utils'))


%% Parameters


Fs = 2;                 % sampling rate
M = 3;                  % number of regions (i.e., number of measurements)
R = 2;                  % number of sources
min_rest = 10;          % minimum rest duration in seconds
max_rest = 15;          % maximum rest duration in seconds
stim_dur = 4;           % stimulus duration in seconds
L_s = 8;                % HRF filter length in seconds
L = L_s*Fs;             % HRF filter length in samples
K = L;                  % number of time-lags, i.e., frontal slices
Lacc = ceil(R*L/(M-R)); % shifting for convolutive mixtures
numReps = 20;           % number of stimulus repetitions
snr_val = 0;            % signal-to-noise ratio of measurements (dB)
numRuns = 100;          % number of monte-carlo iterations
numBTDs = 20;           % number of BTD repetitions
svd_thres = 5;          % tolerance for pseudo-inverse of H


%% Initializations


u = 0:1/Fs:L/Fs; % time-axis of HRFs
H = zeros(L+1,M,R); % stores the simulated HRF forms

costs = zeros(numBTDs,1); % final value of the BTD objective function

all_estim_hrfs = zeros(L+1,M,numBTDs); % all of the estimated HRFs within a 
                                         % Monte-Carlo run

true_hrf_peak_lats = zeros(M,1); % true HRF peak latencies    

final_hrf_error = zeros(numRuns,1); % average (for m = 1,...M) peak latency
                                      % error
                                                                    
ious = zeros(1,numRuns); % intersection-over-union (IoU) values between the 
                           % true and estimated source signals
                                                                    
ious_fixedHrf = zeros(1,numRuns); % for comparison of the IoU values when
                                     % a fixed HRF was assumed

cmap = lines(M); % colormap for visualization

construct_hrf = @(z) (z(2)*(((z(3)^z(1))*(u.^(z(1)-1)).*exp(-z(3)*u))/...
    gamma(z(1)))); % gamma-model for generating the HRFs

% Construct an HRF with fixed parameters to compare EP estimation results
                            % (use the default canonical model parameters)
fixedHrf = 1/gamma(6)*u.^(6-1).*exp(-u)-1/6*1/gamma(16)*u.^(16-1).*exp(-u);


%% Monte-carlo runs


for norun = 1:numRuns


    %% Generate sources


    % Generate the EP (surce signal of interest): a binary signal which is
    % equal to 1 when the stimulus is on, and 0 when off
    
    ep = zeros(round(10*Fs),1); % start with a rest period of 10 seconds

    for i = 1:numReps

        ep = cat(1,ep,ones(stim_dur*Fs,1));
        rest = round(((max_rest-min_rest) * rand + min_rest)*Fs);
        ep = cat(1,ep,zeros(rest,1));

    end

    N = length(ep);

    % For generation of the noise source, refer to [2]
    noise_source = 0.1 + 0.2* randn([N 1]);
    noise_source(N-floor(N/2):N) = noise_source(N-floor(N/2):N) + 0.02;
    
    s = [ep noise_source];


    %% Generate convolutive mixing filters


    % Noise_source will be additive, thus will be convolved with an impulse
    impulse = zeros(L+1,1);
    impulse(floor((L+1)/2)+1) = 1;

    zs = [2+4.75*rand(M,1) rand(M,1) 1.25+3.25*rand(M,1)]; % randomly
    % generate HRF parameters in a predefined range

    for m = 1:M

        for r = 1:R

            if r == 1 % Source of interest

                z = zs(m,:);

                h = construct_hrf(z);
                H(:,m,r) = h;

                [mx,max_idx] = max(h);
                true_hrf_peak_lats(m) = max_idx; % in samples

            end

            if r == 2 % Noise/Artifact source

                H(:,m,r) = m * impulse;

            end

        end

    end


    %% Generate measurements


    x = zeros(N,M);
    sig_pwr = zeros(1,M); noise_pwr = zeros(1,M);

    for r = 1:R

        for m = 1:M

            if r == 1

                temp = conv(s(:,r),H(:,m,r));
                x(:,m) = temp(1:N);
                sig_pwr(m) = var(x(:,m));

            else

                % The added noise power should be according to the 
                % specified SNR value

                noise_pwr(m) = sig_pwr(:,m)/var(s(:,r))/db2pow(snr_val);
                x(:,m) = x(:,m) + sqrt(noise_pwr(m))*s(:,r);

            end

        end

    end

    N = N-10; t_axis = 0:1/Fs:(N-1)/Fs;
    ep = ep(1:N); x = x(1:N,:); s = s(1:N,:);

    stim_on = t_axis(find(diff(ep)==1)+1);
    stim_off = t_axis(find(diff(ep)==-1)+1);

    x_ext = zeros(N-Lacc,M*Lacc); % shifted output signal

    for m = 1:M

        for l = 0:Lacc-1

            x_ext(:,(m-1)*Lacc+l+1) = x(l+1:end-Lacc+l,m);

        end

    end

    x_n = x_ext'; x_n_cut = circshift(x_n(:,floor(L/2)+1:end),-floor(L/2));
    T = scov(x_ext,0:K-1); % tensor of lagged output autocorrelations


    %% Perform BTD


    estimated_hrf_peak_lats = zeros(M,numBTDs); % estimated HRF peak 
    % latencies at each BTD repetition

    for testno = 1:numBTDs

        [sol,cost] = btd_deconv(T,M,R,Lacc,u);

        for m = 1:M

            estim_hrf = fliplr(sol.factors.H1(Lacc*(m-1)+1,1:L+1));

            all_estim_hrfs(:,m,testno) = estim_hrf; 

            [~,max_idx] = max(estim_hrf);

            estimated_hrf_peak_lats(m,testno) = max_idx; % in samples

        end

        costs(testno) = cost;

    end

    final_hrfs = find_solution(costs,all_estim_hrfs,Fs);

    for m = 1:M

        curr_hrf = final_hrfs(:,m);
        [~,max_idx] = max(curr_hrf);
        final_hrf_error(norun) = final_hrf_error(norun) + ...
            abs(true_hrf_peak_lats(m) - max_idx); % Add up the error of all 
                                                                  % regions

    end

    final_hrf_error(norun) = final_hrf_error(norun)/M/Fs; % average over 
    % regions and convert unit to seconds


    % Visualization of the estimated HRFs (gives same plot as in the paper)
    if norun == 2 && snr_val == 0; plot_hrfs(final_hrfs,H,u,cmap); end


    %% Source Signal Estimation
    

    H_T_est = []; % block column of H that is of interest

    for m = 1:M

        H_T_est = cat(1,H_T_est,struct_toeplitz(final_hrfs(:,m),[],...
            [Lacc L+Lacc],zeros(Lacc-1,1),zeros(Lacc-1,1)));

    end

    
    ep_rec = medfilt1(estimate_source...
        (H_T_est,svd_thres,x_n_cut,N,L,Lacc),12); % estimated ep
    
    iou = compute_iou(ep_rec,t_axis,stim_on,stim_off); % ratio (over 1)

    % Normalize the iou ratio with respect to the stimulus duration
    iou = iou * stim_dur; % over 4 (in seconds)

    ious(norun) = mean(iou); % average over stimulus repetitions

    
    % Do the same for the fixed HRF case

    H_T_fixed_est = [];

    for m = 1:M

        H_T_fixed_est = cat(1,H_T_fixed_est,struct_toeplitz(fixedHrf,...
            [],[Lacc L+Lacc],zeros(Lacc-1,1),zeros(Lacc-1,1)));

    end

    ep_rec_fixedHrf = medfilt1(estimate_source...
        (H_T_fixed_est,svd_thres,x_n_cut,N,L,Lacc),12);

    iou_fixedHrf = compute_iou(ep_rec_fixedHrf,t_axis,stim_on,stim_off);

    iou_fixedHrf = iou_fixedHrf * stim_dur; 

    ious_fixedHrf(norun) = mean(iou_fixedHrf); 

    
    % Visualization of the estimated ep (gives same plot as in the paper)
    if norun == 1 && snr_val == 0
        plot_source(ep_rec,stim_on,stim_off,t_axis); 
    end


end


if snr_val == 0 % Display the HRF peak latency score at 0 dB SNR

    fprintf('\nMedian HRF Peak Latency Error: %.1f \n',...
        median(final_hrf_error)); % median across Monte-Carlo iterations
    
    fprintf('Standard Deviation of the HRF Peak Latency Error: %.1f \n',...
        std(final_hrf_error)) % std across Monte-Carlo iterations

end


if ~exist('Results', 'dir'); mkdir('Results'); end

save(['Results/ious_' num2str(snr_val) 'dB.mat'],'ious'); 
save(['Results/ious_fixedHrf' num2str(snr_val) 'dB.mat'],'ious_fixedHrf'); 

