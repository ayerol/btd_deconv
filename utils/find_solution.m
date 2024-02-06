
% find_solution.m returns the final solution of the block-term
% decomposition after clustering of individual runs
%
% INPUT:
%   costs            : objective value of each run
%   all_estim_fitlers: estimated hemodynamic response function of each run
%   Fs               : sampling rate
% 
% OUTPUT:
%   final_sol        : final estimation for the hemodynamic response
                       % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function final_sol = find_solution(costs,all_estim_filters,Fs)


M = size(all_estim_filters,2);
numBTDs = size(all_estim_filters,3);
estimated_hrf_peak_lats = zeros(numBTDs,M);

for testno = 1:numBTDs

    for m = 1:M

        curr_hrf = all_estim_filters(:,m,testno);
        [~,mx_idx] = max(curr_hrf);
        estimated_hrf_peak_lats(testno,m) = mx_idx;

    end

end

test_idx = 1:numBTDs;

% Eliminate the high-cost solutions

imb = imbinarize(costs);

if sum(imb) ~= length(costs)

    estimated_hrf_peak_lats = estimated_hrf_peak_lats(~imb,:);
    test_idx = test_idx(~imb);

end

del = [];

% Eliminate the solutions having a peak latency higher than 4.5 seconds 
% (which is outside the expected range)

for i = 1:size(estimated_hrf_peak_lats,1)

    if sum(estimated_hrf_peak_lats(i,:) > 4.5*Fs) > 0 

        del = cat(2,del,i);

    end

end

estimated_hrf_peak_lats(del,:) = [];
test_idx(del) = [];

% Cluster the remaining solutions

c = cluster(linkage(estimated_hrf_peak_lats),'maxClust',...
    floor(size(estimated_hrf_peak_lats,1)/3));

classes = unique(c);
dists = zeros(1,length(classes));
mean_clusters = zeros(M,length(classes));

for cidx = 1:length(classes)

    curr_cluster = (estimated_hrf_peak_lats(c == classes(cidx),:))';

    if size(curr_cluster,2) == 1 % ignore if there is only one sample in a 
                                                                 % cluster

        dists(cidx) = 1000;
        continue;

    end

    mean_cluster = mean(curr_cluster,2);
    mean_clusters(:,cidx) = mean_cluster;

    % Calculate the intra-cluster distance metric

    dists(cidx) =  1/(sum(c == classes(cidx))) * ...
        max(pdist2(curr_cluster',mean_cluster'));

end

dists(dists == 0) = 1000;
[~,minc] = min(dists);

% Final HRF is equal to the mean of the HRFs belonging to the cluster with
% minimum intra-cluster distance

final_sol = mean(all_estim_filters(:,:,test_idx(c == minc)),3);


end