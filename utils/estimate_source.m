function ep_rec = estimate_source(H,svd_thres,x_n,N,L,Lacc)

[~,S,~] = svd(H); 
tol_pinv = S(svd_thres,svd_thres);
ep_n_estim = pinv(H,tol_pinv)*x_n;

ep_rec = zeros(N,1);

% Ideally, the estimated source matrix should be block-Hankel, meaning the 

% off-diagonal terms should be equal. As such, we take the mean of these 

% terms to reconstruct back a vector that represents the underlying source.

for l = 1:(L+Lacc)

    ep_rec(l) = mean(diag(ep_n_estim(min(l,size(ep_n_estim,1)):-1:1,...
        1:min(l,size(ep_n_estim,2)))));

end

start = L+Lacc;

for l=1:N-(L+Lacc)-1

    ep_rec(start+l) = mean(diag(ep_n_estim(start:-1:start-(L+Lacc)+1,...
        l+1:min(size(ep_n_estim,2),l+L+Lacc))));

end

ep_rec(end) = ep_n_estim(start,end);

end