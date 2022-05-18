function ep_rec = estimate_source(H,svd_thres,x_n,N,L,Lacc,R)

[~,S,~] = svd(H); 
tol_pinv = S(svd_thres,svd_thres);
ep_n_estim = pinv(H,tol_pinv)*x_n;

ep_rec = zeros(N,R);

% Ideally, the estimated source matrix should be block-Hankel, meaning the 

% off-diagonal terms should be equal. As such, we take the mean of these 

% terms to reconstruct back a vector that represents the underlying source.

for r = 1:R

    start = (r-1)*(Lacc+L);
    for l = 1:(L+Lacc)

        ep_rec(l,r) = mean(diag(ep_n_estim(min(start+l,...
            size(ep_n_estim,1)):-1:start+1,1:min(l,size(ep_n_estim,2)))));

    end

    start = start+L+Lacc;
    for l=1:N-(L+Lacc)-1

        ep_rec(start-(r-1)*(L+Lacc)+l,r) = mean(diag(ep_n_estim(start:...
            -1:start-(L+Lacc)+1,l+1:min(size(ep_n_estim,2),l+L+Lacc))));

    end

    ep_rec(end,r) = ep_n_estim(start,end);

end

end