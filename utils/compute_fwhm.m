function fwhm = compute_fwhm(data)

halfMax = max(data) / 2; 

index1 = find(data >= halfMax, 1, 'first');

index2 = find(data >= halfMax, 1, 'last');

fwhm = index2-index1 + 1; % in indices

end