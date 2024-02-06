
% compute_fwhm.m returns the full-width-at-half-maximum of a given function
%
% INPUT:
%   data: hemodynamic response function 
%
% OUTPUT:
%   fwhm: full-width-at-half-maximum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fwhm = compute_fwhm(data)

halfMax = max(data) / 2; 

index1 = find(data >= halfMax, 1, 'first');

index2 = find(data >= halfMax, 1, 'last');

fwhm = index2-index1 + 1; % in indices

end