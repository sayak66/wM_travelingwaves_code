function [D,cc,pv] = phase_distance_map(source, np1,np2 )
% *WAVE*
%
% PHASE CORRELATION DISTANCE    correlation of phase with distance
%                                   (circular-linear), given an input 
%                                   phase map
%
% INPUT
% pl - phase map (r,c)
% source - source point (sc)
% spacing - pixel spacing (sc)
%
% OUTPUT
% cc - circular-linear correlation coefficient, phase correlation w/ distance
% pv - p-value of the correlation (H0: rho == 0, H1: rho != 0)
%

pl = ones(np1,np2);
spacing = 1;

% make matrix of distances from wave center
[r,c] = size(pl); [X,Y] = meshgrid( (1:c)-source(1), (1:r)-source(2) );
D = sqrt( X.^2 + Y.^2 ); D = D .* spacing;

% flatten, remove NaNs
D = D(:); pl = pl(:); D( isnan(pl) ) = []; pl( isnan(pl) ) = [];

D = reshape(D,np1,np2);
