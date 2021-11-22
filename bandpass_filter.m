function [xo] = bandpass_filter(x, f1, f2, filter_order, Fs)
% *WAVE* 
%
% BANDPASS FILTER    filter a datacube between frequencies f1 and f2
%
% INPUT: 
% x - datacube (r,c,t)
% f1 - low-frequency cutoff
% f2 - high-frequency cutoff
% filter_order - filter order (N.B. doubled for forward-reverse filter)
% Fs - sampling frequency
%
% OUTPUT
% xo - output datacube
%

assert( ndims(x) == 3, 'datacube input required' );
xo = zeros( size(x) );

% construct filter
ct = [f1 f2];
ct = ct / (Fs/2);
[b,a] = butter( filter_order, ct ) ;

% proceed with filtering
for rr = 1:size(x,1)
    for cc = 1:size(x,2)
%         if(x(rr,cc,:)==zeros(1,1,size(x,3)))
%             disp('a')
%             xo(rr,cc,:) = zeros(size(x,3),1);
%             continue
%         end
        xo(rr,cc,:) = filtfilt( b, a, squeeze(x(rr,cc,:)) );
    end
end
