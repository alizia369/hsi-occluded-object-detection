function [ norMap ] = normaliseMap( Map )
%NORMALISEMAP Summary of this function goes here
%   Detailed explanation goes here
    maxval_mp=max(Map(:));
    minval_mp=min(Map(:));
    norMap=(Map-minval_mp)./(maxval_mp - minval_mp);

end

