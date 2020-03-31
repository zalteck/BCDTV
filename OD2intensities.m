function [ I ] = OD2intensities( y )
%OD2intensities Convert an optical density image y to a RGB slice image I
%   Detailed explanation goes here

    I = 10.^(-y);
%    I = exp(-y);
end

