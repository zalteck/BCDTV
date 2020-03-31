function [ y ] = intensities2OD( I )
%intensities2OD Convert a RGB slice image I to its corresponding optical
%density (OD) representation y
%   Detailed explanation goes here
    
    I(I<eps(0)) = min(I(I>0));

    y = - log10(I);
%    y = - log(double(I));
end

