% This example shows how to use the method proposed in:
% Fernando Pérez-Bueno, Miguel López-Pérez, Miguel Vega, Javier Mateos, Valery Naranjo, Rafael Molina, Aggelos K. Katsaggelos,
% A TV-based image processing framework for blind color deconvolution and classification of histological images,
% Digital Signal Processing, 2020, 102727, ISSN 1051-2004, https://doi.org/10.1016/j.dsp.2020.102727.
% (http://www.sciencedirect.com/science/article/pii/S1051200420300725)
%
% 
%% Load image and reference vectors
clc,clear all
I = imread('histWB.jpg');
load 'RMImageSet' RM;
[m,n,nc] = size(I);
subplot(231),imshow(I)
title('Original H&E Image')
%% Deconvolution

[CT, M, alpha, beta, gamma] = BCDHETV(im2double(I), RM);
disp('completed')

%% Band visualization (OD space)

ns = size(M,2)
concentrations = reshape(CT',m,n,ns);

%figure()
subplot(232),imshow(concentrations(:,:,1))
title('OD H Band')
subplot(235),imshow(concentrations(:,:,2))
title('OD E Band')


%% Band reconstruction (RGB space)
Hrec_OD = reshape((M(:,1)*CT(1,:))',m,n,nc);
Hrec_RGB = OD2intensities(Hrec_OD);

Erec_OD = reshape((M(:,2)*CT(2,:))',m,n,nc);
Erec_RGB = OD2intensities(Erec_OD);

%figure()
subplot(233),imshow(Hrec_RGB)
title('RGB H Band')
subplot(236),imshow(Erec_RGB)
title('RGB E Band')


