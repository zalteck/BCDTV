function [CT, M, alpha, beta, gamma ] = BCDHETV(I, RM)
%BCDHE Bayesian TV Color Deconvolution of Hematoxilyn, Eosin slides as
%implemented in:
% Fernando Pérez-Bueno, Miguel López-Pérez, Miguel Vega, Javier Mateos, Valery Naranjo, Rafael Molina, Aggelos K. Katsaggelos,
% A TV-based image processing framework for blind color deconvolution and classification of histological images,
% Digital Signal Processing, 2020, 102727, ISSN 1051-2004, https://doi.org/10.1016/j.dsp.2020.102727.
% (http://www.sciencedirect.com/science/article/pii/S1051200420300725)
%
%   

    [m,n,nc] = size(I);
    tamm = m*n;
    ns = size(RM,2); %number of stains

    if (nc ~= 3 )
        error('Input image does not have 3 channels');
    end

    y2d = intensities2OD( I ); % dimension (m,n,nc)
    YT=reshape(y2d,m*n,nc)';

    clear I y2d 
    
    % Initial values
    CT = RM \ YT;
    CT(CT < eps) = eps;
    M = RM;

    for s=1:ns
        SigmaC = zeros(m*n,ns);
        SigmaM(s) = 0 ;
    end

    % Some stuff
    
    term=1.e-05;
    nitermin=3;
    nitermax=50;
    epsW=mean(CT(:))*1.0e-6;
    
    DhtDh = [-1 2 -1];
    DvtDv = DhtDh';
    
    DhtDh = psf2otf(DhtDh, [m, n]) ;
    DvtDv = psf2otf(DvtDv, [m, n]) ;

    

    CT0 = CT;
    iter = 1; convH = term +1.0; convE = term +1.0;
    %Iterations
    while ( (iter <= nitermin) || (((convH > term) || (convE > term)) && (iter <= nitermax)) )

        % Parameters update
        % beta
        beta = beta_update(YT,CT,SigmaC,M,SigmaM);

        % alpha
        [alpha, WT] = alpha_update_TV(CT,DhtDh,DvtDv,SigmaC,epsW); 

        % gamma
        gamma = gamma_update(M,SigmaM,RM);

        fprintf('iter: %3d\t beta: %f\n',iter,beta)
        fprintf('H\t alpha: %f\t gamma: %f\n',alpha(1),gamma(1))
        fprintf('E\t alpha: %f\t gamma: %f\n',alpha(2),gamma(2))

        % Color vector update
        [M, SigmaM] = color_vector_update(YT,CT,SigmaC,M,RM,beta,gamma);

        % Concentration update
        [CT,SigmaC] = conc_update_TV(YT,CT,M,SigmaM,DhtDh,DvtDv,beta,alpha,m,n,WT); 

        convH = sum((CT(1,:)- CT0(1,:)).*(CT(1,:)- CT0(1,:))) / sum(CT0(1,:).*CT0(1,:));
        convE = sum((CT(2,:)- CT0(2,:)).*(CT(2,:)- CT0(2,:))) / sum(CT0(2,:).*CT0(2,:));
        CT0 = CT;

        fprintf('H\t conv: %e\n',convH)
        fprintf('E\t conv: %e\n',convE)
        M


        iter = iter +1;

    end

    CT(CT < eps) = eps;

end


