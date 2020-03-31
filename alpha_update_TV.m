function [alpha WT] = alpha_update_TV(CT,DhtDh,DvtDv,SigmaC,epsW)

    ns = size(CT,1);
    tamm = size(CT,2);
    
    alpha = zeros(ns,1);
    for s=1:ns
        tmp = SigmaC(s) .* (DhtDh + DvtDv);
        traza = sum(tmp(:));
        
        [Dhy, Dvy] = circ_gradient(reshape(CT(s,:),size(DhtDh)));
        
        v = Dhy.^2 + Dvy.^2 + traza/tamm;
        v(v<0)=0;
        
        tmp =   v.^ (0.5) ;
        W = 1./ ( epsW + tmp); % These are the Weights
        WT(s,:) = W(:)';
        
        normaptraza = 2 * sum(tmp(:));                
        
        alpha(s) = tamm/(normaptraza + eps);
    end

end
