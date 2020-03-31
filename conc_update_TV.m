function [CT,SigmaC]  = conc_update_TV(YT,CT,M,SigmaM,DhtDh,DvtDv,beta,alpha,m,n,WT)

    ns = size(CT,1);
    tamm = size(YT,2);
    
    eps_cs=1.0e-07;
    itmax_cs=50;

    SigmaC = zeros(tamm,ns);
    
  
    for s=1:ns  
        [zminus, ~] = computingEsZs(YT,CT,M);
        W = reshape(WT(s,:)',m,n);
        expM2 = M(:,s)' * M(:,s) + 3.0 * SigmaM(s);
        auxSigmaC = 1.0 ./( beta * expM2 + alpha(s) * mean(W(:))*(DhtDh + DvtDv));
        SigmaC(:,s) = auxSigmaC(:);
        
        inv_cov_fix = beta * expM2;
        indep_term = beta * zminus(:,s);
        
        
        fprintf('\tGradiente conjugado s = %d\n', s);
        
        inv_cov = get_mic_handle(inv_cov_fix,alpha(s),W,m,n);
        
        [cs, flag_cs,relres_cs,iter_cs ]=pcg(inv_cov, indep_term, eps_cs, itmax_cs,[], [], CT(s,:)' );
        fprintf('\tFlag: %d, RelRes: %d, Iters: %d\t\n',flag_cs,relres_cs,iter_cs);
        
        CT(s,:) = cs(:)';
        %CT(CT < eps) = eps; %%%%% Yo no forzaba no negatividad en cada vuelta. La forzaba al final, en todo caso
    end
end

function h = get_mic_handle(inv_cov_fix,alpha,W,nr,nc)
            h = @mic;
            function SigmaCscs = mic(cs)
                SigmaCscs = multiply_by_invcov(cs, inv_cov_fix,alpha,W,nr,nc);
            end
end

function SigmaCscs = multiply_by_invcov(cs,inv_cov_fix,alpha,W,nr,nc)

        F1 = inv_cov_fix* cs;
        
        %% Prior
        
        cs2d=reshape(cs,nr,nc);       
        
        [Dhcs, Dvcs] = circ_gradient(cs2d);

        [F2, ~] = Tcirc_gradient( W.* Dhcs );

        [temp,F3] = Tcirc_gradient( W.* Dvcs );
        
        priorTerm = reshape(alpha * (F2 + F3) , nr*nc, 1);

        SigmaCscs = F1 + priorTerm;
    
    end

