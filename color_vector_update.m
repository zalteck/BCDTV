function [M, SigmaM] = color_vector_update(YT,CT,SigmaC,M,RM,beta,gamma)


nc = size(YT,1);
ns = size(CT,1);
tamm = size(YT,2);

SigmaM = zeros(ns,1);

for s=1:ns
    [~, eminus] = computingEsZs(YT,CT,M);
    expC2 = (CT(s,:) * CT(s,:)' + sum(SigmaC(:,s)));
    SigmaM(s) = 1.0 ./( beta * expC2 + gamma(s));
    M(:,s) = SigmaM(s) .* (beta * (CT(s,:)*eminus(:,:,s))'  + gamma(s)*RM(:,s)); %%% CUIDADO AQUI. ESTO HAY QUE COMPROBARLO CON LA TEORIA
    % M(M<0)=0;
    scalefactor = norm(M(:,s));
    M(:,s) = M(:,s) / scalefactor;
    SigmaM(s) = SigmaM(s) / (scalefactor.^2);  %%%% IMPORTANTE: Si se normaliza la media hay que normalizar tambi�n la varianza
end

end
