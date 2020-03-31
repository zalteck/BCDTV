function gamma = gamma_update(M,SigmaM,RM)


    aux = M-RM + eps;
    aux2 = aux' * aux;
    auxgamma = 3.0 ./ diag(aux2 + diag(3.0 * SigmaM));

    gamma = auxgamma;

end