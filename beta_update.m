function beta = beta_update(YT,CT,SigmaC,M,SigmaM)

    tamm = size(YT,2);


    norma = trace((YT-M*CT)*(YT-M*CT)');
    trazadiagonal = diag(sum(SigmaC));
    traza1 = trace(M' * M * trazadiagonal);
    traza2 = trace((CT * CT' + trazadiagonal) * diag(3 * SigmaM));
    auxbeta = 3*tamm / ( norma + traza1 + traza2);

    beta = auxbeta;

end

