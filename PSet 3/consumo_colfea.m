function cfea = consumo_colfea(k, a, parametro)
 
soma = 0;
    for ss = 1:parametro.ptos
        soma = soma + a(ss)*fun_psi(ss, k, parametro);
    end 

cfea = soma;
end

