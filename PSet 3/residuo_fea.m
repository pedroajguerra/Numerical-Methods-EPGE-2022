function res_fea = residuo_fea(k, a, j, parametro)
    c_old_fea = consumo_colfea(k, a(j,:), parametro);
    g_fea = parametro.grid_z(j)*(k^parametro.alfa) + (1 - parametro.delta)*k - c_old_fea;
    parcela_1_fea = zeros(parametro.n,1);
    parcela_2_fea = zeros(parametro.n,1);
    for i = 1:parametro.n
        parcela_1_fea(i) = 1 - parametro.delta + parametro.alfa*parametro.grid_z(i)*g_fea.^(parametro.alfa - 1);
        pol_c_fea = consumo_colfea(g_fea, a(i,:), parametro);
        parcela_2_fea(i) = (pol_c_fea./c_old_fea).^(-parametro.mu);
    end   
    res_fea = parametro.beta.*parametro.P(j,:)*(parcela_1_fea.*parcela_2_fea) - 1;    
end

