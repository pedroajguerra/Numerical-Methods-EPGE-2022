function res = residuo(gama, parametro, d, j, k_aux)
    c_old = consumo(gama(j,:), d, parametro, k_aux);
    g = parametro.grid_z(j)*(k_aux^parametro.alfa) + (1 - parametro.delta)*k_aux - c_old;
    parcela_1 = zeros(parametro.n,1);
    parcela_2 = zeros(parametro.n,1);
    for i = 1:parametro.n
        parcela_1(i) = 1 - parametro.delta + parametro.alfa*parametro.grid_z(i)*g.^(parametro.alfa - 1);
        pol_c = consumo(gama(i,:), d, parametro, g);
        parcela_2(i) = (pol_c./c_old).^(-parametro.mu);
    end   
    res = parametro.beta.*parametro.P(j,:)*(parcela_1.*parcela_2) - 1;
end

