function c = consumo(gama, d, parametro, k_aux)
    % Normalizando os valores do capital    
    %grid_k_aux = parametro.grid_k;
    k_max = max(parametro.grid_k);
    k_min = min(parametro.grid_k);
    teste = 2*((k_aux-k_min)/(k_max-k_min)) - 1;
 
    % Criando a função consumo
    soma = 0;
    dim = d + 1;
    for j = 1:dim
        soma = soma + gama(j)*cheb(j-1,teste);
    end
    
    c = soma;
end

