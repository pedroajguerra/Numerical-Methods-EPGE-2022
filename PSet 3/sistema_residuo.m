function sp = sistema_residuo(gama,d, parametro) % função que aplica os gamas nos residuos e que eu vou querer encontrar esses gamas que zeram os resíduos
    dim = d + 1;
    matriz_residuo = ones(parametro.n, dim);
    [k_node, k_convertido] = chebroot(d,parametro); % interpolation nodes; k_node = raiz e k_convertido é trazido de volta p grid original
    for ii = 1:parametro.n
       for ij = 1:dim
           matriz_residuo(ii,ij) = residuo(gama, parametro, d, ii, k_convertido(ij));
       end
    end
    sp = reshape(matriz_residuo,[],1); % reshape para transformar num vetor (dica do item 6 - Ana Paula) 
end

