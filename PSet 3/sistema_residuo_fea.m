function sp_fea = sistema_residuo_fea(a, parametro)
    matriz_residuo_fea = ones(parametro.n, parametro.ptos);
    for ii = 1:parametro.n
       for ij = 1:parametro.ptos
           matriz_residuo_fea(ii,ij) = residuo_fea(parametro.grid_k_fea(ij), a, ii, parametro);
       end
    end
    sp_fea = reshape(matriz_residuo_fea,[],1); % reshape para transformar num vetor (dica do item 6 - Ana Paula)
end

