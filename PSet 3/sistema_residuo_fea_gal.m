function sp_fea_gal = sistema_residuo_fea_gal(a, parametro)

matriz_residuo_fea_gal = zeros(parametro.int, parametro.n);
[raiz, ~]=chebroot(parametro.int-1,parametro);
e1=zeros(parametro.int,parametro.n);
e2=zeros(parametro.int,parametro.n);

    for ii=1:parametro.n
        for ij=1:parametro.int-1
            k_aux=(raiz(ij)+1).*(parametro.grid_k_fea(ij+1)-parametro.grid_k_fea(ij))./2 + parametro.grid_k_fea(ij); 
            parcela1=residuo_fea(k_aux,a,ii,parametro);
            parcela2=sqrt(1-raiz(ij).^2).*(parametro.grid_k_fea(ij+1)-k_aux)./(parametro.grid_k_fea(ij+1)-parametro.grid_k_fea(ij));
            e1(ij,ii)=(pi*(parametro.grid_k_fea(ij+1)-parametro.grid_k_fea(ij))/(2*(parametro.int-1))).*(parcela1')*parcela2;
        end
    end
    
    for ii=1:parametro.n
        for ij=2:parametro.int
            k_aux=(raiz(ij)+1).*(parametro.grid_k_fea(ij)-parametro.grid_k_fea(ij-1))./2 + parametro.grid_k_fea(ij-1); 
            parcela1=residuo_fea(k_aux,a,ii,parametro);
            parcela2=sqrt(1-raiz(ij).^2).*(k_aux-parametro.grid_k_fea(ij-1))./(parametro.grid_k_fea(ij)-parametro.grid_k_fea(ij-1));
            e2(ij,ii)=(pi*(parametro.grid_k_fea(ij)-parametro.grid_k_fea(ij-1))/(2*(parametro.int-1))).*(parcela1')*parcela2; 
        end
    end

matriz_residuo_fea_gal = e1+e2;
sp_fea_gal = reshape(matriz_residuo_fea_gal, [], 1); % reshape para transformar num vetor (dica do item 6 - Ana Paula)
    
end

