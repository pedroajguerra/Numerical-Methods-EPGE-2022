function psi = fun_psi(i,k, parametro)
    if i == 1
        if k >= parametro.grid_k_fea(i) && parametro.grid_k_fea(i+1) >= k
            psi = (parametro.grid_k_fea(i+1) - k)/(parametro.grid_k_fea(i+1) - parametro.grid_k_fea(i)); % divisao no slide 61
        else        
            psi = 0;
        end
        
    elseif i == parametro.ptos    
        if k >= parametro.grid_k_fea(i-1) && parametro.grid_k_fea(i) >= k    
            psi = (k - parametro.grid_k_fea(i-1))/(parametro.grid_k_fea(i) - parametro.grid_k_fea(i-1));
        else        
            psi = 0;        
        end
        
    else        
        if k >= parametro.grid_k_fea(i-1) && parametro.grid_k_fea(i) >= k            
            psi = (k - parametro.grid_k_fea(i-1))/(parametro.grid_k_fea(i) - parametro.grid_k_fea(i-1));            
        elseif k > parametro.grid_k_fea(i) && parametro.grid_k_fea(i+1) >= k
           psi = (parametro.grid_k_fea(i+1) - k)/(parametro.grid_k_fea(i+1) - parametro.grid_k_fea(i));           
        else 
           psi = 0;           
        end        
    end
end

