function [raiz, raiz_gridk] = chebroot(d, parametro) 
    dim = d + 1;
    upper_bound_k = max(parametro.grid_k);
    lower_bound_k = min(parametro.grid_k);
    raiz = zeros(1,dim);
    raiz_gridk = zeros(1,dim); % vou usar isso para retomar para o grid de capital original 
    for i = 1:dim
       raiz(i) = - cos(((2*i - 1)*pi)/(2*dim));
       raiz_gridk(i) = ((upper_bound_k - lower_bound_k)/2)*(raiz(i) + 1) + lower_bound_k;
    end
end

