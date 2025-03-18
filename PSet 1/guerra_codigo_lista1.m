clear
%%% Numerical Methods (prof. Cézar Santos)
%%% Problem Set 1 - Pedro Augusto Januzzi Guerra

% Parameters:
rho = 0.95;     
sigma = 0.007;
mu = 0; %mean of error % note this is also the mean of the AR(1) process
m = 3;

%%% Item (1)

% Create grid of z
n = 9;
upper_bound = m*sigma*(sqrt(1-rho^2))^-1;
lower_bound = -upper_bound;

grid_z = linspace(lower_bound, upper_bound,n);

% Calculate transition matrix
delta_z = (upper_bound - lower_bound)*((n-1)^-1);
P = zeros(n,n);

% Calculate the first and last columns of matrix
P(:,1) = normcdf((grid_z(1)-rho*grid_z+delta_z/2)/sigma);
P(:,n) = 1 - normcdf((grid_z(n)-rho*grid_z-delta_z/2)/sigma);

% Calculate other entries of the matrix
for j=2:(n-1)
    P(:,j)= normcdf((grid_z(j)+(delta_z/2)-rho*grid_z)/sigma) - normcdf((grid_z(j)-(delta_z/2)-rho*grid_z)/sigma);
end

% % Alternative method using normcdf but without standardizing the distribution
% P = zeros(n,n);
% P(:,1) = normcdf((grid_z(1)-rho*grid_z+delta_z/2),mu,sigma);
% P(:,n) = 1 - normcdf((grid_z(n)-rho*grid_z-delta_z/2),mu,sigma);
% 
% %  Calculate other entries of the matrix
% for j=2:(n-1)
%     P(:,j)= normcdf((grid_z(j)+(delta_z/2)-rho*grid_z),mu,sigma) - normcdf((grid_z(j)-(delta_z/2)-rho*grid_z),mu,sigma);
% end

%%% Item (2)

% Create grid of z
sigma_z = sqrt(sigma^2*((1-rho^2)^-1));
upper_bound_rou = sigma_z*sqrt(n-1);
lower_bound_rou = -upper_bound_rou;
grid_z_rou = linspace(lower_bound_rou, upper_bound_rou,n); % rou refers to Rouwenhorst

p = (1+rho)/2; % variable used to compute P recursively

P_rou = [p 1-p; 1-p p]; %matrix P2
%P_rou3 = p*[P_rou zeros(2,1); zeros(1,2) 0] + (1-p)*[zeros(2,1) P_rou; zeros(1,2) 0] + (1-p)*[zeros(1,2) 0; P_rou zeros(2,1)] + p*[zeros(1,2) 0; zeros(2,1) P_rou]

% Computando a matriz P_rou recursivamente
for i=3:n
    P_rou = p*[P_rou zeros(i-1,1); zeros(1,i-1) 0] + (1-p)*[zeros(i-1,1) P_rou; zeros(1,i-1) 0] + ...
        (1-p)*[zeros(1,i-1) 0; P_rou zeros(i-1,1)] + p*[zeros(1,i-1) 0; zeros(i-1,1) P_rou];  
end

% Para normalizar as linhas:
for x=1:n
    P_rou(x,:) = P_rou(x,:)/sum(P_rou(x,:));
end

% Matriz de transição
P_rou;
    
%%% Item (3)

rng(1); % botando um seed para padronizar o processo gerador de dados aleatórios
t = 10000; % número de períodos
choque = normrnd(mu,sigma,1,t); %gerando choques aleatórios normalmente distribuídos (média 0, desvpad 0.007)

% Simulando o processo AR(1)
z_t = zeros(1,t);

for i=2:t
    z_t(i) = rho*z_t(i-1)+choque(i);
end

% Processo AR(1)
z_t;

% Gráfico AR(1)
plot(z_t, 'b')
title ('Processo AR(1) simulado para 10000 períodos', 'FontSize', 34)
xlabel('Períodos', 'FontSize', 28)
ylabel('z_t', 'FontSize', 28)

% Simulando o processo discretizado via método de Tauchen

% Pegando o ponto no grid_z em que z = 0 (i.e., a mediana)

mediana = median(grid_z, 'all');

% Encontrando o estado inicial

estado = find(grid_z==mediana);

% Criando um vetor para armazenar os estados em que estou (começando da mediana = 0)
z_discretizado = [mediana zeros(1,t-1)];

for i=2:t
    cdf_estado = P(estado,1); %considerando que estou no estado inicial igual a 5, qual a probabilidade de eu ir para cada estado
    for j = 1:n
        if normcdf(choque(i)/sigma) < cdf_estado
           estado = j;
           z_discretizado(i) = grid_z(estado); % vetor que acumula para qual estado eu fui dado o choque recebido (z_t discretizado)
           break
        else 
            cdf_estado = P(estado,j+1) + cdf_estado;
        end  
    end    
end    

% Processo discretizado (Tauchen)
z_discretizado;

% Criando o gráfico do proceso discretizado (Tauchen)
plot(z_discretizado, 'r')
title ('Processo discretizado via Tauchen', 'FontSize', 34)
xlabel('Períodos', 'FontSize', 28)
ylabel('z discretizado', 'FontSize', 28)

% Plotando a simulação do AR(1) e do processo discretizado via Tauchen juntos
plot(z_t,'b')
hold on
title ('AR(1) com processo discretizado via Tauchen', 'FontSize', 34)
plot(z_discretizado,'r')
xlabel('Períodos', 'FontSize', 28)
ylabel('z', 'FontSize', 28)
hold off
legend('AR(1)', 'Discretizado (Tauchen)')

% Simulando o processo discretizado via método de Rouwenhorst

mediana_rou = median(grid_z_rou, 'all');
estado_rou = find(grid_z_rou == mediana_rou);

% Criando um vetor para armazenar os estados em que estou(começando da mediana = 0)
z_discretizado_rou = [mediana_rou zeros(1,t-1)];

for i=2:t
    cdf_estado_rou = P_rou(estado_rou,1); %considerando que estou no estado inicial igual a 5, qual probabilidade de eu ir para outro estado
    for j = 1:n
        if normcdf(choque(i)/sigma) < cdf_estado_rou
           estado_rou = j;
           z_discretizado_rou(i) = grid_z_rou(estado_rou); % vetor que acumula para qual estado eu fui dado o choque recebido (z_t discretizado)
           break
        else 
            cdf_estado_rou = P_rou(estado_rou,j+1) + cdf_estado_rou;
        end  
    end    
end    

% Processo discretizado(Rouwenhorst)
z_discretizado_rou;

% Criando o gráfico do proceso discretizado (Tauchen)
plot(z_discretizado_rou,'r')
title ('Processo discretizado via Rouwenhorst', 'FontSize', 34)
xlabel('Períodos', 'FontSize', 28)
ylabel('z discretizado', 'FontSize', 28)

% Plotando a simulação do AR(1) e do processo discretizado via Rouwenhorst juntos
plot(z_t,'b')
hold on
title ('AR(1) com processo discretizado via Rouwenhorst', 'FontSize', 34)
plot(z_discretizado_rou,'r')
xlabel('Períodos', 'FontSize', 28)
ylabel('z', 'FontSize', 28)
hold off
legend('AR(1)', 'Discretizado (Rouwenhorst)')

%%% Item (4)

% Rodando a regressão do z_discretizado no seu próprio lag (método de Tauchen)

reg_tauchen = fitlm(lagmatrix(z_discretizado,1),z_discretizado,'Intercept',false);

% Rodando a regressão do z_discretizado no seu próprio lag (método de Rouwenhorst)

reg_rou = fitlm(lagmatrix(z_discretizado_rou,1),z_discretizado_rou,'Intercept',false);

% O coeficiente da regressão é muito próximo de 0.95 (rho), então estão bem próximos do processo gerador de dados real.
% Repare que usando o método de Rouwenhorst, a aproximação do processo real é melhor.
% Se testarmos rho = 0.7, o método de Tauchen fica mais próximo.

% Calculando o intervalo de confiança de 95% para o coeficiente da regressão usando Tauchen

ci_tauchen = reg_tauchen.coefCI

% Calculando o intervalo de confiança de 95% para o coeficiente da regressão usando Rouwenhorst

ci_rou = reg_rou.coefCI

% Repare que, nos dois casos, 0.95 está contido no intervalo.
%%
clear
%%% Item (5)

rho = 0.99;
sigma = 0.007;
mu = 0; %média do erro
m = 3;

%%% Item (1)

% Criando o grid de z
n = 9;
upper_bound = m*sigma*(sqrt(1-rho^2))^-1;
lower_bound = -upper_bound;

grid_z = linspace(lower_bound, upper_bound,n);

% Calculando a matriz de transição
delta_z = (upper_bound - lower_bound)*((n-1)^-1);
P = zeros(n,n);

% Calculando a primeira e última colunas da matriz
P(:,1) = normcdf((grid_z(1)-rho*grid_z+delta_z/2)/sigma);
P(:,n) = 1 - normcdf((grid_z(n)-rho*grid_z-delta_z/2)/sigma);

% Calculando as outras entradas da matriz
for j=2:(n-1)
    P(:,j)= normcdf((grid_z(j)+(delta_z/2)-rho*grid_z)/sigma) - normcdf((grid_z(j)-(delta_z/2)-rho*grid_z)/sigma);
end

% Matriz de transição
P;

%%% Item (2)

% Definindo variáveis que serão usadas
sigma_z = sqrt(sigma^2*((1-rho^2)^-1));
upper_bound_new = sigma_z*sqrt(n-1);
lower_bound_new = -upper_bound_new;
grid_z_rou = linspace(lower_bound_new, upper_bound_new,n); % rou é referência a Rouwenhorst

p = (1+rho)/2; % variável usada para computar a matriz P de forma recursiva

P_rou = [p 1-p; 1-p p]; %matriz P2
%P_rou3 = p*[P_rou zeros(2,1); zeros(1,2) 0] + (1-p)*[zeros(2,1) P_rou; zeros(1,2) 0] + (1-p)*[zeros(1,2) 0; P_rou zeros(2,1)] + p*[zeros(1,2) 0; zeros(2,1) P_rou]

for i=3:n
    P_rou = p*[P_rou zeros(i-1,1); zeros(1,i-1) 0] + (1-p)*[zeros(i-1,1) P_rou; zeros(1,i-1) 0] + ...
        (1-p)*[zeros(1,i-1) 0; P_rou zeros(i-1,1)] + p*[zeros(1,i-1) 0; zeros(i-1,1) P_rou];  
end

% Para normalizar as linhas:
for x=1:n
    P_rou(x,:) = P_rou(x,:)/sum(P_rou(x,:));
end

% Matriz de transição
P_rou;
    
%%% Item (3)

rng(1); % botando um seed para padronizar o processo gerador de dados aleatórios
t = 10000; % número de períodos
choque = normrnd(mu,sigma,1,t); %gerando choques aleatórios normalmente distribuídos (média 0, desvpad 0.007)

z_t = zeros(1,t);

for i=2:t
    z_t(i) = rho*z_t(i-1)+choque(i);
end

% Processo AR(1)
z_t;

% Gráfico AR(1)
plot(z_t,'b')
title ('Processo AR(1) simulado para 10000 períodos', 'FontSize', 34)
xlabel('Períodos', 'FontSize', 28)
ylabel('z_t', 'FontSize', 28)

% Simulando o processo discretizado via método de Tauchen

% Pegando o ponto no grid_z em que z = 0 (i.e., a mediana)

mediana = median(grid_z, 'all');

% Encontrando o estado inicial

estado = find(grid_z==mediana);

% Criando um vetor para armazenar os estados em que estou(começando da mediana = 0)
z_discretizado = [mediana zeros(1,t-1)];

for i=2:t
    cdf_estado = P(estado,1); %considerando que estou no estado inicial igual a 5, qual probabilidade de eu ir para outro estado
    for j = 1:n
        if normcdf(choque(i)/sigma) < cdf_estado
           estado = j;
           z_discretizado(i) = grid_z(estado); % vetor que acumula para qual estado eu fui dado o choque recebido (z_t discretizado)
           break
        else 
            cdf_estado = P(estado,j+1) + cdf_estado;
        end  
    end    
end    

% Processo discretizado (Tauchen)
z_discretizado;

% Criando o gráfico do proceso discretizado (Tauchen)
plot(z_discretizado,'r')
title ('Processo discretizado via Tauchen', 'FontSize', 34)
xlabel('Períodos', 'FontSize', 28)
ylabel('z discretizado', 'FontSize', 28)

% Plotando a simulação do AR(1) e do processo discretizado via Tauchen juntos
plot(z_t,'b')
hold on
title ('AR(1) com processo discretizado via Tauchen', 'FontSize', 34)
plot(z_discretizado,'r')
xlabel('Períodos', 'FontSize', 28)
ylabel('z', 'FontSize', 28)
hold off
legend('AR(1)', 'Discretizado (Tauchen)')

% Simulando o processo discretizado via método de Rouwenhorst

mediana_rou = median(grid_z_rou, 'all');
estado_rou = find(grid_z_rou == mediana_rou);

% Criando um vetor para armazenar os estados em que estou(começando da mediana = 0)
z_discretizado_rou = [mediana_rou zeros(1,t-1)];

for i=2:t
    cdf_estado_rou = P_rou(estado_rou,1); %considerando que estou no estado inicial igual a 5, qual probabilidade de eu ir para outro estado
    for j = 1:n
        if normcdf(choque(i)/sigma) < cdf_estado_rou
           estado_rou = j;
           z_discretizado_rou(i) = grid_z_rou(estado_rou); % vetor que acumula para qual estado eu fui dado o choque recebido (z_t discretizado)
           break
        else 
            cdf_estado_rou = P_rou(estado_rou,j+1) + cdf_estado_rou;
        end  
    end    
end    

% Processo discretizado(Rouwenhorst)
z_discretizado_rou;

% Criando o gráfico do proceso discretizado (Tauchen)
plot(z_discretizado_rou,'r')
title ('Processo discretizado via Rouwenhorst', 'FontSize', 34)
xlabel('Períodos', 'FontSize', 28)
ylabel('z discretizado', 'FontSize', 28)

% Plotando a simulação do AR(1) e do processo discretizado via Rouwenhorst juntos
plot(z_t,'b')
hold on
title ('AR(1) com processo discretizado via Rouwenhorst', 'FontSize', 34)
plot(z_discretizado_rou,'r')
xlabel('Períodos', 'FontSize', 28)
ylabel('z', 'FontSize', 28)
hold off
legend('AR(1)', 'Discretizado (Rouwenhorst)')

% Repare que eles não estão muito próximos. 
% Se refizermos o exercício, considerando 101 pontos no grid_z, o processo discretizado fica muito próximo do AR(1).

%%% Item (4)

% Rodando a regressão do z_discretizado no seu próprio lag (método de Tauchen)

reg_tauchen = fitlm(lagmatrix(z_discretizado,1),z_discretizado,'Intercept',false);

% Rodando a regressão do z_discretizado no seu próprio lag (método de Rouwenhorst)

reg_rou = fitlm(lagmatrix(z_discretizado_rou,1),z_discretizado_rou,'Intercept',false);

% O coeficiente da regressão é muito próximo de 0.99 (rho), então estão bem próximos do processo gerador de dados real.
% Repare que usando o método de Rouwenhorst, a aproximação do processo real é melhor.

% Calculando o intervalo de confiança de 95% para o coeficiente da regressão usando Tauchen

ci_tauchen = reg_tauchen.coefCI;

% Calculando o intervalo de confiança de 95% para o coeficiente da regressão usando Rouwenhorst

ci_rou = reg_rou.coefCI;

% Repare que, no caso de Tauchen, 0.99 está fora do intervalo de confiança.
% Uma possível explicação para isso é que, conforme rho -> 1, o método de Tauchen não funciona muito bem.
% Já no caso de Rouwenhorst, 0.99 está contido no intervalo de confiança. 

% Nesse caso, conseguimos visualizar graficamente que utilizamos poucos pontos no grid.
% Uma alternativa é aumentar o grid. Realizo esse teste colocando 101 pontos no grid e processo discretizado se aproxima ainda mais do processo contínuo.
% Além disso, com mais pontos no grid, rho = 0.99 está dentro do intervalo
% de confiança de 95% usando o processo discretizado via Tauchen.