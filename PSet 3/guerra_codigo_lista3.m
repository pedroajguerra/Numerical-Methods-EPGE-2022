%%
clear
clc
%%% M�todos Num�ricos (prof. C�zar Santos)
%%% Lista 3 - Pedro Augusto Januzzi Guerra

% Item 1

% Par�metros utilizados:
rho = 0.95; % persist�ncia do choque de produtividade
sigma = 0.007; % desvio-padr�o do choque
m = 3; % scaling parameter
beta = 0.987; % fator de desconto
mu = 2; % coeficiente de avers�o ao risco
alfa = 1/3; % participa��o do capital na fun��o de produ��o
delta = 0.012; % deprecia��o do capital

% Come�arei a quest�o discretizando o choque via M�todo de Tauchen
% Como nada foi dito quanto ao tamanho dos grids, vou considerar an�logo ao que fizemos na lista anterior

% Criando o grid de z
n = 7;
upper_bound = m*sigma*(sqrt(1-rho^2))^-1;
lower_bound = -upper_bound;
grid_z = linspace(lower_bound, upper_bound,n);

% Calculando a matriz de transi��o
delta_z = (upper_bound - lower_bound)*((n-1)^-1);
P = zeros(n,n);

% Calculando a primeira e �ltima coluna da matriz
P(:,1) = normcdf((grid_z(1)-rho*grid_z+delta_z/2)/sigma);
P(:,n) = 1 - normcdf((grid_z(n)-rho*grid_z-delta_z/2)/sigma);

% Calculando as outras entradas da matriz
for j=2:(n-1)
    P(:,j)= normcdf((grid_z(j)+(delta_z/2)-rho*grid_z)/sigma) - normcdf((grid_z(j)-(delta_z/2)-rho*grid_z)/sigma);
end

% Matriz de transi��o
P;

% Tirando a exponencial do grid dos choques
grid_z = exp(grid_z);
 
% Encontrando o capital de estado estacion�rio pela f�rmula do item (2)
k_ss = (1/alfa*((1/beta) + delta - 1))^(1/(alfa - 1));

% Criando o grid do capital
n_k = 500;
lower_bound_k = 0.75*k_ss;
upper_bound_k = 1.25*k_ss;

grid_k = linspace(lower_bound_k, upper_bound_k, n_k);

%% Redefinindo os par�metros para colocar nas fun��es que criei

parametro.mu = mu;
parametro.rho = rho;
parametro.sigma = sigma;
parametro.m = m ;
parametro.beta = beta;
parametro.alfa = alfa;
parametro.delta = delta;
parametro.grid_k = grid_k;
parametro.grid_z = grid_z;
parametro.n = n;
parametro.n_k = n_k;
parametro.P = P;
%% POLIN�MIO DE CHEBYSHEV

% Criei um arquivo de fun��o para calcular o polin�mio de Chebyshev
% Nomeei a fun��o de "cheb" e os argumentos s�o o grau do polin�mio e o ponto que estamos calculando
% Agora, checarei se a fun��o est� correta plotando o gr�fico visto nos slides:

% Criando um grid de pontos para testar
grid_chev = linspace(-1,1, n_k);

% Criando uma matriz 4 x 500, em que cada linha representar� uma ordem do polin�mio de chebyshev
cb = zeros(5,n_k);

% Fazendo um loop para preencher a matriz acima
for j = 1:5
    n_dim = j - 1; % fiz isso para n�o dar erro de argumento incorreto quando for preencher a matriz cbv
    for i = 1:n_k
        cb(j,i) = cheb(n_dim,grid_chev(i));
    end
end
%% Gr�ficos - Chebyshev
plot(grid_chev, cb(1,:),'k-', 'LineWidth',2)
hold on
plot(grid_chev, cb(2,:), 'g-', 'LineWidth',1.5)
plot(grid_chev, cb(3,:), 'b-', 'LineWidth',1.5)
plot(grid_chev, cb(4,:), 'r-', 'LineWidth',1.5)
plot(grid_chev, cb(5,:),'c', 'LineWidth',1.5)
xlabel('x'), ylabel('Tn(x)')
legend('n=0','n=1','n=2','n=3','n=4')
hold off

%% Encontrando gamas �timos

% Chutes iniciais para d = 1
d_aux = 1;
gama_chute = ones(n, d_aux+1);

% Ajustes para o fsolve (sugest�o do Rafael)
tol_fsolve = 10^(-10);
options = optimoptions(@fsolve, 'FunctionTolerance', tol_fsolve, 'StepTolerance', tol_fsolve, 'Maxiterations',100,'OptimalityTolerance', tol_fsolve);
%options = optimset('Display','off');

% Fazendo o for sugerido no item (5) - Ana Paula
d = 5; % poderia usar 4 porque, com 5, a 6a coluna est� zerando j� 
tic
for i = 1:d
    if i == 1
       zero_funcao = @(x) sistema_residuo(x, i, parametro);
       gama_otimo = fsolve(zero_funcao, gama_chute, options);       
    else
       for j = 1:n
          gama_otimo(j, i+1)= 0; % indo do d=1 para o d=2,3,4,...
       end
       zero_funcao = @(x) sistema_residuo(x, i, parametro);
       gama_otimo = fsolve(zero_funcao,gama_otimo, options);
    end
end
demora = toc

%% Fun��o Pol�tica para o Consumo no Grid K

% Definindo uma matriz que ser� preenchida com as escolhas �timas de c
pol_c = zeros(n,n_k);

for j = 1:n
    for i = 1:n_k
       pol_c(j,i) = consumo(gama_otimo(j,:),d,parametro,grid_k(i));
    end
end
pol_c = pol_c';

%% Encontrando a fun��o pol�tica do capital

z_matriz = repmat(grid_z, n_k, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz';

g = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - pol_c;

%% Fun��o Valor

% Definindo valores iniciais para colocar no loop
dif = 1;
iter = 0;

% Definindo a toler�ncia para fazer a converg�ncia da fun��o valor
tol = 10^-5;

% Chute inicial para a fun��o Valor
V_old = ((pol_c.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old = zeros(n_k,n);

% Definindo uma matriz que vai acumular as atualiza��es da fun��o valor
V_nova = V_old;

%% Itera��o Fun��o Valor
tic
while (dif > tol)
    V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela fun��o valor que o for originou
    for j = 1:n 
        for i = 1:n_k 
                % g que encontrei j� � minha fun��o pol�tica do capital
                % criando um "novo" grid ("aproximado") de capital como a dist�ncia de cada elemento at� o k_prime escolhido
                grid_knovo = abs(grid_k - g(i,j)); 
                % o k_prime escolhido � tal que:
                pol_k = min(grid_knovo);
                % encontrando a posi��o de k_prime no grid
                posicao = find(pol_k == grid_knovo);
                % calculando o consumo
                c = grid_z(j)*(grid_k(i)^alfa) + (1 - delta)*grid_k(i) - grid_k(posicao); 
                % calculando a utilidade dado o consumo definido acima
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; % fun��o utilidade CRRA
                else
                    u = -Inf; % condi��o para impedir consumo negativo
                end
               % Criando mais um for para calcular o valor esperado que aparece no final do lado direito da Eq. de Bellman 
               valor_esperado = 0;
               for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado = valor_esperado + P(j,s)*V_old(posicao,s);
               end
               vetor_valor = u + beta*valor_esperado;                 
        end
        V_nova(i,j) = vetor_valor;
        %V_nova(i,j) = u + beta*valor_esperado;
    end    
    dif = norm(V_nova - V_old); % tirando a norma
    iter = iter + 1; % contando as itera��es
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif]) % para printar o n�mero da itera��o e o a diferen�a entre V_nova e V_old
end
demora = toc % tempo para rodar o while

%% Erros da Equa��o de Euler

% Encontrando Erros da Equa��o de Euler

% Criando uma matriz para acumular os EEE
erro_euler = zeros(n_k,n);

for j = 1:n
    for i = 1:n_k
        c = pol_c(i,j);
        k_pol = g(i,j);
        parcela_1 = zeros(1,n); % primeira parcela dentro do valor esperado da minha equa��o de Euler 
        parcela_2 = zeros(1,n); % segunda parcela dentro do valor esperado da minha equa��o de Euler
        for a = 1:n
            parcela_1(a) = 1 - delta + alfa*grid_z(a)*(g(i,j)^(alfa - 1));
            c_pol = consumo(gama_otimo(a,:), d, parametro, k_pol);
            parcela_2(a) = c_pol^(-mu);
        end 
        argumento_esperanca = parcela_1.*parcela_2; % argumento que aparece dentro do operador esperan�a;
        esperanca = P(j,:)*argumento_esperanca'; % tirando o valor esperado
        u_inversa = (beta*esperanca)^(-1/mu); % multiplicando beta pelo valor esperado e tirando a inversa da fun��o utilidade
        arg = abs(1 - u_inversa/c); % tirando o m�dulo
        erro_euler(i,j) = log10(arg); % erro da equa��o de euler
    end    
end

% Calculando o m�ximo erro de Euler
m1 = max(erro_euler);
%% Gr�ficos

plot(V_nova);
title('Fun��o Valor');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,V_nova);
title('Fun��o Valor');
xlabel('Estoque de Capital');
ylabel('Choques');
zlabel('Fun��o Valor');

plot(pol_c);
title('Fun��o Pol�tica Consumo');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,pol_c);
title('Fun��o Pol�tica Consumo');
xlabel('Estoque de Capital');
ylabel('Choques');
zlabel('Consumo');

plot(g);
title('Fun��o Pol�tica Capital');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,g);
title('Fun��o Pol�tica Capital');
xlabel('Estoque de Capital');
ylabel('Choques');
zlabel('Capital');

plot(erro_euler);
title('Erros da Equa��o de Euler');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,erro_euler);
title('Erros da Equa��o de Euler');
xlabel('Capital');
ylabel('Choques');
zlabel('EEE');

%% Quest�o 2 - Colocacao

% Definindo quantidade de pontos a ser utilizada
ptos = 11;

% Criando "novo grid de capital"
grid_k_fea = linspace(lower_bound_k, upper_bound_k, ptos);

% Colocando na estrutura de parametro
parametro.ptos = ptos;
parametro.grid_k_fea = grid_k_fea;

% Definindo chute inicial para a_0
a = ones(n, ptos);

for i = 1:n
    for j = 1:ptos
        a(i,j) = j;
    end
end

zero_funcao_fea = @(x) sistema_residuo_fea(x, parametro);
a_otimo = fsolve(zero_funcao_fea, a);

%% Fun��o pol�tica para o consumo

% Definindo uma matriz que ser� preenchida com as escolhas �timas de c
pol_c_fea = zeros(n,n_k);

for j = 1:n
    for i = 1:n_k
       pol_c_fea(j,i) = consumo_colfea(grid_k(i),a_otimo(j,:),parametro);
    end
end
pol_c_fea = pol_c_fea';
%% Encontrando a fun��o pol�tica do capital

z_matriz = repmat(grid_z, n_k, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz';

g_fea = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - pol_c_fea;
%% Fun��o Valor

% Definindo valores iniciais para colocar no loop
dif_fea = 1;
iter_fea = 0;

% Definindo a toler�ncia para fazer a converg�ncia da fun��o valor
tol_fea = 10^-5;

% Chute inicial para a fun��o Valor
V_old_fea = ((pol_c_fea.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old_fea = zeros(n_k,n);

% Definindo uma matriz que vai acumular as atualiza��es da fun��o valor
V_nova_fea = V_old_fea;

%% Itera��o Fun��o Valor
tic
while (dif_fea > tol_fea)
    V_old_fea = V_nova_fea; % depois que os for's abaixo rodarem, atualizo a minha V_old pela fun��o valor que o for originou
    for j = 1:n 
        for i = 1:n_k 
                % g que encontrei j� � minha fun��o pol�tica do capital
                % criando um "novo" grid ("aproximado") de capital como a dist�ncia de cada elemento at� o k_prime escolhido
                grid_knovo_fea = abs(grid_k - g_fea(i,j)); 
                % o k_prime escolhido � tal que:
                pol_k_fea = min(grid_knovo_fea);
                % encontrando a posi��o de k_prime no grid
                posicao_fea = find(pol_k_fea == grid_knovo_fea);
                % calculando o consumo
                c_fea = grid_z(j)*(grid_k(i)^alfa) + (1 - delta)*grid_k(i) - grid_k(posicao_fea); 
                % calculando a utilidade dado o consumo definido acima
                if c_fea > 0
                    u_fea = (c_fea^(1-mu)-1)/(1-mu) ; % fun��o utilidade CRRA
                else
                    u_fea = -Inf; % condi��o para impedir consumo negativo
                end
               % Criando mais um for para calcular o valor esperado que aparece no final do lado direito da Eq. de Bellman 
               valor_esperado_fea = 0;
               for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado_fea = valor_esperado_fea + P(j,s)*V_old_fea(posicao_fea,s);
               end
               vetor_valor_fea = u_fea + beta*valor_esperado_fea;                 
        end
        V_nova_fea(i,j) = vetor_valor_fea;
        %V_nova_fea(i,j) = u + beta*valor_esperado;
    end    
    dif_fea = norm(V_nova_fea - V_old_fea); % tirando a norma
    iter_fea = iter_fea + 1; % contando as itera��es
    fprintf('Itera��es %4i %6.8f\n ', [iter_fea, dif_fea]) % para printar o n�mero da itera��o e o a diferen�a entre V_nova e V_old
end
demora_fea = toc % tempo para rodar o while

%% Erros da Equa��o de Euler

% Encontrando Erros da Equa��o de Euler

% Criando uma matriz para acumular os EEE
erro_euler_fea = zeros(n_k,n);

for j = 1:n
    for i = 1:n_k
        c_fea = pol_c_fea(i,j);
        k_pol_fea = g_fea(i,j);
        parcela_1_fea = zeros(1,n); % primeira parcela dentro do valor esperado da minha equa��o de Euler 
        parcela_2_fea = zeros(1,n); % segunda parcela dentro do valor esperado da minha equa��o de Euler
        for a = 1:n
            parcela_1_fea(a) = 1 - delta + alfa*grid_z(a)*(g_fea(i,j)^(alfa - 1));
            c_pol_fea = consumo_colfea(k_pol_fea,a_otimo(a,:),parametro);
            parcela_2_fea(a) = c_pol_fea^(-mu);
        end 
        argumento_esperanca_fea = parcela_1_fea.*parcela_2_fea; % argumento que aparece dentro do operador esperan�a;
        esperanca_fea = P(j,:)*argumento_esperanca_fea'; % tirando o valor esperado
        u_inversa_fea = (beta*esperanca_fea)^(-1/mu); % multiplicando beta pelo valor esperado e tirando a inversa da fun��o utilidade
        arg_fea = abs(1 - u_inversa_fea/c_fea); % tirando o m�dulo
        erro_euler_fea(i,j) = log10(arg_fea); % erro da equa��o de euler
    end    
end

% Calculando o maximo erro de Euler

m2 = max(erro_euler_fea);

%% Gr�ficos

plot(V_nova_fea);
title('Fun��o Valor');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,V_nova_fea);
title('Fun��o Valor');
xlabel('Estoque de Capital');
ylabel('Choques');
zlabel('Fun��o Valor');

plot(pol_c_fea);
title('Fun��o Pol�tica Consumo');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,pol_c_fea);
title('Fun��o Pol�tica Consumo');
xlabel('Estoque de Capital');
ylabel('Choques');
zlabel('Consumo');

plot(g_fea);
title('Fun��o Pol�tica Capital');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,g_fea);
title('Fun��o Pol�tica Capital');
xlabel('Estoque de Capital');
ylabel('Choques');
zlabel('Capital');

plot(erro_euler_fea);
title('Erros da Equa��o de Euler');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,erro_euler_fea);
title('Erros da Equa��o de Euler');
xlabel('Capital');
ylabel('Choques');
zlabel('EEE');

%% Quest�o 2 - Galerkin

% Igual a secao anterior::::
% Definindo quantidade de pontos a ser utilizada
ptos = 11;

% Criando "novo" grid de capital
grid_k_fea = linspace(lower_bound_k, upper_bound_k, ptos);

% Numero de intervalos para "integral"
int = 11;

% Colocando na estrutura de parametro
parametro.ptos = ptos;
parametro.grid_k_fea = grid_k_fea;
parametro.int = int;

% Como a estrutura das funcoes psi, consumo e residuo sao as mesmas, nao
% criarei novas funcoes para elas

% Definindo chute inicial para a_0
a = ones(n, ptos);

for i = 1:n
   a(i,:) = linspace(2,4,ptos); 
end

%for i = 1:n
 %   for j = 1:ptos
  %      a(i,j) = j;
   % end
%end

zero_funcao_fea_gal = @(x) sistema_residuo_fea_gal(x, parametro);
a_otimo_gal = fsolve(zero_funcao_fea_gal, a);

%% Fun��o pol�tica para o consumo

% Definindo uma matriz que ser� preenchida com as escolhas �timas de c
pol_c_fea_gal = zeros(n,n_k);

for j = 1:n
    for i = 1:n_k
       pol_c_fea_gal(j,i) = consumo_colfea(grid_k(i),a_otimo_gal(j,:),parametro);
    end
end
pol_c_fea_gal = pol_c_fea_gal';
%% Encontrando a fun��o pol�tica do capital

z_matriz = repmat(grid_z, n_k, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz';

g_fea_gal = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - pol_c_fea_gal;
%% Fun��o Valor

% Definindo valores iniciais para colocar no loop
dif_fea_gal = 1;
iter_fea_gal = 0;

% Definindo a toler�ncia para fazer a converg�ncia da fun��o valor
tol_fea_gal = 10^-5;

% Chute inicial para a fun��o Valor
V_old_fea_gal = ((pol_c_fea_gal.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old_fea_gal = zeros(n_k,n);

% Definindo uma matriz que vai acumular as atualiza��es da fun��o valor
V_nova_fea_gal = V_old_fea_gal;

%% Itera��o Fun��o Valor
tic
while (dif_fea_gal > tol_fea_gal)
    V_old_fea_gal = V_nova_fea_gal; % depois que os for's abaixo rodarem, atualizo a minha V_old pela fun��o valor que o for originou
    for j = 1:n 
        for i = 1:n_k 
                % g que encontrei j� � minha fun��o pol�tica do capital
                % criando um "novo" grid ("aproximado") de capital como a dist�ncia de cada elemento at� o k_prime escolhido
                grid_knovo_fea_gal = abs(grid_k - g_fea_gal(i,j)); 
                % o k_prime escolhido � tal que:
                pol_k_fea_gal = min(grid_knovo_fea_gal);
                % encontrando a posi��o de k_prime no grid
                posicao_fea_gal = find(pol_k_fea_gal == grid_knovo_fea_gal);
                % calculando o consumo
                c_fea_gal = grid_z(j)*(grid_k(i)^alfa) + (1 - delta)*grid_k(i) - grid_k(posicao_fea_gal); 
                % calculando a utilidade dado o consumo definido acima
                if c_fea_gal > 0
                    u_fea_gal = (c_fea_gal^(1-mu)-1)/(1-mu) ; % fun��o utilidade CRRA
                else
                    u_fea_gal = -Inf; % condi��o para impedir consumo negativo
                end
               % Criando mais um for para calcular o valor esperado que aparece no final do lado direito da Eq. de Bellman 
               valor_esperado_fea_gal = 0;
               for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado_fea_gal = valor_esperado_fea_gal + P(j,s)*V_old_fea_gal(posicao_fea_gal,s);
               end
               vetor_valor_fea_gal = u_fea_gal + beta*valor_esperado_fea_gal;                 
        end
        V_nova_fea_gal(i,j) = vetor_valor_fea_gal;
    end    
    dif_fea_gal = norm(V_nova_fea_gal - V_old_fea_gal); % tirando a norma
    iter_fea_gal = iter_fea_gal + 1; % contando as itera��es
    fprintf('Itera��es %4i %6.8f\n ', [iter_fea_gal, dif_fea_gal]) % para printar o n�mero da itera��o e o a diferen�a entre V_nova e V_old
end
demora_fea_gal = toc % tempo para rodar o while

%% Erros da Equa��o de Euler

% Encontrando Erros da Equa��o de Euler

% Criando uma matriz para acumular os EEE
erro_euler_fea_gal = zeros(n_k,n);

for j = 1:n
    for i = 1:n_k
        c_fea_gal = pol_c_fea_gal(i,j);
        k_pol_fea_gal = g_fea_gal(i,j);
        parcela_1_fea_gal = zeros(1,n); % primeira parcela dentro do valor esperado da minha equa��o de Euler 
        parcela_2_fea_gal = zeros(1,n); % segunda parcela dentro do valor esperado da minha equa��o de Euler
        for a = 1:n
            parcela_1_fea_gal(a) = 1 - delta + alfa*grid_z(a)*(g_fea_gal(i,j)^(alfa - 1));
            c_pol_fea_gal = consumo_colfea(k_pol_fea_gal,a_otimo_gal(a,:),parametro);
            parcela_2_fea_gal(a) = c_pol_fea_gal^(-mu);
        end 
        argumento_esperanca_fea_gal = parcela_1_fea_gal.*parcela_2_fea_gal; % argumento que aparece dentro do operador esperan�a;
        esperanca_fea_gal = P(j,:)*argumento_esperanca_fea_gal'; % tirando o valor esperado
        u_inversa_fea_gal = (beta*esperanca_fea_gal)^(-1/mu); % multiplicando beta pelo valor esperado e tirando a inversa da fun��o utilidade
        arg_fea_gal = abs(1 - u_inversa_fea_gal/c_fea_gal); % tirando o m�dulo
        erro_euler_fea_gal(i,j) = log10(arg_fea_gal); % erro da equa��o de euler
    end    
end

% Calculando o maximo erro de Euler
m3 = max(erro_euler_fea_gal);

%% Gr�ficos

plot(V_nova_fea_gal);
title('Fun��o Valor');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,V_nova_fea_gal);
title('Fun��o Valor');
xlabel('Estoque de Capital');
ylabel('Choques');
zlabel('Fun��o Valor');

plot(pol_c_fea_gal);
title('Fun��o Pol�tica Consumo');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,pol_c_fea_gal);
title('Fun��o Pol�tica Consumo');
xlabel('Estoque de Capital');
ylabel('Choques');
zlabel('Consumo');

plot(g_fea_gal);
title('Fun��o Pol�tica Capital');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,g_fea_gal);
title('Fun��o Pol�tica Capital');
xlabel('Estoque de Capital');
ylabel('Choques');
zlabel('Capital');

plot(erro_euler_fea_gal);
title('Erros da Equa��o de Euler');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,erro_euler_fea_gal);
title('Erros da Equa��o de Euler');
xlabel('Capital');
ylabel('Choques');
zlabel('EEE');









