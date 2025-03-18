%%
clear
clc
%%% Métodos Numéricos (prof. Cézar Santos)
%%% Lista 2 - Pedro Augusto Januzzi Guerra

% % Parâmetros utilizados:
% rho = 0.95; % persistência do choque de produtividade
% sigma = 0.007; % desvio-padrão do choque
% m = 3; % scaling parameter
% beta = 0.987; % fator de desconto
% mu = 2; % coeficiente de aversão ao risco
% alfa = 1/3; % participação do capital na função de produção
% delta = 0.012; % depreciação do capital

rho = 0.95; % persistência do choque de produtividade
sigma = 0.01; % desvio-padrão do choque
m = 3; % scaling parameter
beta = 0.95; % fator de desconto
mu = 2; % coeficiente de aversão ao risco
alfa = 1/3; % participação do capital na função de produção
delta = 0.05; % depreciação do capital

%%% Item (3)

% Começarei a questão discretizando o choque via Método de Tauchen
% Criando o grid de z
n = 7;
upper_bound = m*sigma*(sqrt(1-rho^2))^-1;
lower_bound = -upper_bound;
grid_z = linspace(lower_bound, upper_bound,n);

% Calculando a matriz de transição
delta_z = (upper_bound - lower_bound)*((n-1)^-1);
P = zeros(n,n);

% Calculando a primeira e última coluna da matriz
P(:,1) = normcdf((grid_z(1)-rho*grid_z+delta_z/2)/sigma);
P(:,n) = 1 - normcdf((grid_z(n)-rho*grid_z-delta_z/2)/sigma);

% Calculando as outras entradas da matriz
for j=2:(n-1)
    P(:,j)= normcdf((grid_z(j)+(delta_z/2)-rho*grid_z)/sigma) - normcdf((grid_z(j)-(delta_z/2)-rho*grid_z)/sigma);
end

% Matriz de transição
P;

% Tirando a exponencial do grid dos choques
grid_z = exp(grid_z);
 
% Encontrando o capital de estado estacionário pela fórmula do item (2)
k_ss = (1/alfa*((1/beta) + delta - 1))^(1/(alfa - 1));

% Criando o grid do capital
n_k = 500;
lower_bound_k = 0.75*k_ss;
upper_bound_k = 1.25*k_ss;

grid_k = linspace(lower_bound_k, upper_bound_k, n_k);

% Note que, para cada (k,z), teremos um k' que maximiza V(k,z).
% Logo, nossa função valor será uma matriz de n_k por n
% Resolverei, inicialmente, por força bruta:

% Definindo um chute inicial para a função valor com base na sugestão de notas de aula de Ben Moll e sugestão da Ana Paula (monitora)
for i = 1:n
    for j = 1:n_k
       c(j,i) = grid_z(i)*(grid_k(j)^alfa) + (1-delta)*grid_k(j) - grid_k(j); 
    end
end
V_old = ((c.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old = zeros(n_k,n); % Testei um chute inicial V(.)=0, porém é um chute ruim e leva mais tempo para rodar.

% Definindo uma matriz que vai acumular as atualizações da função valor
V_nova = V_old;

% Definindo uma matriz de zeros para a policy function, g.
g = zeros(n_k,n);

% Definindo valores iniciais para colocar no loop
dif = 1;
iter = 0;

% Definindo a tolerância para fazer a convergência da função valor
tol = 10^-5;

%% Força Bruta
tic
while (dif > tol)
    V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela função valor que o for originou
    for j = 1:n % selecionando choque a choque
        for i = 1:n_k % dado o choque escolhido, percorrendo todos os possíveis valores de k
            for a = 1:n_k % dado o par (choque, capital), percorrendo novamente o grid_k para selecionar o k'
                c = grid_z(j)*(grid_k(i)^alfa) + (1-delta)*grid_k(i) - grid_k(a); % restrição que define o consumo, em que k(a) seria o k' na notação que estamos habituados
                % calculando a utilidade dado o consumo definido acima
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; % função utilidade CRRA
                else
                    u = -Inf; % condição para impedir consumo negativo
                end
               % Criando mais um for para calcular o valor esperado que aparece no final do lado direito da Eq. de Bellman 
               valor_esperado = 0;
               for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado; %vetor que me gera todas as funções valores aplicadas nos k'                
               %vetor_valor(a) = u + beta*dot(P(j,:),V_old(a,:)); % alternativa sem loop, mas demorou muito
            end
            % V_nova(i,j) = max(vetor_valor); % pegando a maior função valor aplicada nos k' e tornando essa minha V_nova
            % indice = find(vetor_valor == V_nova(i,j)); % encontrando o índice da função valor que gerou minha V_nova
            [V_nova(i,j), indice] = max(vetor_valor); % tentei isso, mas não funcionou. Se sobrar tempo, tentar resolver assim também.
            % V_nova(i,j)=V_n;
            g(i,j) = grid_k(indice); % função política de capital (k').
        end    
    end 
    % dif = norm(V_nova - V_old); % tirando a norma
    dif = max(max(abs(V_nova-V_old)));
    iter = iter + 1; % contando as iterações
    fprintf('Iterações %4i %6.8f\n ', [iter, dif]) % para printar o número da iteração e o a diferença entre V_nova e V_old
end
demora = toc % tempo para rodar o while

%% CONSUMO 
% Resta encontrar a função política para o consumo
z_matriz = repmat(grid_z, n_k, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz'; %Transpondo k_matriz para se tornar 500x7

% Encontrando a função política para o consumo, c
c = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - g;

%% MONOTONICIDADE
% Acima, usei força bruta. Vamos testar usando monotonicidade para tornar o código mais eficiente.
% A ideia é simples: suponha que cada agente "mora" em um ponto no grid de
% capital. É razoável assumir que as decisões de poupança de agentes mais
% ricos terão magnitudes maiores que aquelas tomadas por agentes mais
% pobres. Assim, a ideia da monotonicidade é: dado que um agente X escolheu
% um nível de k' que mora no grid_k, o agente X+1 escolherá um k'(y) maior
% ou igual ao k'(x).

% As explicações linha a linha são as mesmas que no caso força bruta
tic
while (dif > tol)
    V_old = V_nova; 
    for j = 1:n 
        for i = 1:n_k 
            if i==1
                g_inicial=1; % definindo um ponto inicial para procurar a função política
            end      
            for a = g_inicial:n_k 
                c = grid_z(j)*(grid_k(i)^alfa) + (1-delta)*grid_k(i) - grid_k(a); 
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; 
                else
                    u = -Inf; 
                end
               valor_esperado = 0;
               for s = 1:n 
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado;                 
            end
            V_nova(i,j) = max(vetor_valor); 
            indice = find(vetor_valor == V_nova(i,j));
            g(i,j) = grid_k(indice);
            g_inicial = indice; % dado que k_1 > k_2, então g(k_1) > g(k_2). Assim, atualizo o ponto no grid_k que vou começar a procurar o k' maximizador
        end    
    end 
    dif = norm(V_nova - V_old);
    iter = iter + 1;
    fprintf('Iterações %4i %6.8f\n ', [iter, dif])
end
demora = toc

%%
% Resta encontrar a função política para o consumo
z_matriz = repmat(grid_z, n_k, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz'; %Transpondo k_matriz para se tornar 500x7

% Encontrando a função política para o consumo
c = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - g;

%% Gráficos

plot(V_nova);
title('Função Valor');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,V_nova);

plot(c);
title('Função Política Consumo');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,c);

plot(g);
title('Função Política Capital');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,g);
%% ERROS DE EQUAÇÃO DE EULER

% Encontrando Erros da Equação de Euler

% Criando uma matriz para acumular os EEE

erro_euler = zeros(n_k,n);

for j = 1:n
    for i = 1:n_k
        parcela_1 = zeros(1,n); % primeira parcela dentro do valor esperado da minha equação de Euler 
        parcela_2 = zeros(1,n); % segunda parcela dentro do valor esperado da minha equação de Euler
        for a = 1:n
            parcela_1(a) = 1 - delta + alfa*grid_z(a)*(g(i,j)^(alfa - 1));
            idx = find(g(i,j) == grid_k);
            parcela_2(a) = c(idx,a)^(-mu);
        end 
        argumento_esperanca = parcela_1.*parcela_2; % argumento que aparece dentro do operador esperança;
        esperanca = P(j,:)*argumento_esperanca'; % tirando o valor esperado
        u_inversa = (beta*esperanca)^(-1/mu); % multiplicando beta pelo valor esperado e tirando a inversa da função utilidade
        arg = abs(1 - u_inversa/c(i,j)); % tirando o módulo
        erro_euler(i,j) = log10(arg); % erro da equação de euler
    end    
end

%% Gráficos - EEE

plot(erro_euler);
title('Erro de Euler');
xlabel('Estoque de Capital');
mesh(erro_euler);
mesh(k_matriz, z_matriz,erro_euler);
%% ACELERADOR COM FORÇA BRUTA - ITEM (4)

% Vou utilizar o acelerador
tic

while (dif > tol)
    if mod(iter, 10) == 0 || iter <= 100 % testar com limites maximos de iteração diferentes e com "&&" 
        % a função mod(iter,10) vai pegar o número da última iteração e
        % dividir por 10. Se o resto for igual a 0 (divisão exata), calculo
        % o maximo. Caso contrário, não, mas coloquei p isso acontecer só
        % depois da centésima iteração.
        V_old = V_nova; 
        for j = 1:n 
            for i = 1:n_k 
                for a = 1:n_k 
                    c = grid_z(j)*(grid_k(i)^alfa) + (1-delta)*grid_k(i) - grid_k(a);                     
                    if c > 0
                        u = (c^(1-mu)-1)/(1-mu);
                    else
                        u = -Inf; 
                    end
                    valor_esperado = 0;
                    for s = 1:n 
                        valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
                    end
                    vetor_valor(a) = u + beta*valor_esperado;                 
                end
                V_nova(i,j) = max(vetor_valor); 
                indice = find(vetor_valor == V_nova(i,j));
                g(i,j) = grid_k(indice); 
            end    
        end
    else
        V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela função valor que o for originou
        for j = 1:n % selecionando choque a choque
            for i = 1:n_k % dado o choque escolhido, percorrendo todos os possíveis valores de k
                c = grid_z(j)*(grid_k(i)^alfa) + (1-delta)*grid_k(i) - g(i,j);
                u = (c^(1-mu)-1)/(1-mu);
                valor_esperado = 0;
                for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado = valor_esperado + P(j,s)*V_old(find(g(i,j)==grid_k),s);
                end
                V_nova(i,j) = u + beta*valor_esperado;
            end
        end
    end
    dif = norm(V_nova - V_old);
    iter = iter + 1;
    fprintf('Iterações %4i %6.8f\n ', [iter, dif])
end

demora = toc 

%% ACELERADOR COM MONOTONICIDADE
tic

while (dif > tol)
    if mod(iter, 10) == 0 || iter <= 10 
        V_old = V_nova; 
        for j = 1:n 
            for i = 1:n_k 
                if i==1
                    g_inicial=1;
                end    
                for a = g_inicial:n_k 
                    c = grid_z(j)*(grid_k(i)^alfa) + (1-delta)*grid_k(i) - grid_k(a);                     
                    if c > 0
                        u = (c^(1-mu)-1)/(1-mu);
                    else
                        u = -Inf; 
                    end
                    valor_esperado = 0;
                    for s = 1:n 
                        valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
                    end
                    vetor_valor(a) = u + beta*valor_esperado;                 
                end
                V_nova(i,j) = max(vetor_valor); 
                indice = find(vetor_valor == V_nova(i,j));
                g(i,j) = grid_k(indice); 
                g_inicial = indice;
            end    
        end
    else
        V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela função valor que o for originou
        for j = 1:n % selecionando choque a choque
            for i = 1:n_k % dado o choque escolhido, percorrendo todos os possíveis valores de k
                c = grid_z(j)*(grid_k(i)^alfa) + (1-delta)*grid_k(i) - g(i,j);
                u = (c^(1-mu)-1)/(1-mu);
                valor_esperado = 0;
                for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado = valor_esperado + P(j,s)*V_old(find(g(i,j)==grid_k),s);
                end
                V_nova(i,j) = u + beta*valor_esperado;
            end
        end
    end
    dif = norm(V_nova - V_old);
    iter = iter + 1;
    fprintf('Iterações %4i %6.8f\n ', [iter, dif])
end

demora = toc 

%% CONSUMO
% Resta encontrar a função política para o consumo
z_matriz = repmat(grid_z, n_k, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz'; %Transpondo k_matriz para se tornar 500x7

% Encontrando a função política para o consumo, c
c = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - g
%% Gráficos

plot(V_nova);
title('Função Valor');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,V_nova);

plot(c);
title('Função Política Consumo');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,c);

plot(g);
title('Função Política Capital');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,g);

%% MULTIGRID - ITEM (5)

% A ideia do multigrid é iterar a função valor para um grid mais grosso a
% fim de arrumar um chute inicial melhor. Uma coisa importante de se notar
% é a quantidade de iterações que levamos para ter aa convergência.
% Utilizando força bruta e V_0 = 0, a iteração da função valor por força
% bruta e com 500 pontos no grid leva 1161 iterações e 407,62 segundos.
% Quando utilizo esse chute inicial para o grid mais grosso (100 pontos),
% leva 17 segundos para realizar a iteração. Utilizando esse resultado como
% chute inicial para o grid mais grosso (500 pontos), são necessárias 508
% iterações e 177 segundos, ou seja, no total levo 194 segundos, enquanto
% sem multigrid levo 407. Com 5000 pontos no grid, força bruta é muito
% lento. Monotonicidade vai ajudar, mas seria necessário acelerar ainda
% mais o código (utilizando concavidade + monotonicidade, por exemplo, mas
% nao sei se vai dar tempo de fazer isso).

clear
clc

% Parâmetros utilizados:
rho = 0.95; % persistência do choque de produtividade
sigma = 0.007; % desvio-padrão do choque
m = 3; % scaling parameter
beta = 0.987; % fator de desconto
mu = 2; % coeficiente de aversão ao risco
alfa = 1/3; % participação do capital na função de produção
delta = 0.012; % depreciação do capital

% Começarei a questão discretizando o choque via Método de Tauchen
% Criando o grid de z
n = 7;
upper_bound = m*sigma*(sqrt(1-rho^2))^-1;
lower_bound = -upper_bound;
grid_z = linspace(lower_bound, upper_bound,n);

% Calculando a matriz de transição
delta_z = (upper_bound - lower_bound)*((n-1)^-1);
P = zeros(n,n);

% Calculando a primeira e última coluna da matriz
P(:,1) = normcdf((grid_z(1)-rho*grid_z+delta_z/2)/sigma);
P(:,n) = 1 - normcdf((grid_z(n)-rho*grid_z-delta_z/2)/sigma);

% Calculando as outras entradas da matriz
for j=2:(n-1)
    P(:,j)= normcdf((grid_z(j)+(delta_z/2)-rho*grid_z)/sigma) - normcdf((grid_z(j)-(delta_z/2)-rho*grid_z)/sigma);
end

% Matriz de transição
P;

% Tirando a exponencial do grid dos choques
grid_z = exp(grid_z);
 
% Encontrando o capital de estado estacionário pela fórmula do item (2)
k_ss = (1/alfa*((1/beta) + delta - 1))^(1/(alfa - 1));

% Possíveis tamanhos para o grid do capital
n_k1 = 100;
n_k2 = 500;
n_k3 = 5000;

% Criando o grid do capital
lower_bound_k = 0.75*k_ss;
upper_bound_k = 1.25*k_ss;

grid_k1 = linspace(lower_bound_k, upper_bound_k, n_k1);

% Definindo um chute inicial para a função valor com base na sugestão de notas de aula de Ben Moll e sugestão da Ana Paula (monitora)
for i = 1:n
    for j = 1:n_k1
       c(j,i) = grid_z(i)*(grid_k1(j)^alfa) + (1-delta)*grid_k1(j) - grid_k1(j); 
    end
end
V_old = ((c.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old = zeros(n_k1,n);

% Definindo uma matriz que vai acumular as atualizações da função valor
V_nova = V_old;

% Definindo uma matriz de zeros para a policy function, g.
g = zeros(n_k1,n);

% Definindo valores iniciais para colocar no loop
dif = 1;
iter = 0;

% Definindo a tolerância para fazer a convergência da função valor
tol = 10^-5;

% Farei, incialmente, por força bruta
tic
while (dif > tol)
    V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela função valor que o for originou
    for j = 1:n % selecionando choque a choque
        for i = 1:n_k1 % dado o choque escolhido, percorrendo todos os possíveis valores de k
            for a = 1:n_k1 % dado o par (choque, capital), percorrendo novamente o grid_k para selecionar o k'
                c = grid_z(j)*(grid_k1(i)^alfa) + (1-delta)*grid_k1(i) - grid_k1(a); % restrição que define o consumo, em que k(a) seria o k' na notação que estamos habituados
                % calculando a utilidade dado o consumo definido acima
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; % função utilidade CRRA
                else
                    u = -Inf; % condição para impedir consumo negativo
                end
               % Criando mais um for para calcular o valor esperado que aparece no final do lado direito da Eq. de Bellman 
               valor_esperado = 0;
               for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado; %vetor que me gera todas as funções valores aplicadas nos k'                
               %vetor_valor(a) = u + beta*dot(P(j,:),V_old(a,:)); % alternativa sem loop, mas demorou muito
            end
            V_nova(i,j) = max(vetor_valor); % pegando a maior função valor aplicada nos k' e tornando essa minha V_nova
            indice = find(vetor_valor == V_nova(i,j)); % encontrando o índice da função valor que gerou minha V_nova
            %[V_nova, index] = max(vetor_valor) % tentei isso, mas não funcionou. Se sobrar tempo, tentar resolver assim também.
            g(i,j) = grid_k1(indice); % função política de capital (k').
        end    
    end 
    dif = norm(V_nova - V_old); % tirando a norma
    iter = iter + 1; % contando as iterações
    fprintf('Iterações %4i %6.8f\n ', [iter, dif]) % para printar o número da iteração e o a diferença entre V_nova e V_old
end
demora = toc % tempo para rodar o while

V_nova_k1 = V_nova; % apenas para armazenar os valores dessa iteração
g_k1 = g; % apenas para armazenar os valores dessa iteração

%% Próximo passo: grid com 500 pontos (ainda com força bruta)

% Criando o segundo grid de capital, agora com 500 pontos
grid_k2 = linspace(lower_bound_k, upper_bound_k, n_k2);

% Informações prévias que precisamos retomar
for i = 1:n
    for j = 1:n_k2
      c(j,i) = grid_z(i)*(grid_k2(j)^alfa) + (1-delta)*grid_k2(j) - grid_k2(j); 
    end
end
V_old = ((c.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old = zeros(n_k2,n);
g = zeros(n_k2,n);

% Interpolando utilizando spline
for j = 1:n
    V_old(:,j) = spline(grid_k1, V_nova_k1(:,j), grid_k2);
end

V_nova = V_old;

dif = 1;
iter = 0;
tol = 10^-5;

tic
while (dif > tol)
    V_old = V_nova; 
    for j = 1:n 
        for i = 1:n_k2 
            for a = 1:n_k2
                c = grid_z(j)*(grid_k2(i)^alfa) + (1-delta)*grid_k2(i) - grid_k2(a); 
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; 
                else
                    u = -Inf; 
                end
               valor_esperado = 0;
               for s = 1:n 
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado; 
            end
            V_nova(i,j) = max(vetor_valor); 
            indice = find(vetor_valor == V_nova(i,j)); 
            g(i,j) = grid_k2(indice); 
        end    
    end 
    dif = norm(V_nova - V_old); 
    iter = iter + 1; 
    fprintf('Iterações %4i %6.8f\n ', [iter, dif]) 
end
demora = toc 

V_nova_k2 = V_nova;
g_k3 = g;

%% Próximo passo: grid com 5000 pontos (ainda com força bruta)

% Criando o segundo grid de capital, agora com 5000 pontos
grid_k3 = linspace(lower_bound_k, upper_bound_k, n_k3);

% Informações prévias que precisamos retomar
for i = 1:n
    for j = 1:n_k3
       c(j,i) = grid_z(i)*(grid_k3(j)^alfa) + (1-delta)*grid_k3(j) - grid_k3(j); 
    end
end
V_old = ((c.^(1-mu)-1)./(1-mu))./(1-beta);
% V_old = zeros(n_k3,n)
g = zeros(n_k3,n);

% Interpolando
for j = 1:n
    V_old(:,j) = spline(grid_k2, V_nova_k2(:,j), grid_k3);
end

V_nova = V_old;

dif = 1;
iter = 0;
tol = 10^-5;

tic
while (dif > tol)
    V_old = V_nova; 
    for j = 1:n 
        for i = 1:n_k3 
            for a = 1:n_k3 
                c = grid_z(j)*(grid_k3(i)^alfa) + (1-delta)*grid_k3(i) - grid_k3(a); 
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; 
                else
                    u = -Inf; 
                end
               valor_esperado = 0;
               for s = 1:n 
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado; 
            end
            V_nova(i,j) = max(vetor_valor); 
            indice = find(vetor_valor == V_nova(i,j)); 
            g(i,j) = grid_k3(indice); 
        end    
    end 
    dif = norm(V_nova - V_old); 
    iter = iter + 1; 
    fprintf('Iterações %4i %6.8f\n ', [iter, dif]) 
end
demora = toc 

%% MULTIGRID COM MONOTONICIDADE

clear
clc

% Parâmetros utilizados:
rho = 0.95; % persistência do choque de produtividade
sigma = 0.007; % desvio-padrão do choque
m = 3; % scaling parameter
beta = 0.987; % fator de desconto
mu = 2; % coeficiente de aversão ao risco
alfa = 1/3; % participação do capital na função de produção
delta = 0.012; % depreciação do capital

% Começarei a questão discretizando o choque via Método de Tauchen
% Criando o grid de z
n = 7;
upper_bound = m*sigma*(sqrt(1-rho^2))^-1;
lower_bound = -upper_bound;
grid_z = linspace(lower_bound, upper_bound,n);

% Calculando a matriz de transição
delta_z = (upper_bound - lower_bound)*((n-1)^-1);
P = zeros(n,n);

% Calculando a primeira e última coluna da matriz
P(:,1) = normcdf((grid_z(1)-rho*grid_z+delta_z/2)/sigma);
P(:,n) = 1 - normcdf((grid_z(n)-rho*grid_z-delta_z/2)/sigma);

% Calculando as outras entradas da matriz
for j=2:(n-1)
    P(:,j)= normcdf((grid_z(j)+(delta_z/2)-rho*grid_z)/sigma) - normcdf((grid_z(j)-(delta_z/2)-rho*grid_z)/sigma);
end

% Matriz de transição
P;

% Tirando a exponencial do grid dos choques
grid_z = exp(grid_z);
 
% Encontrando o capital de estado estacionário pela fórmula do item (2)
k_ss = (1/alfa*((1/beta) + delta - 1))^(1/(alfa - 1));

% Possíveis tamanhos para o grid do capital
n_k1 = 100;
n_k2 = 500;
n_k3 = 5000;

% Criando o grid do capital
lower_bound_k = 0.75*k_ss;
upper_bound_k = 1.25*k_ss;

grid_k1 = linspace(lower_bound_k, upper_bound_k, n_k1);

% Definindo um chute inicial para a função valor com base na sugestão de notas de aula de Ben Moll e sugestão da Ana Paula (monitora)
for i = 1:n
    for j = 1:n_k1
       c(j,i) = grid_z(i)*(grid_k1(j)^alfa) + (1-delta)*grid_k1(j) - grid_k1(j); 
    end
end
V_old = ((c.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old = zeros(n_k1,n);

% Definindo uma matriz que vai acumular as atualizações da função valor
V_nova = V_old;

% Definindo uma matriz de zeros para a policy function, g.
g = zeros(n_k1,n);

% Definindo valores iniciais para colocar no loop
dif = 1;
iter = 0;

% Definindo a tolerância para fazer a convergência da função valor
tol = 10^-5;

% Farei, incialmente, por força bruta
tic
while (dif > tol)
    V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela função valor que o for originou
    for j = 1:n % selecionando choque a choque
        for i = 1:n_k1 % dado o choque escolhido, percorrendo todos os possíveis valores de k
            if i == 1 
                g_inicial = 1;
            end
            for a = g_inicial:n_k1 % dado o par (choque, capital), percorrendo novamente o grid_k para selecionar o k'
                c = grid_z(j)*(grid_k1(i)^alfa) + (1-delta)*grid_k1(i) - grid_k1(a); % restrição que define o consumo, em que k(a) seria o k' na notação que estamos habituados
                % calculando a utilidade dado o consumo definido acima
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; % função utilidade CRRA
                else
                    u = -Inf; % condição para impedir consumo negativo
                end
               % Criando mais um for para calcular o valor esperado que aparece no final do lado direito da Eq. de Bellman 
               valor_esperado = 0;
               for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado; %vetor que me gera todas as funções valores aplicadas nos k'                
               %vetor_valor(a) = u + beta*dot(P(j,:),V_old(a,:)); % alternativa sem loop, mas demorou muito
            end
            V_nova(i,j) = max(vetor_valor); % pegando a maior função valor aplicada nos k' e tornando essa minha V_nova
            indice = find(vetor_valor == V_nova(i,j)); % encontrando o índice da função valor que gerou minha V_nova
            %[V_nova, index] = max(vetor_valor) % tentei isso, mas não funcionou. Se sobrar tempo, tentar resolver assim também.
            g(i,j) = grid_k1(indice); % função política de capital (k').
            g_inicial = indice;
        end    
    end 
    dif = norm(V_nova - V_old); % tirando a norma
    iter = iter + 1; % contando as iterações
    fprintf('Iterações %4i %6.8f\n ', [iter, dif]) % para printar o número da iteração e o a diferença entre V_nova e V_old
end
demora = toc % tempo para rodar o while

V_nova_k1 = V_nova;
g_k1 = g;

%% Próximo passo: grid com 500 pontos (agora com monotonicidade)

% Criando o segundo grid de capital, agora com 500 pontos
grid_k2 = linspace(lower_bound_k, upper_bound_k, n_k2);

% Informações prévias que precisamos retomar
for i = 1:n
    for j = 1:n_k2
      c(j,i) = grid_z(i)*(grid_k2(j)^alfa) + (1-delta)*grid_k2(j) - grid_k2(j); 
    end
end
V_old = ((c.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old = zeros(n_k2,n);
g = zeros(n_k2,n);

% Interpolando utilizando spline
for j = 1:n
    V_old(:,j) = spline(grid_k1, V_nova_k1(:,j), grid_k2);
end

V_nova = V_old;

dif = 1;
iter = 0;
tol = 10^-5;

tic
while (dif > tol)
    V_old = V_nova; 
    for j = 1:n 
        for i = 1:n_k2
            if i == 1
                g_inicial = 1;
            end
            for a = g_inicial:n_k2
                c = grid_z(j)*(grid_k2(i)^alfa) + (1-delta)*grid_k2(i) - grid_k2(a); 
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; 
                else
                    u = -Inf; 
                end
               valor_esperado = 0;
               for s = 1:n 
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado; 
            end
            V_nova(i,j) = max(vetor_valor); 
            indice = find(vetor_valor == V_nova(i,j)); 
            g(i,j) = grid_k2(indice); 
            g_inicial = indice;
        end    
    end 
    dif = norm(V_nova - V_old); 
    iter = iter + 1; 
    fprintf('Iterações %4i %6.8f\n ', [iter, dif]) 
end
demora = toc 

V_nova_k2 = V_nova;
g_k2 = g;

%% Próximo passo: grid com 5000 pontos (agora com monotonicidade)
n_k3 = 5000;
% Criando o segundo grid de capital, agora com 5000 pontos
grid_k3 = linspace(lower_bound_k, upper_bound_k, n_k3);

% Informações prévias que precisamos retomar
for i = 1:n
    for j = 1:n_k3
       c(j,i) = grid_z(i)*(grid_k3(j)^alfa) + (1-delta)*grid_k3(j) - grid_k3(j); 
    end
end
V_old = ((c.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old = zeros(n_k3,n)
g = zeros(n_k3,n);

% Interpolando
for j = 1:n
    V_old(:,j) = spline(grid_k2, V_nova_k2(:,j), grid_k3);
end

V_nova = V_old;

dif = 1;
iter = 0;
tol = 10^-5;

tic
while (dif > tol)
    V_old = V_nova; 
    for j = 1:n 
        for i = 1:n_k3
            if i == 1
                g_inicial = 1;
            end
            for a = 1:n_k3 
                c = grid_z(j)*(grid_k3(i)^alfa) + (1-delta)*grid_k3(i) - grid_k3(a); 
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; 
                else
                    u = -Inf; 
                end
               valor_esperado = 0;
               for s = 1:n 
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado; 
            end
            V_nova(i,j) = max(vetor_valor); 
            indice = find(vetor_valor == V_nova(i,j)); 
            g(i,j) = grid_k3(indice); 
        end    
    end 
    dif = norm(V_nova - V_old); 
    iter = iter + 1; 
    fprintf('Iterações %4i %6.8f\n ', [iter, dif]) 
end
demora = toc 

V_nova_k3 = V_nova;
g_k3 = g;
%% Consumo
% Resta encontrar a função política para o consumo
z_matriz = repmat(grid_z, n_k3, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k3, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz'; %Transpondo k_matriz para se tornar 500x7

% Encontrando a função política para o consumo
c = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - g_k3;

%% ERROS DE EQUAÇÃO DE EULER

% Encontrando Euler Equation Errors

% Criando uma matriz para acumular os EEE

erro_euler_mult = zeros(n_k3,n);

for j = 1:n
    for i = 1:n_k3
        parcela_1 = zeros(1,n); % primeira parcela dentro do valor esperado da minha equação de Euler 
        parcela_2 = zeros(1,n); % segunda parcela dentro do valor esperado da minha equação de Euler
        for a = 1:n
            parcela_1(a) = 1 - delta + alfa*grid_z(a)*(g_k3(i,j)^(alfa - 1));
            idx = find(g_k3(i,j) == grid_k3);
            parcela_2(a) = c(idx,a)^(-mu);
        end 
        argumento_esperanca = parcela_1.*parcela_2;
        esperanca = P(j,:)*argumento_esperanca';
        u_inversa = (beta*esperanca)^(-1/mu);
        arg = abs(1 - u_inversa/c(i,j));
        erro_euler_mult(i,j) = log10(arg); 
    end    
end

%% Gráficos
z_matriz = repmat(grid_z, n_k3, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k3, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz';

plot(V_nova);
title('Função Valor - Multigrid');
xlabel('Estoque de Capital');

mesh(V_nova);

mesh(erro_euler_mult);

plot(c);
title('Função Política Consumo - Multigrid');
xlabel('Estoque de Capital');

mesh(c);

plot(g);
title('Função Política Capital - Multigrid');
xlabel('Estoque de Capital');

mesh(g);

plot(erro_euler_mult);
title('Erro de Euler - Multigrid');
xlabel('Estoque de Capital');

%surf(grid_k3, grid_z, V_nova,'FaceColor','interp');
%title('Função Valor - Multigrid');
%colormap(hot);

%surf(grid_k3, grid_z, c,'FaceColor','interp');
%title('Função Valor - Multigrid');
%colormap(hot);

%surf(grid_k3, grid_z, g,'FaceColor','interp');
%title('Função Valor - Multigrid');
%colormap(hot);

%surf(grid_k3, grid_z, erro_euler_mult,'FaceColor','interp');
%title('Função Valor - Multigrid');
%colormap(hot);
%% GRID ENDÓGENO - ITEM (6)

% Repetindo as informações iniciais como nos itens anteriores
clear
clc

% Parâmetros utilizados:
rho = 0.95; % persistência do choque de produtividade
sigma = 0.007; % desvio-padrão do choque
m = 3; % scaling parameter
beta = 0.987; % fator de desconto
mu = 2; % coeficiente de aversão ao risco
alfa = 1/3; % participação do capital na função de produção
delta = 0.012; % depreciação do capital

% Criando o grid de z
n = 7;
upper_bound = m*sigma*(sqrt(1-rho^2))^-1;
lower_bound = -upper_bound;
grid_z = linspace(lower_bound, upper_bound,n);

% Calculando a matriz de transição
delta_z = (upper_bound - lower_bound)*((n-1)^-1);
P = zeros(n,n);

% Calculando a primeira e última coluna da matriz
P(:,1) = normcdf((grid_z(1)-rho*grid_z+delta_z/2)/sigma);
P(:,n) = 1 - normcdf((grid_z(n)-rho*grid_z-delta_z/2)/sigma);

% Calculando as outras entradas da matriz
for j=2:(n-1)
    P(:,j)= normcdf((grid_z(j)+(delta_z/2)-rho*grid_z)/sigma) - normcdf((grid_z(j)-(delta_z/2)-rho*grid_z)/sigma);
end

% Matriz de transição
P;

% Tirando a exponencial do grid dos choques
grid_z = exp(grid_z);
 
% Encontrando o capital de estado estacionário pela fórmula do item (2)
k_ss = (1/alfa*((1/beta) + delta - 1))^(1/(alfa - 1));

% Criando o grid do capital (não é o endógeno, obviamente)
n_k = 500;
lower_bound_k = 0.75*k_ss;
upper_bound_k = 1.25*k_ss;

grid_k = linspace(lower_bound_k, upper_bound_k, n_k);

%% Passo a passo dos slides

% Tentando seguir as instruções dos slides, começarei dando um chute inicial na função política do consumo
% Para esse chute, utilizarei o lado esquerdo da equação de Euler
% apresentada no relatório, já considerando a restrição substituída dentro
% de u'(c).

% Criando uma matriz com zeros de dimensão 500x7

c_old = zeros(n_k,n);

for j = 1:n
    for i = 1:n_k
        c_old(i,j) = grid_z(j).*(grid_k(i).^alfa) + (1 - delta).*grid_k(i) - grid_k(i);
    end
end

% Matriz chute para c preenchida:
c_old;

% Iterando a função consumo
% Definindo uma matriz que vai acumular as atualizações da função valor
c_nova = c_old;

% Definindo valores iniciais para colocar no loop
dif = 1;
iter = 0;

% Definindo a tolerância para fazer a convergência da função valor
tol = 10^-5;

%% ITERAÇÃO

tic
while (dif > tol)
    c_old = c_nova;
    for j = 1:n
        grid_k_egm = zeros(1,n_k);
        for i = 1:n_k
            parcela_1 = zeros(1,n);
            parcela_2 = zeros(1,n);
            for a = 1:n 
                parcela_1(a) = 1 - delta + alfa*grid_z(a)*(grid_k(i)^(alfa - 1));
                parcela_2(a) = c_old(i,a)^(-mu);
            end
            argumento_esperanca = parcela_1.*parcela_2;
            esperanca = P(j,:)*argumento_esperanca';
            u_inversa = (beta*esperanca)^(-1/mu);
            u_inv_soma_k_linha = u_inversa + grid_k(i);
            funcao = @(k)((1-delta)*k + grid_z(j)*(k^alfa)-u_inv_soma_k_linha);
            if i == 1
                k0 = grid_k(i);                   % Chute inicial para achar o 0
            else
                k0 = grid_k(i-1);                % Começamos no elemento anterior do grid
            end
            grid_k_egm(i) = fzero(funcao, k0);
            k_egm(i) = grid_k(i);
        end
        k_exogeno = spline(grid_k_egm, k_egm, grid_k);
        for s = 1:n_k
            c_nova(s,j) = grid_z(j)*(grid_k(s)^alfa) + (1 - delta)*grid_k(s) - k_exogeno(s);
        end
    end
    dif = norm(c_nova - c_old); % tirando a norma
    iter = iter + 1; % contando as iterações
    fprintf('Iterações %4i %6.8f\n ', [iter, dif])
end    
demora = toc


    