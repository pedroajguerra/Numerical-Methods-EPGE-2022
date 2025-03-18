%%
clear
clc
%%% M�todos Num�ricos (prof. C�zar Santos)
%%% Lista 2 - Pedro Augusto Januzzi Guerra

% % Par�metros utilizados:
% rho = 0.95; % persist�ncia do choque de produtividade
% sigma = 0.007; % desvio-padr�o do choque
% m = 3; % scaling parameter
% beta = 0.987; % fator de desconto
% mu = 2; % coeficiente de avers�o ao risco
% alfa = 1/3; % participa��o do capital na fun��o de produ��o
% delta = 0.012; % deprecia��o do capital

rho = 0.95; % persist�ncia do choque de produtividade
sigma = 0.01; % desvio-padr�o do choque
m = 3; % scaling parameter
beta = 0.95; % fator de desconto
mu = 2; % coeficiente de avers�o ao risco
alfa = 1/3; % participa��o do capital na fun��o de produ��o
delta = 0.05; % deprecia��o do capital

%%% Item (3)

% Come�arei a quest�o discretizando o choque via M�todo de Tauchen
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

% Note que, para cada (k,z), teremos um k' que maximiza V(k,z).
% Logo, nossa fun��o valor ser� uma matriz de n_k por n
% Resolverei, inicialmente, por for�a bruta:

% Definindo um chute inicial para a fun��o valor com base na sugest�o de notas de aula de Ben Moll e sugest�o da Ana Paula (monitora)
for i = 1:n
    for j = 1:n_k
       c(j,i) = grid_z(i)*(grid_k(j)^alfa) + (1-delta)*grid_k(j) - grid_k(j); 
    end
end
V_old = ((c.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old = zeros(n_k,n); % Testei um chute inicial V(.)=0, por�m � um chute ruim e leva mais tempo para rodar.

% Definindo uma matriz que vai acumular as atualiza��es da fun��o valor
V_nova = V_old;

% Definindo uma matriz de zeros para a policy function, g.
g = zeros(n_k,n);

% Definindo valores iniciais para colocar no loop
dif = 1;
iter = 0;

% Definindo a toler�ncia para fazer a converg�ncia da fun��o valor
tol = 10^-5;

%% For�a Bruta
tic
while (dif > tol)
    V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela fun��o valor que o for originou
    for j = 1:n % selecionando choque a choque
        for i = 1:n_k % dado o choque escolhido, percorrendo todos os poss�veis valores de k
            for a = 1:n_k % dado o par (choque, capital), percorrendo novamente o grid_k para selecionar o k'
                c = grid_z(j)*(grid_k(i)^alfa) + (1-delta)*grid_k(i) - grid_k(a); % restri��o que define o consumo, em que k(a) seria o k' na nota��o que estamos habituados
                % calculando a utilidade dado o consumo definido acima
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; % fun��o utilidade CRRA
                else
                    u = -Inf; % condi��o para impedir consumo negativo
                end
               % Criando mais um for para calcular o valor esperado que aparece no final do lado direito da Eq. de Bellman 
               valor_esperado = 0;
               for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado; %vetor que me gera todas as fun��es valores aplicadas nos k'                
               %vetor_valor(a) = u + beta*dot(P(j,:),V_old(a,:)); % alternativa sem loop, mas demorou muito
            end
            % V_nova(i,j) = max(vetor_valor); % pegando a maior fun��o valor aplicada nos k' e tornando essa minha V_nova
            % indice = find(vetor_valor == V_nova(i,j)); % encontrando o �ndice da fun��o valor que gerou minha V_nova
            [V_nova(i,j), indice] = max(vetor_valor); % tentei isso, mas n�o funcionou. Se sobrar tempo, tentar resolver assim tamb�m.
            % V_nova(i,j)=V_n;
            g(i,j) = grid_k(indice); % fun��o pol�tica de capital (k').
        end    
    end 
    % dif = norm(V_nova - V_old); % tirando a norma
    dif = max(max(abs(V_nova-V_old)));
    iter = iter + 1; % contando as itera��es
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif]) % para printar o n�mero da itera��o e o a diferen�a entre V_nova e V_old
end
demora = toc % tempo para rodar o while

%% CONSUMO 
% Resta encontrar a fun��o pol�tica para o consumo
z_matriz = repmat(grid_z, n_k, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz'; %Transpondo k_matriz para se tornar 500x7

% Encontrando a fun��o pol�tica para o consumo, c
c = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - g;

%% MONOTONICIDADE
% Acima, usei for�a bruta. Vamos testar usando monotonicidade para tornar o c�digo mais eficiente.
% A ideia � simples: suponha que cada agente "mora" em um ponto no grid de
% capital. � razo�vel assumir que as decis�es de poupan�a de agentes mais
% ricos ter�o magnitudes maiores que aquelas tomadas por agentes mais
% pobres. Assim, a ideia da monotonicidade �: dado que um agente X escolheu
% um n�vel de k' que mora no grid_k, o agente X+1 escolher� um k'(y) maior
% ou igual ao k'(x).

% As explica��es linha a linha s�o as mesmas que no caso for�a bruta
tic
while (dif > tol)
    V_old = V_nova; 
    for j = 1:n 
        for i = 1:n_k 
            if i==1
                g_inicial=1; % definindo um ponto inicial para procurar a fun��o pol�tica
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
            g_inicial = indice; % dado que k_1 > k_2, ent�o g(k_1) > g(k_2). Assim, atualizo o ponto no grid_k que vou come�ar a procurar o k' maximizador
        end    
    end 
    dif = norm(V_nova - V_old);
    iter = iter + 1;
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif])
end
demora = toc

%%
% Resta encontrar a fun��o pol�tica para o consumo
z_matriz = repmat(grid_z, n_k, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz'; %Transpondo k_matriz para se tornar 500x7

% Encontrando a fun��o pol�tica para o consumo
c = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - g;

%% Gr�ficos

plot(V_nova);
title('Fun��o Valor');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,V_nova);

plot(c);
title('Fun��o Pol�tica Consumo');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,c);

plot(g);
title('Fun��o Pol�tica Capital');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,g);
%% ERROS DE EQUA��O DE EULER

% Encontrando Erros da Equa��o de Euler

% Criando uma matriz para acumular os EEE

erro_euler = zeros(n_k,n);

for j = 1:n
    for i = 1:n_k
        parcela_1 = zeros(1,n); % primeira parcela dentro do valor esperado da minha equa��o de Euler 
        parcela_2 = zeros(1,n); % segunda parcela dentro do valor esperado da minha equa��o de Euler
        for a = 1:n
            parcela_1(a) = 1 - delta + alfa*grid_z(a)*(g(i,j)^(alfa - 1));
            idx = find(g(i,j) == grid_k);
            parcela_2(a) = c(idx,a)^(-mu);
        end 
        argumento_esperanca = parcela_1.*parcela_2; % argumento que aparece dentro do operador esperan�a;
        esperanca = P(j,:)*argumento_esperanca'; % tirando o valor esperado
        u_inversa = (beta*esperanca)^(-1/mu); % multiplicando beta pelo valor esperado e tirando a inversa da fun��o utilidade
        arg = abs(1 - u_inversa/c(i,j)); % tirando o m�dulo
        erro_euler(i,j) = log10(arg); % erro da equa��o de euler
    end    
end

%% Gr�ficos - EEE

plot(erro_euler);
title('Erro de Euler');
xlabel('Estoque de Capital');
mesh(erro_euler);
mesh(k_matriz, z_matriz,erro_euler);
%% ACELERADOR COM FOR�A BRUTA - ITEM (4)

% Vou utilizar o acelerador
tic

while (dif > tol)
    if mod(iter, 10) == 0 || iter <= 100 % testar com limites maximos de itera��o diferentes e com "&&" 
        % a fun��o mod(iter,10) vai pegar o n�mero da �ltima itera��o e
        % dividir por 10. Se o resto for igual a 0 (divis�o exata), calculo
        % o maximo. Caso contr�rio, n�o, mas coloquei p isso acontecer s�
        % depois da cent�sima itera��o.
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
        V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela fun��o valor que o for originou
        for j = 1:n % selecionando choque a choque
            for i = 1:n_k % dado o choque escolhido, percorrendo todos os poss�veis valores de k
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
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif])
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
        V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela fun��o valor que o for originou
        for j = 1:n % selecionando choque a choque
            for i = 1:n_k % dado o choque escolhido, percorrendo todos os poss�veis valores de k
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
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif])
end

demora = toc 

%% CONSUMO
% Resta encontrar a fun��o pol�tica para o consumo
z_matriz = repmat(grid_z, n_k, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz'; %Transpondo k_matriz para se tornar 500x7

% Encontrando a fun��o pol�tica para o consumo, c
c = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - g
%% Gr�ficos

plot(V_nova);
title('Fun��o Valor');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,V_nova);

plot(c);
title('Fun��o Pol�tica Consumo');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,c);

plot(g);
title('Fun��o Pol�tica Capital');
xlabel('Estoque de Capital');

mesh(k_matriz,z_matriz,g);

%% MULTIGRID - ITEM (5)

% A ideia do multigrid � iterar a fun��o valor para um grid mais grosso a
% fim de arrumar um chute inicial melhor. Uma coisa importante de se notar
% � a quantidade de itera��es que levamos para ter aa converg�ncia.
% Utilizando for�a bruta e V_0 = 0, a itera��o da fun��o valor por for�a
% bruta e com 500 pontos no grid leva 1161 itera��es e 407,62 segundos.
% Quando utilizo esse chute inicial para o grid mais grosso (100 pontos),
% leva 17 segundos para realizar a itera��o. Utilizando esse resultado como
% chute inicial para o grid mais grosso (500 pontos), s�o necess�rias 508
% itera��es e 177 segundos, ou seja, no total levo 194 segundos, enquanto
% sem multigrid levo 407. Com 5000 pontos no grid, for�a bruta � muito
% lento. Monotonicidade vai ajudar, mas seria necess�rio acelerar ainda
% mais o c�digo (utilizando concavidade + monotonicidade, por exemplo, mas
% nao sei se vai dar tempo de fazer isso).

clear
clc

% Par�metros utilizados:
rho = 0.95; % persist�ncia do choque de produtividade
sigma = 0.007; % desvio-padr�o do choque
m = 3; % scaling parameter
beta = 0.987; % fator de desconto
mu = 2; % coeficiente de avers�o ao risco
alfa = 1/3; % participa��o do capital na fun��o de produ��o
delta = 0.012; % deprecia��o do capital

% Come�arei a quest�o discretizando o choque via M�todo de Tauchen
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

% Poss�veis tamanhos para o grid do capital
n_k1 = 100;
n_k2 = 500;
n_k3 = 5000;

% Criando o grid do capital
lower_bound_k = 0.75*k_ss;
upper_bound_k = 1.25*k_ss;

grid_k1 = linspace(lower_bound_k, upper_bound_k, n_k1);

% Definindo um chute inicial para a fun��o valor com base na sugest�o de notas de aula de Ben Moll e sugest�o da Ana Paula (monitora)
for i = 1:n
    for j = 1:n_k1
       c(j,i) = grid_z(i)*(grid_k1(j)^alfa) + (1-delta)*grid_k1(j) - grid_k1(j); 
    end
end
V_old = ((c.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old = zeros(n_k1,n);

% Definindo uma matriz que vai acumular as atualiza��es da fun��o valor
V_nova = V_old;

% Definindo uma matriz de zeros para a policy function, g.
g = zeros(n_k1,n);

% Definindo valores iniciais para colocar no loop
dif = 1;
iter = 0;

% Definindo a toler�ncia para fazer a converg�ncia da fun��o valor
tol = 10^-5;

% Farei, incialmente, por for�a bruta
tic
while (dif > tol)
    V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela fun��o valor que o for originou
    for j = 1:n % selecionando choque a choque
        for i = 1:n_k1 % dado o choque escolhido, percorrendo todos os poss�veis valores de k
            for a = 1:n_k1 % dado o par (choque, capital), percorrendo novamente o grid_k para selecionar o k'
                c = grid_z(j)*(grid_k1(i)^alfa) + (1-delta)*grid_k1(i) - grid_k1(a); % restri��o que define o consumo, em que k(a) seria o k' na nota��o que estamos habituados
                % calculando a utilidade dado o consumo definido acima
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; % fun��o utilidade CRRA
                else
                    u = -Inf; % condi��o para impedir consumo negativo
                end
               % Criando mais um for para calcular o valor esperado que aparece no final do lado direito da Eq. de Bellman 
               valor_esperado = 0;
               for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado; %vetor que me gera todas as fun��es valores aplicadas nos k'                
               %vetor_valor(a) = u + beta*dot(P(j,:),V_old(a,:)); % alternativa sem loop, mas demorou muito
            end
            V_nova(i,j) = max(vetor_valor); % pegando a maior fun��o valor aplicada nos k' e tornando essa minha V_nova
            indice = find(vetor_valor == V_nova(i,j)); % encontrando o �ndice da fun��o valor que gerou minha V_nova
            %[V_nova, index] = max(vetor_valor) % tentei isso, mas n�o funcionou. Se sobrar tempo, tentar resolver assim tamb�m.
            g(i,j) = grid_k1(indice); % fun��o pol�tica de capital (k').
        end    
    end 
    dif = norm(V_nova - V_old); % tirando a norma
    iter = iter + 1; % contando as itera��es
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif]) % para printar o n�mero da itera��o e o a diferen�a entre V_nova e V_old
end
demora = toc % tempo para rodar o while

V_nova_k1 = V_nova; % apenas para armazenar os valores dessa itera��o
g_k1 = g; % apenas para armazenar os valores dessa itera��o

%% Pr�ximo passo: grid com 500 pontos (ainda com for�a bruta)

% Criando o segundo grid de capital, agora com 500 pontos
grid_k2 = linspace(lower_bound_k, upper_bound_k, n_k2);

% Informa��es pr�vias que precisamos retomar
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
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif]) 
end
demora = toc 

V_nova_k2 = V_nova;
g_k3 = g;

%% Pr�ximo passo: grid com 5000 pontos (ainda com for�a bruta)

% Criando o segundo grid de capital, agora com 5000 pontos
grid_k3 = linspace(lower_bound_k, upper_bound_k, n_k3);

% Informa��es pr�vias que precisamos retomar
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
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif]) 
end
demora = toc 

%% MULTIGRID COM MONOTONICIDADE

clear
clc

% Par�metros utilizados:
rho = 0.95; % persist�ncia do choque de produtividade
sigma = 0.007; % desvio-padr�o do choque
m = 3; % scaling parameter
beta = 0.987; % fator de desconto
mu = 2; % coeficiente de avers�o ao risco
alfa = 1/3; % participa��o do capital na fun��o de produ��o
delta = 0.012; % deprecia��o do capital

% Come�arei a quest�o discretizando o choque via M�todo de Tauchen
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

% Poss�veis tamanhos para o grid do capital
n_k1 = 100;
n_k2 = 500;
n_k3 = 5000;

% Criando o grid do capital
lower_bound_k = 0.75*k_ss;
upper_bound_k = 1.25*k_ss;

grid_k1 = linspace(lower_bound_k, upper_bound_k, n_k1);

% Definindo um chute inicial para a fun��o valor com base na sugest�o de notas de aula de Ben Moll e sugest�o da Ana Paula (monitora)
for i = 1:n
    for j = 1:n_k1
       c(j,i) = grid_z(i)*(grid_k1(j)^alfa) + (1-delta)*grid_k1(j) - grid_k1(j); 
    end
end
V_old = ((c.^(1-mu)-1)./(1-mu))./(1-beta);
%V_old = zeros(n_k1,n);

% Definindo uma matriz que vai acumular as atualiza��es da fun��o valor
V_nova = V_old;

% Definindo uma matriz de zeros para a policy function, g.
g = zeros(n_k1,n);

% Definindo valores iniciais para colocar no loop
dif = 1;
iter = 0;

% Definindo a toler�ncia para fazer a converg�ncia da fun��o valor
tol = 10^-5;

% Farei, incialmente, por for�a bruta
tic
while (dif > tol)
    V_old = V_nova; % depois que os for's abaixo rodarem, atualizo a minha V_old pela fun��o valor que o for originou
    for j = 1:n % selecionando choque a choque
        for i = 1:n_k1 % dado o choque escolhido, percorrendo todos os poss�veis valores de k
            if i == 1 
                g_inicial = 1;
            end
            for a = g_inicial:n_k1 % dado o par (choque, capital), percorrendo novamente o grid_k para selecionar o k'
                c = grid_z(j)*(grid_k1(i)^alfa) + (1-delta)*grid_k1(i) - grid_k1(a); % restri��o que define o consumo, em que k(a) seria o k' na nota��o que estamos habituados
                % calculando a utilidade dado o consumo definido acima
                if c > 0
                    u = (c^(1-mu)-1)/(1-mu) ; % fun��o utilidade CRRA
                else
                    u = -Inf; % condi��o para impedir consumo negativo
                end
               % Criando mais um for para calcular o valor esperado que aparece no final do lado direito da Eq. de Bellman 
               valor_esperado = 0;
               for s = 1:n % repare que estou percorrendo o grid de choques
                    valor_esperado = valor_esperado + P(j,s)*V_old(a,s);
               end
               vetor_valor(a) = u + beta*valor_esperado; %vetor que me gera todas as fun��es valores aplicadas nos k'                
               %vetor_valor(a) = u + beta*dot(P(j,:),V_old(a,:)); % alternativa sem loop, mas demorou muito
            end
            V_nova(i,j) = max(vetor_valor); % pegando a maior fun��o valor aplicada nos k' e tornando essa minha V_nova
            indice = find(vetor_valor == V_nova(i,j)); % encontrando o �ndice da fun��o valor que gerou minha V_nova
            %[V_nova, index] = max(vetor_valor) % tentei isso, mas n�o funcionou. Se sobrar tempo, tentar resolver assim tamb�m.
            g(i,j) = grid_k1(indice); % fun��o pol�tica de capital (k').
            g_inicial = indice;
        end    
    end 
    dif = norm(V_nova - V_old); % tirando a norma
    iter = iter + 1; % contando as itera��es
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif]) % para printar o n�mero da itera��o e o a diferen�a entre V_nova e V_old
end
demora = toc % tempo para rodar o while

V_nova_k1 = V_nova;
g_k1 = g;

%% Pr�ximo passo: grid com 500 pontos (agora com monotonicidade)

% Criando o segundo grid de capital, agora com 500 pontos
grid_k2 = linspace(lower_bound_k, upper_bound_k, n_k2);

% Informa��es pr�vias que precisamos retomar
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
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif]) 
end
demora = toc 

V_nova_k2 = V_nova;
g_k2 = g;

%% Pr�ximo passo: grid com 5000 pontos (agora com monotonicidade)
n_k3 = 5000;
% Criando o segundo grid de capital, agora com 5000 pontos
grid_k3 = linspace(lower_bound_k, upper_bound_k, n_k3);

% Informa��es pr�vias que precisamos retomar
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
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif]) 
end
demora = toc 

V_nova_k3 = V_nova;
g_k3 = g;
%% Consumo
% Resta encontrar a fun��o pol�tica para o consumo
z_matriz = repmat(grid_z, n_k3, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k3, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz'; %Transpondo k_matriz para se tornar 500x7

% Encontrando a fun��o pol�tica para o consumo
c = z_matriz.*(k_matriz.^alfa) + (1-delta).*k_matriz - g_k3;

%% ERROS DE EQUA��O DE EULER

% Encontrando Euler Equation Errors

% Criando uma matriz para acumular os EEE

erro_euler_mult = zeros(n_k3,n);

for j = 1:n
    for i = 1:n_k3
        parcela_1 = zeros(1,n); % primeira parcela dentro do valor esperado da minha equa��o de Euler 
        parcela_2 = zeros(1,n); % segunda parcela dentro do valor esperado da minha equa��o de Euler
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

%% Gr�ficos
z_matriz = repmat(grid_z, n_k3, 1); %Criando uma matriz repetindo grid_z para se tornar 500x7
k_matriz = repmat(grid_k3, n, 1); %Criando uma matriz repetindo grid_z para se tornar 7x500
k_matriz = k_matriz';

plot(V_nova);
title('Fun��o Valor - Multigrid');
xlabel('Estoque de Capital');

mesh(V_nova);

mesh(erro_euler_mult);

plot(c);
title('Fun��o Pol�tica Consumo - Multigrid');
xlabel('Estoque de Capital');

mesh(c);

plot(g);
title('Fun��o Pol�tica Capital - Multigrid');
xlabel('Estoque de Capital');

mesh(g);

plot(erro_euler_mult);
title('Erro de Euler - Multigrid');
xlabel('Estoque de Capital');

%surf(grid_k3, grid_z, V_nova,'FaceColor','interp');
%title('Fun��o Valor - Multigrid');
%colormap(hot);

%surf(grid_k3, grid_z, c,'FaceColor','interp');
%title('Fun��o Valor - Multigrid');
%colormap(hot);

%surf(grid_k3, grid_z, g,'FaceColor','interp');
%title('Fun��o Valor - Multigrid');
%colormap(hot);

%surf(grid_k3, grid_z, erro_euler_mult,'FaceColor','interp');
%title('Fun��o Valor - Multigrid');
%colormap(hot);
%% GRID END�GENO - ITEM (6)

% Repetindo as informa��es iniciais como nos itens anteriores
clear
clc

% Par�metros utilizados:
rho = 0.95; % persist�ncia do choque de produtividade
sigma = 0.007; % desvio-padr�o do choque
m = 3; % scaling parameter
beta = 0.987; % fator de desconto
mu = 2; % coeficiente de avers�o ao risco
alfa = 1/3; % participa��o do capital na fun��o de produ��o
delta = 0.012; % deprecia��o do capital

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

% Criando o grid do capital (n�o � o end�geno, obviamente)
n_k = 500;
lower_bound_k = 0.75*k_ss;
upper_bound_k = 1.25*k_ss;

grid_k = linspace(lower_bound_k, upper_bound_k, n_k);

%% Passo a passo dos slides

% Tentando seguir as instru��es dos slides, come�arei dando um chute inicial na fun��o pol�tica do consumo
% Para esse chute, utilizarei o lado esquerdo da equa��o de Euler
% apresentada no relat�rio, j� considerando a restri��o substitu�da dentro
% de u'(c).

% Criando uma matriz com zeros de dimens�o 500x7

c_old = zeros(n_k,n);

for j = 1:n
    for i = 1:n_k
        c_old(i,j) = grid_z(j).*(grid_k(i).^alfa) + (1 - delta).*grid_k(i) - grid_k(i);
    end
end

% Matriz chute para c preenchida:
c_old;

% Iterando a fun��o consumo
% Definindo uma matriz que vai acumular as atualiza��es da fun��o valor
c_nova = c_old;

% Definindo valores iniciais para colocar no loop
dif = 1;
iter = 0;

% Definindo a toler�ncia para fazer a converg�ncia da fun��o valor
tol = 10^-5;

%% ITERA��O

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
                k0 = grid_k(i-1);                % Come�amos no elemento anterior do grid
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
    iter = iter + 1; % contando as itera��es
    fprintf('Itera��es %4i %6.8f\n ', [iter, dif])
end    
demora = toc


    