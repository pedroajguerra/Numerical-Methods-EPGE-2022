%%
clear
clc
%%% Metodos Numéricos (prof. Cezar Santos)
%%% Lista 4 - Pedro Augusto Januzzi Guerra

%% Item (a)

%%% Parametros utilizados: fixos e variaveis (aqueles que vamos mudar nos ultimos itens

% Fixos:
beta = 0.96; % fator de desconto
m = 3; % scaling parameter

% Variaveis:
rho = 0.90; % persistencia do choque de produtividade
sigma = 0.01; % desvio-padrao do choque
gama = 1.0001; % coeficiente de aversao ao risco

%%% Comecarei a questao discretizando o choque via Metodo de Tauchen

% Criando o grid de z
n = 9;
upper_bound = m*sigma*(sqrt(1-rho^2))^-1;
lower_bound = -upper_bound;
grid_z = linspace(lower_bound, upper_bound,n);

% Calculando a matriz de transicao
delta_z = (upper_bound - lower_bound)*((n-1)^-1);
P = zeros(n,n);

% Calculando a primeira e ultima coluna da matriz
P(:,1) = normcdf((grid_z(1)-rho*grid_z+delta_z/2)/sigma);
P(:,n) = 1 - normcdf((grid_z(n)-rho*grid_z-delta_z/2)/sigma);

% Calculando as outras entradas da matriz
for iz=2:(n-1)
    P(:,iz)= normcdf((grid_z(iz)+(delta_z/2)-rho*grid_z)/sigma) - normcdf((grid_z(iz)-(delta_z/2)-rho*grid_z)/sigma);
end

% Matriz de transicao
P;

% Tirando a exponencial do grid dos choques
grid_z = exp(grid_z);
%% Item (b)

%%% Como não foi especificado informacao sobre o tamanho do grid, utilizarei como parametro o trabalho que fizemos para a materia de macro 3
%%% Utilizarei, portanto, 201 pontos no grid, limite inferior igual a -1 e superior igual a 4
%%% Apos rodar o codigo, tentarei rodar outros valores como robustez

% Definindo a taxa de juros de equilibrio deterministico
r = 1/beta - 1;

% Parametros para o grid de ativos:

n_a = 201; % quantidade de pontos no grid
a_min = -1; % limite de endividamento
a_max = 4;

grid_a = linspace(a_min, a_max, n_a); % grid de ativos

%% Informacoes extras para o loop

%%% Vou utilizar monotonicidade: mesmo codigo da lista 2

% Definindo um chute inicial para a função valor com base na restricao e assumindo que a' = a
for ia = 1:n
    for iz = 1:n_a
       c_aux(iz,ia) = grid_z(ia) + r*grid_a(iz); 
    end
end
V_old = ((c_aux.^(1-gama)-1)./(1-gama))./(1-beta);

% Definindo uma matriz que vai acumular as atualizações da função valor
V_nova = V_old;

% Definindo uma matriz de zeros para a policy function, g.
g = zeros(n_a,n);

% Definindo valores iniciais para colocar no loop
dif = 1;
iter = 0;

% Definindo a tolerancia para fazer a convergencia da funcao valor
tol = 10^-5;

%% Iteracao da funcao valor utilizando monotonicidade
tic
while (dif > tol)
    V_old = V_nova; 
    for iz = 1:n 
        for ia = 1:n_a 
            if ia==1
                g_inicial=1; % definindo um ponto inicial para procurar a função política
            end
                        
            % Fazer o for a partir da segunda posicao
            for ja = g_inicial:n_a 
                
                c = grid_z(iz) + (1+r)*grid_a(ia) - grid_a(ja); 
                
                if c >= 0
                    u = (c^(1-gama)-1)/(1-gama) ; 
                else
                    u = -Inf; 
                end
               valor_esperado = 0;
               for s = 1:n 
                    valor_esperado = valor_esperado + P(iz,s)*V_old(ja,s);
               end
               vetor_valor(ja) = u + beta*valor_esperado;  
            end
            V_nova(ia,iz) = max(vetor_valor); 
            indice = find(vetor_valor == V_nova(ia,iz));
            g(ia,iz) = grid_a(indice);
            g_inicial = indice;
        end    
    end 
    dif = max(max(abs((V_nova - V_old))));
    iter = iter + 1;
    fprintf('Iterações %4i %6.8f\n ', [iter, dif])
end
demora = toc

%%
% Resta encontrar a funcao politica para o consumo
z_matriz = repmat(grid_z, n_a, 1); % Criando uma matriz repetindo grid_z para se tornar 201x9
a_matriz = repmat(grid_a, n, 1); % Criando uma matriz repetindo grid_a para se tornar 9x201
a_matriz = a_matriz'; % Transpondo a_matriz para se tornar 201x9

% Encontrando a função política para o consumo
c = z_matriz + (1+r).*a_matriz - g;

%% ERROS DE EQUAÇÃO DE EULER 

% Encontrando Erros da Equacao de Euler

% Criando uma matriz para acumular os EEE

erro_euler = zeros(n_a,n);

tic
for j = 1:n
    for i = 1:n_a
        parcela_1 = zeros(1,n); % primeira parcela dentro do valor esperado da minha equacao de Euler 
        parcela_2 = zeros(1,n); % segunda parcela dentro do valor esperado da minha equacao de Euler
        for a = 1:n
            parcela_1(a) = 1 + r; % essa eh a condicao que muda em relacao as listas anteriores
            idx = find(g(i,j) == grid_a);
            parcela_2(a) = c(idx,a)^(-gama);
        end 
        argumento_esperanca = parcela_1.*parcela_2; % argumento que aparece dentro do operador esperanca;
        esperanca = P(j,:)*argumento_esperanca'; % tirando o valor esperado
        u_inversa = (beta*esperanca)^(-1/gama); % multiplicando beta pelo valor esperado e tirando a inversa da funcao utilidade
        arg = abs(1 - u_inversa/c(i,j)); % tirando o modulo
        erro_euler(i,j) = log10(arg); % erro da equacao de euler
    end    
end
demora = toc
%% DISTRIBUIÇÃO INVARIANTE (ainda estou no item (b))

% Chute inicial para a distribuicao: uniforme
pi_old = ones(n_a,n)/(n_a*n); 

% Definindo uma matriz que vai acumular as atualizacoes da funcao valor
pi_nova = pi_old;

% Parametros para o loop
dif_pi = 1;
iter_pi = 0;

tic
while dif_pi > tol % vou usar a mesma tolerancia utilizada para iteracao da funcao valor
    pi_old = pi_nova;
    for iz = 1:n
        for ia = 1:n_a
            indicadora = (g == grid_a(ia));
            aux = pi_old.*indicadora;
            pi_nova(ia,iz) = sum(aux)*P(:,iz);
        end
    end
    dif_pi = max(max(abs((pi_nova - pi_old))));
    iter_pi = iter_pi + 1;
    fprintf('Iterações %4i %6.8f\n ', [iter_pi, dif_pi])
end
demora = toc


excesso = 0; % calculando excesso de demanda ou oferta sob r = 1/beta - 1
tic
for iz = 1:n
   for ia = 1:n_a
       excesso = excesso + g(ia,iz)*pi_nova(ia,iz);
   end
end
tempo = toc

%% Graficos

plot(grid_a,V_nova);
title('Função Valor');
xlabel('Estoque de Ativos');

mesh(a_matriz,z_matriz,V_nova);
title('Função Valor');
xlabel('Estoque de Ativos');
ylabel('Choques');
zlabel('Função Valor');

plot(c);
title('Função Política Consumo');
xlabel('Estoque de Ativos');

mesh(a_matriz,z_matriz,c);
title('Função Política Consumo');
xlabel('Estoque de Ativos');
ylabel('Choques');
zlabel('Consumo');

plot(g);
title('Função Política Ativos');
xlabel('Estoque de Ativos');

mesh(a_matriz,z_matriz,g);
title('Função Política Ativos');
xlabel('Estoque de Ativos');
ylabel('Choques');
zlabel('Ativos');

plot(erro_euler);
title('Erros da Equação de Euler');
xlabel('Estoque de Ativos');

mesh(a_matriz,z_matriz,erro_euler);
title('Erros da Equação de Euler');
xlabel('Ativos');
ylabel('Choques');
zlabel('EEE');


plot(pi_nova);
title('Distribuição Invariante');
xlabel('Estoque de Ativos');


mesh(a_matriz, z_matriz, pi_nova);
title('Distribuição Invariante');
xlabel('Ativos');
ylabel('Choques');
zlabel('');

%% ITEM (C) 

% As proximas secoes deste codigo se dedicarao a resolver o problema do
% item (c), isto é, nao mais considerando uma taxa de juros dada. Para
% isso, utilizarei a funcao criada, denominada de "inputs", que me gerara
% todos os dados necessarios para calcular a taxa de juros que faz o
% excesso de demanda ou oferta zerar.

%% Estrutura de parametros

% essa estrutura foi utilizada dentro da funcao que crio para utilizar no
% calculo da demanda liquida por titulos. Contudo, incluir essa estrutura
% tornou o codigo mais lento.

parametro.gama = gama;
parametro.rho = rho;
parametro.sigma = sigma;
parametro.m = m ;
parametro.beta = beta;
parametro.grid_a = grid_a;
parametro.grid_z = grid_z;
parametro.n = n;
parametro.n_a = n_a;
parametro.P = P;

%% Encontrar a taxa de juros de equilibrio

%%% Antes, realizei a iteracao da funcao valor com uma taxa de juros
%%% exogena igual a 1/beta - 1. Como temos um modelo em equilibrio geral,
%%% queremos encontrar a taxa de juros que equilibra o mercado de titulos.
%%% Para isso, criei uma funcao chamada inputs_demanda que resolve todo o
%%% problema acima (iteracao funcao valor + distribuicao invariante) e me
%%% retorna a funcao politica dos ativos junto da distribuicao invariante.
%%% Apos isso, realizarei o metodo da bissecao nesta secao a fim de
%%% encontrar a taxa de juros de equilibrio.
%%% Para explicacoes mais detalhadas, checar o relatorio.

%%% Precisamos encontrar a e b tais que f(a) < 0 e f(b) > 0 para que possamos fazer uso do teorema do valor intermediario
%%% Com r = 1/beta - 1, o excesso de oferta eh positivo, logo o excesso de
%%% demanda eh negativo (fui checando)
%%% Com r = 0.04, o excesso de demanda eh positivo (fui checando) 
%%% Entao, utilizarei esses valores como limites iniciais
%%% Uma observacao importante eh que esses testes funcionam para todos os
%%% itens (mudancas de parametros), exceto para o item (f) da lista, para o
%%% qual precisei ajustar e utilizar a = 0.038. 
%%% Como descobri isso? Testando e observando como ficavam as diferencas no
%%% loop utilizado no metodo da bissecao


a = 0.04;
%a = 0.0385; % para o item (f)  
b = 1/beta-1;

dif_bis = 1;
iter = 0;
tol_bis = 10^-5;

tic
while dif_bis > tol_bis 
   [~,~,g_a,~,pi_a] = inputs_demanda(a);
   excesso_a = 0;
   for iz = 1:n
       for ia = 1:n_a
           excesso_a = excesso_a + g_a(ia,iz)*pi_a(ia,iz);
       end
   end
   
   [~,~,g_b,~,pi_b] = inputs_demanda(b);
   excesso_b = 0;
   for izz = 1:n
       for iaa = 1:n_a
           excesso_b = excesso_b + g_b(iaa,izz)*pi_b(iaa,izz);
       end
   end
    
   c = (a+b)/2;
   [V_c,c_c,g_c,eee_c,pi_c] = inputs_demanda(c);
   excesso_c = 0;
   for izzz = 1:n
       for iaaa = 1:n_a
           excesso_c = excesso_c + g_c(iaaa,izzz)*pi_c(iaaa,izzz);
       end
   end

   if excesso_c > 0
       b = c;
       dif_bis = abs(c-a);
   else
       a = c;
       dif_bis = abs(c-b);
   end
   iter = iter + 1;
   %dif_bis = abs(excesso_c);
   fprintf('Iterações %4i %6.8f\n ', [iter, dif_bis]) % para printar o numero da iteracao e o a diferenca 
end
demora = toc

% A taxa de juros que equilibra o mercado de titulos eh:
r_equilibrio = c;

%% Graficos

plot(grid_a,V_c);
title('Função Valor');
xlabel('Estoque de Ativos');

mesh(a_matriz,z_matriz,V_c);
title('Função Valor');
xlabel('Estoque de Ativos');
ylabel('Choques');
zlabel('Função Valor');

plot(grid_a,c_c);
title('Função Política Consumo');
xlabel('Estoque de Ativos');

mesh(a_matriz,z_matriz,c_c);
title('Função Política Consumo');
xlabel('Estoque de Ativos');
ylabel('Choques');
zlabel('Consumo');

plot(grid_a,g_c);
title('Função Política Ativos com zoom');
xlabel('Estoque de Ativos');

mesh(a_matriz,z_matriz,g_c);
title('Função Política Ativos');
xlabel('Estoque de Ativos');
ylabel('Choques');
zlabel('Ativos');

plot(grid_a,eee_c);
title('Erros da Equação de Euler');
xlabel('Estoque de Ativos');

mesh(a_matriz,z_matriz,eee_c);
title('Erros da Equação de Euler');
xlabel('Ativos');
ylabel('Choques');
zlabel('EEE');

mesh(a_matriz, z_matriz, pi_c);
title('Distribuição Invariante');
xlabel('Ativos');
ylabel('Choques');
zlabel('');


