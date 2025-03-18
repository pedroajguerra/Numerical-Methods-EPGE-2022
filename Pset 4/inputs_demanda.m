function [V, c, pol_g, eee, pi] = inputs_demanda(r)

%%% Parametros utilizados: fixos e variaveis (aqueles que vamos mudar nos ultimos itens

% Fixos:
beta = 0.96; % fator de desconto
m = 3; % scaling parameter

% Variaveis:
rho = 0.9; % persistencia do choque de produtividade
sigma = 0.01; % desvio-padrao do choque
gama = 1.0001; % coeficiente de aversao ao risco
%gama = 5;
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

n_a = 201; % quantidade de pontos no grid
a_min = -1; % limite de endividamente
a_max = 4;

grid_a = linspace(a_min, a_max, n_a); % grid de ativos

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

% Definindo a tolerância para fazer a convergência da função valor
tol = 10^-5;

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
end

% Resta encontrar a funcao politica para o consumo
z_matriz = repmat(grid_z, n_a, 1); % Criando uma matriz repetindo grid_z para se tornar 201x9
a_matriz = repmat(grid_a, n, 1); % Criando uma matriz repetindo grid_a para se tornar 9x201
a_matriz = a_matriz'; % Transpondo a_matriz para se tornar 500x7

% Encontrando a função política para o consumo
c = z_matriz + (1+r).*a_matriz - g;

% Encontrando Erros da Equacao de Euler

% Criando uma matriz para acumular os EEE

erro_euler = zeros(n_a,n);

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

% Encontrando a distribuicao invariante

% Chute inicial para a distribuicao
pi_old = ones(n_a,n)/(n_a*n);

% Definindo uma matriz que vai acumular as atualizacoes da funcao valor
pi_nova = pi_old;

% Parametros para o loop
dif_pi = 1;

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
end


V = V_nova;
c = c;
pol_g = g;
eee = erro_euler;
pi = pi_nova;

end

