function [V, c, pol_g, eee, pi] = inputs(parametro,r)
%%% Vou utilizar monotonicidade: mesmo codigo da lista 2

% Definindo um chute inicial para a função valor com base na restricao e assumindo que a' = a
for ia = 1:parametro.n
    for iz = 1:parametro.n_a
       c_aux(iz,ia) = parametro.grid_z(ia) + r*parametro.grid_a(iz); 
    end
end
V_old = ((c_aux.^(1-parametro.gama)-1)./(1-parametro.gama))./(1-parametro.beta);

% Definindo uma matriz que vai acumular as atualizações da função valor
V_nova = V_old;

% Definindo uma matriz de zeros para a policy function, g.
g = zeros(parametro.n_a,parametro.n);

% Definindo valores iniciais para colocar no loop
dif = 1;

% Definindo a tolerancia para fazer a convergencia da funcao valor
tol = 10^-5;

while (dif > tol)
    V_old = V_nova; 
    for iz = 1:parametro.n 
        for ia = 1:parametro.n_a 
            if ia==1
                g_inicial=1; % definindo um ponto inicial para procurar a função política
            end
                        
            % Fazer o for a partir da segunda posicao
            for ja = g_inicial:parametro.n_a 
                
                c = parametro.grid_z(iz) + (1+r)*parametro.grid_a(ia) - parametro.grid_a(ja); 
                
                if c >= 0
                    u = (c^(1-parametro.gama)-1)/(1-parametro.gama) ; 
                else
                    u = -Inf; 
                end
               valor_esperado = 0;
               for s = 1:parametro.n 
                    valor_esperado = valor_esperado + parametro.P(iz,s)*V_old(ja,s);
               end
               vetor_valor(ja) = u + parametro.beta*valor_esperado;  
            end
            V_nova(ia,iz) = max(vetor_valor); 
            indice = find(vetor_valor == V_nova(ia,iz));
            g(ia,iz) = parametro.grid_a(indice);
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

% Chute inicial para a distribuicao invariante
pi_old = ones(parametro.n_a,parametro.n)/(parametro.n_a*parametro.n);

% Definindo uma matriz que vai acumular as atualizacoes da funcao valor
pi_nova = pi_old;

% Parametros para o loop
dif_pi = 1;


while dif_pi > tol % vou usar a mesma tolerancia utilizada para iteracao da funcao valor
    pi_old = pi_nova;
    for iz = 1:parametro.n
        for ia = 1:parametro.n_a
            indicadora = (g == parametro.grid_a(ia));
            aux = pi_old.*indicadora;
            pi_nova(ia,iz) = sum(aux)*parametro.P(:,iz);
        end
    end
    dif_pi = max(max(abs((pi_nova - pi_old))));
end

V = V_nova;
c = c_matriz;
pol_g = g;
eee = erro_euler;
pi = pi_nova;

end
