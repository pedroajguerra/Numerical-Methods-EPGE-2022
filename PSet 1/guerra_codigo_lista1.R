### Lista 1 - M�todos Num�ricos
### Pedro Augusto Januzzi Guerra

install.packages("stargazer")
install.packages("pracma")
install.packages("dplyr")

library(stargazer)
library(pracma)
library(dplyr)

rho <- 0.95
sigma <- 0.007
mu = 0 # m�dia do erro
m = 3

### Item 1

# Criando o grid de z

n <- 9
upper_bound <- m*sigma*(sqrt(1-rho^2)^-1)
lower_bound <- -upper_bound

grid_z <- linspace(lower_bound, upper_bound, n)

# Calculando a matriz de transi��o

delta_z <- (upper_bound - lower_bound)*((n-1)^-1)
P <- matrix(0, n, n)

# Calculando a primeira e �ltima coluna da matriz

P[,1] <- pnorm((grid_z[1]-rho*grid_z+delta_z/2)/sigma)
P[,n] <- 1 - pnorm((grid_z[n]-rho*grid_z-delta_z/2)/sigma)

# Calculando as outras entradas da matriz

for (i in 2:(n-1)){
  P[,i] <- pnorm((grid_z[i] + delta_z/2 - rho*grid_z)/sigma) - pnorm((grid_z[i] - delta_z/2 - rho*grid_z)/sigma)
  
}

# Matriz de transi��o

P

### Item 2

# Definindo vari�veis 

sigma_z <- sqrt(sigma^2*((1-rho^2)^-1))
upper_bound_rou <- sigma_z*sqrt(n-1)
lower_bound_rou <- -upper_bound_rou

grid_z_rou <- linspace(lower_bound_rou, upper_bound_rou, n)     

p  = (1+rho)/2  

P_rou <- matrix(c(p ,1-p, 1-p, p), 2, 2) # matriz P2

# Computando a matriz P_rou recursivamente

for (i in 3:n){
  # Criando matrizes auxiliares com todas as entradas iguais a zero
  aux_1 <- matrix(0,i,i)
  aux_2 <- matrix(0,i,i)
  aux_3 <- matrix(0,i,i)
  aux_4 <- matrix(0,i,i)
  
  # Substituindo algumas entradas das matrizes auxiliares pela matriz de transi��o P referente ao estado anterior
  aux_1[1:i-1,1:i-1] <- P_rou
  aux_2[1:i-1,2:i] <- P_rou
  aux_3[2:i,1:i-1] <- P_rou
  aux_4[2:i,2:i] <- P_rou
  
  P_rou <- p*aux_1 + (1-p)*aux_2 + (1-p)*aux_3 + p*aux_4
}

# Para normalizar as linhas: 

for (i in 1:n){
  P_rou[i,] = P_rou[i,]/sum(P_rou[i,])
}

# Matriz de transi��o

P_rou

### Item 3

set.seed(1) # botando um seed para padronizar o processo gerador de dados aleat�rios
t <- 10000 # n�mero de per�odos
choque <- rnorm(t,mu,sigma)

# Simulando o processo AR(1)

z_t <- rep(0,t)

for (i in 2:t){
  z_t[i] <- rho*z_t[i-1] + choque[i]
}

# Processo AR(1)
z_t

# Gr�fico AR(1)
plot(z_t, type = 'l', col = 'blue', 
     main = 'Processo AR(1) simulado para 10000 per�odos', ylim = c(-0.08,0.08),
     xlab = "Per�odos")  

# Simulando o processo discretizado via m�todo de Tauchen

# Pegando o ponto no grid_z em que z = 0 (i.e., a mediana)

mediana <- median(grid_z)

# Encontrando o estado inicial

estado <- match(mediana, grid_z)

# Criando um vetor para armazenar os estados em que estou (come�ando da mediana = 0)

z_discretizado <- c(mediana, rep(0,t-1))

for (i in 2:t){
  cdf_estado = P[estado,1] # considerando que estou no estado inicial igual a 5, qual a probabilidade de ir para cada estado
  for (j in 1:n){
    if (pnorm(choque[i]/sigma) < cdf_estado) {
      estado <- j
      z_discretizado[i] <- grid_z[estado]
      break
    } else {
      cdf_estado <- P[estado,j+1] + cdf_estado
    }
  }
}

# Processo discretizado (Tauchen)
z_discretizado

# Criando o gr�fico do processo discretizado (Tauchen)
plot(z_discretizado, type = 'l',col='red', 
     main = 'Processo discretizado via Tauchen',
     ylim = c(-0.08,0.08),xlab = "Per�odos")

# Plotando a simula��o do AR(1) e do processo discretizado via Tauchen juntos
plot(z_t, type = 'l', col = 'blue', ylab = "z", xlab = "Per�odos", 
     ylim = c(-0.08,0.08)) 
par(new=TRUE)
plot(z_discretizado, type = 'l',col='red', 
     main = 'AR(1) com processo discretizado via Tauchen', ylim = c(-0.08,0.08),
     ylab = "", xlab = "Per�odos")

# Simulando o processo discretizado via m�todo de Rouwenhorst

mediana_rou <- median(grid_z_rou)

# Encontrando o estado inicial

estado_rou <- match(mediana_rou, grid_z_rou)

# Criando um vetor para armazenar os estados em que estou (come�ando da mediana = 0)

z_discretizado_rou <- c(mediana_rou, rep(0,t-1))

for (i in 2:t){
  cdf_estado_rou = P_rou[estado_rou,1] # considerando que estou no estado inicial igual a 5, qual a probabilidade de ir para cada estado
  for (j in 1:n){
    if (pnorm(choque[i]/sigma) < cdf_estado_rou) {
      estado_rou <- j
      z_discretizado_rou[i] <- grid_z_rou[estado_rou]
      break
    } else {
      cdf_estado_rou <- P_rou[estado_rou,j+1] + cdf_estado_rou
    }
  }
}

# Processo discretizado (Rouwenhorst)
z_discretizado_rou

# Criando o gr�fico do processo discretizado (Rouwenhorst)
plot(z_discretizado_rou, type = 'l',col='red', 
     main = 'Processo discretizado via Rouwenhorst',
     ylim = c(-0.08,0.08),xlab = "Per�odos")

# Plotando a simula��o do AR(1) e do processo discretizado via Rouwenhorst juntos
plot(z_t, type = 'l', col='blue', ylab = "z", xlab = "Per�odos", ylim = c(-0.08,0.08))
par(new=TRUE)
plot(z_discretizado_rou, type = 'l',col='red', 
     main = 'AR(1) com processo discretizado via Rouwenhorst', ylab = "", xlab = "",
     ylim = c(-0.08,0.08))
### Item 4

# Rodando a regress�o do z_discretizado no seu pr�prio lag (m�todo de Tauchen)

reg_tauchen <- lm(z_discretizado ~ -1 + lag(z_discretizado, 1))
summary(reg_tauchen)

# Rodando a regress�o do z_discretizado no seu pr�prio lag (m�todo de Rouwenhorst)

reg_rou <- lm(z_discretizado_rou ~ -1 + lag(z_discretizado_rou, 1))
summary(reg_rou)

# Calculando o intervalo de confian�a de 95% para o coeficiente da regress�o usando Tauchen

ci_tauchen <- confint(reg_tauchen)
ci_tauchen

# Calculando o intervalo de confian�a de 95% para o coeficiente da regress�o usando Rouwenhorst

ci_rou <- confint(reg_rou)
ci_rou

'Se considerarmos rho = 0.95, vemos que o intervalo de confian�a para o coeficiente
da regress�o (tanto usando Tauchen quanto Rouwenhorst) cont�m 0.95'

'J� no caso de rho = 0.99, podemos ver que Tauchen n�o funciona muito bem, 
o que � evidenciado pelo fato de que 0.99 n�o est� contido no intervalo de 
confian�a de 95%. No caso da regress�o usando Rouwenhorst, rho = 0.99 est� 
contido no intervalo de confian�a.'

### Item 5

'Basta repetir todo o c�digo acima, apenas alterando o rho para 0.99. 
Nesse caso, conseguimos visualizar graficamente que utilizamos poucos pontos no grid.
Uma alternativa � aumentar o grid. Realizo esse teste colocando 101 pontos no grid 
e processo discretizado se aproxima ainda mais do processo cont�nuo.'

