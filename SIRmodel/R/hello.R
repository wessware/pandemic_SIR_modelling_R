#defining the three functions for SIR
dSdt <- function(t,S, I) {
  return(-beta * S * I /N)
}
dIdt <- function(t, S, I) {
  return(beta * S * I /N - gamma* I)
}
Drdt <- function(t, I) {
  return(gamma * I)
}
#define the Runge Kutta fourth order method for the SIR

RK4SIR <- function(n, beta, gamma, S0, I0, R0 = 0, dt = 1) {
  N <<- S0 + I0 + R0 #the fixed population

  S <- c(S0, rep(0, n))
  I <- c(I0, rep(0, n))
  R <- c(R0, rep(0, n))
  for (i in 1:n) {
    Si <- S[i]
    Ii <- I[i]

    #Ri <- R[i]

    S.k1 <- dSdt(i, Si, Ii)
    I.k1 <- DIdt(i, Si, Ii)

    #R.k1 <- dRdt(i, Ii)

    S.k2 <- dSdt(i + dt / 2, Si + dt / 2 * S.k1, Ii + dt / 2 * I.k1)
    I.k2 <- dIdt(i + dt / 2, Si + dt / 2 * S.k1, Ii + dt / 2 * I.k1)
    #     R.k2 <- dRdt(i + dt / 2, Ii + dt / 2 * I.k1)

    S.k3 <- dSdt(i + dt / 2, Si + dt / 2 * S.k2, Ii + dt / 2 * I.k2)
    I.k3 <- dIdt(i + dt / 2, Si + dt / 2 * S.k2, Ii + dt / 2 * I.k2)
    #     R.k3 <- dRdt(i + dt / 2, Ii + dt / 2 * I.k2)

    S.k4 <- dSdt(i + dt, Si + dt * S.k3, Ii + dt * I.k3)
    I.k4 <- dIdt(i + dt, Si + dt * S.k3, Ii + dt * I.k3)
    #     R.k4 <- dRdt(i + dt, Ii + dt * I.k3)

    S[i + 1] <- Si + dt / 6 * (S.k1 + 2 * S.k2 + 2 * S.k3 + S.k4)
    I[i + 1] <- Ii + dt / 6 * (I.k1 + 2 * I.k2 + 2 * I.k3 + I.k4)
    #     R[i + 1] <- Ri + dt / 6 * (R.k1 + 2 * R.k2 + 2 * R.k3 + R.k4)
  }
  R <- N - S - I
  return(data.frame(n = 0:n, S = S, I = I, R = R))

S0 <- 7900000
I0 <- 10
beta <- 1 / 2
gamma <- 1 / 3
n <- 200

r <- RK4SIR(n, beta, gamma, S0, I0)
}


