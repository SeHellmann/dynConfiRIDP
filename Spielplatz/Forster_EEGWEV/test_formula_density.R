N <- 800


# t0=0
# z=0.5
# d=0
# sz=0
# sv=0
# st0=0
# tau=1
# w=0.5
# muvis=NULL
# sigvis=0
# svis=1
# lambda = 0
# s=1
# simult_conf = FALSE
# precision=1e-5
# z_absolute = FALSE
# stop_on_error=TRUE
# stop_on_zero=FALSE


#         "n", "a", "v", "t0", "z", "d", "sz", "sv", "st0", "tau", "w", "muvis", "sigvis", "svis", "lambda", "s")
df <- rdynaViTE(N,    1,  0.5, 0.1, 0.5,   0,  0.1,     1, 0,       1,    0.5, 0.5,         0,     1,      0.5)
unique(df$response)
parmatrix <- matrix(c(rep(rep(c(1,0.6), each=N/2/2),2), # a
                      rep(c(0.7, 0.3), each=N/2), # v
                      rep(0.1, N), # t0
                      rep(0.5, N), # z
                      rep(0, N), # d
                      rep(0.3, N), # sz
                      rep(1, N), # sv
                      rep(0.3, N), # st0
                      rep(1, N), # tau
                      rep(0.5, N), # w
                      rep(0, N), # sigvis
                      rep(1, N), # svis
                      rep(0.5, N), # lambda
                      rep(c(-4, 4), each=N) # th1 und th2
), nrow=N, byrow = FALSE)
colnames(parmatrix) <- c("a", "v", "t0", "z", "d", "sz", "sv", "st0", "tau", "w", "sigvis", "svis", "lambda", "th1", "th2")
dWEV_parammatrix(rt=df$rt, response=df$response, parammatrix = parmatrix, precision=2)
parmatrix <- matrix(c(rep(rep(c(1,0.6), each=N/2/2),2),
                      rep(c(0.7, 0.3), each=N/2),
                      rep(c(0, 2), each=N)), nrow=N, byrow = FALSE)

colnames(parmatrix) <- c("a", "v", "th1", "th2")
parmatrix
head(parmatrix)
dWEV_parammatrix(rt=df$rt, response=df$response, parammatrix = parmatrix, precision = 2)
ddynaViTE(rt=df$rt, response=df$response, th1 =0, th2 = 2, a=1, v=0.7, t0=0.1, z=0.5, d=0, sz=0.3, sv=1,
     st0=0.3, tau=1, w=0.5, sigvis=0, svis=1, lambda=0.5,  precision = 2)
rt=df$rt
response=df$response
parammatrix = parmatrix
simult_conf = FALSE
precision=1e-5
z_absolute = FALSE
stop_on_error=TRUE
stop_on_zero=FALSE






library(tidyverse)
rt <- rep(seq(0, 4, length.out=200))
response <- c(1,2)
expand.grid(rt=rt, response=response)



a =  1
v=0.5
th1 = 0
th2=2

