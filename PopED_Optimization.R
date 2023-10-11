#################### paper R code #################### 
set.seed(1013)
install.packages("PopED")
library(PopED)



### example 1 ###
## set-up for example 1
{
  # 1) model function - using the above solution, we define the model function
  ff1 <- function(model_switch,xt,parameters,poped.db){
    with(as.list(parameters),{
      y=xt
      N = floor(xt/TAU)+1
      y=(DOSE*Favail/V)*(KA/(KA - CL/V)) * 
        (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) - 
           exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
      return(list( y=y,poped.db=poped.db))
    })
  }
  
  # 2) parameter definition function
  sfg1 <- function(x,a,bpop,b,bocc){
    parameters=c( V=bpop[1]*exp(b[1]),
                  KA=bpop[2]*exp(b[2]),
                  CL=bpop[3]*exp(b[3]),
                  Favail=bpop[4],
                  DOSE=a[1],
                  TAU=a[2])
    return( parameters ) 
  }
  
  # 3) residual function
  feps1 <- function(model_switch,xt,parameters,epsi,poped.db){
    returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
    y <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    
    y = y*(1+epsi[,1])+epsi[,2]
    
    return(list( y= y,poped.db =poped.db )) 
  }
  
  
  # 4) design space
  db11_dispersed <- create.poped.database(ff_fun="ff1",
                                          fg_fun="sfg1",
                                          fError_fun="feps1",
                                          bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                          notfixed_bpop=c(1,1,1,0),
                                          d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                          sigma=c(0.04,5e-6),
                                          notfixed_sigma=c(1,0),
                                          m=1,
                                          groupsize=40,
                                          xt=c(50, 100, 150, 200, 250),
                                          minxt=c(0,0,0,0,0),
                                          maxxt=c(250,250,250,250,250),
                                          bUseGrouped_xt=0,
                                          a=list(c(DOSE=20,TAU=24)),
                                          maxa=c(DOSE=100,TAU=24),
                                          mina=c(DOSE=0,TAU=12))
  db11_skewed <- create.poped.database(ff_fun="ff1",
                                       fg_fun="sfg1",
                                       fError_fun="feps1",
                                       bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                       notfixed_bpop=c(1,1,1,0),
                                       d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                       sigma=c(0.04,5e-6),
                                       notfixed_sigma=c(1,0),
                                       m=1,
                                       groupsize=40,
                                       xt=c(6, 12, 24, 216, 240), #this is 6,12h and 1,9,10days after
                                       minxt=c(0,0,0,0,0),
                                       maxxt=c(250,250,250,250,250),
                                       bUseGrouped_xt=0,
                                       a=list(c(DOSE=20,TAU=24)),
                                       maxa=c(DOSE=100,TAU=24),
                                       mina=c(DOSE=0,TAU=12))
  
  db11_windowadj <- create.poped.database(ff_fun="ff1",
                                          fg_fun="sfg1",
                                          fError_fun="feps1",
                                          bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                          notfixed_bpop=c(1,1,1,0),
                                          d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                          sigma=c(0.04,5e-6),
                                          notfixed_sigma=c(1,0),
                                          m=1,
                                          groupsize=40,
                                          xt=c(50, 100, 150, 200, 250),
                                          minxt=c(0,50,100,150,200),
                                          maxxt=c(50,100,150,200,250),
                                          bUseGrouped_xt=0,
                                          a=list(c(DOSE=20,TAU=24)),
                                          maxa=c(DOSE=100,TAU=24),
                                          mina=c(DOSE=0,TAU=12))
  
  db11_optimized  <- create.poped.database(ff_fun="ff1",
                                           fg_fun="sfg1",
                                           fError_fun="feps1",
                                           bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                           notfixed_bpop=c(1,1,1,0),
                                           d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                           sigma=c(0.04,5e-6),
                                           notfixed_sigma=c(1,0),
                                           m=1,
                                           groupsize=40,
                                           xt=c(0.568,  14.15 , 14.15,  14.15,240),
                                           minxt=c(0,0,0,0,0),
                                           maxxt=c(250,250,250,250,250),
                                           bUseGrouped_xt=0,
                                           a=list(c(DOSE=20,TAU=24)),
                                           maxa=c(DOSE=100,TAU=24),
                                           mina=c(DOSE=0,TAU=12))
  
  
  db11_optimized_adj <- create.poped.database(ff_fun="ff1",
                                              fg_fun="sfg1",
                                              fError_fun="feps1",
                                              bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                              notfixed_bpop=c(1,1,1,0),
                                              d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                              sigma=c(0.04,5e-6),
                                              notfixed_sigma=c(1,0),
                                              m=1,
                                              groupsize=40,
                                              xt=c(15.31,  51.02,  146.9,    168,  242.9),
                                              minxt=c(0,50,100,150,200),
                                              maxxt=c(50,100,150,200,250),
                                              bUseGrouped_xt=0,
                                              a=list(c(DOSE=20,TAU=24)),
                                              maxa=c(DOSE=100,TAU=72),
                                              mina=c(DOSE=0,TAU=0))
  db12_freq1 <- create.poped.database(ff_fun="ff1",
                                      fg_fun="sfg1",
                                      fError_fun="feps1",
                                      bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                      notfixed_bpop=c(1,1,1,0),
                                      d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                      sigma=c(0.02,5e-6),
                                      notfixed_sigma=c(1,0),
                                      m=1,
                                      groupsize=40,
                                      xt=c(15.31,  51.02, 144,194.8, 218.8),
                                      minxt=c(0,50,100,150,200),
                                      maxxt=c(50,100,150,200, 250),
                                      bUseGrouped_xt=1,
                                      a=list(c(DOSE=10,TAU=12)),
                                      maxa=c(DOSE=200,TAU=72),
                                      mina=c(DOSE=0,TAU=0))
  db12_freq2 <- create.poped.database(ff_fun="ff1",
                                      fg_fun="sfg1",
                                      fError_fun="feps1",
                                      bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                      notfixed_bpop=c(1,1,1,0),
                                      d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                      sigma=c(0.02,5e-6),
                                      notfixed_sigma=c(1,0),
                                      m=1,
                                      groupsize=40,
                                      xt=c(15.31,  51.02, 144,194.8, 218.8),
                                      minxt=c(0,50,100,150,200),
                                      maxxt=c(50,100,150,200, 250),
                                      bUseGrouped_xt=1,
                                      a=list(c(DOSE=5,TAU=6)),
                                      maxa=c(DOSE=200,TAU=72),
                                      mina=c(DOSE=0,TAU=0))
  db12_less1 <- create.poped.database(ff_fun="ff1",
                                      fg_fun="sfg1",
                                      fError_fun="feps1",
                                      bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                      notfixed_bpop=c(1,1,1,0),
                                      d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                      sigma=c(0.02,5e-6),
                                      notfixed_sigma=c(1,0),
                                      m=1,
                                      groupsize=40,
                                      xt=c(15.31,  51.02, 144,194.8, 218.8),
                                      minxt=c(0,50,100,150,200),
                                      maxxt=c(50,100,150,200, 250),
                                      bUseGrouped_xt=1,
                                      a=list(c(DOSE=30,TAU=36)),
                                      maxa=c(DOSE=200,TAU=72),
                                      mina=c(DOSE=0,TAU=0))
  
  db12_less2 <- create.poped.database(ff_fun="ff1",
                                      fg_fun="sfg1",
                                      fError_fun="feps1",
                                      bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9), 
                                      notfixed_bpop=c(1,1,1,0),
                                      d=c(V=0.09,KA=0.09,CL=0.25^2), 
                                      sigma=c(0.02,5e-6),
                                      notfixed_sigma=c(1,0),
                                      m=1,
                                      groupsize=40,
                                      xt=c(15.31,  51.02, 144,194.8, 218.8),
                                      minxt=c(0,50,100,150,200),
                                      maxxt=c(50,100,150,200, 250),
                                      bUseGrouped_xt=1,
                                      a=list(c(DOSE=40,TAU=48)),
                                      maxa=c(DOSE=200,TAU=72),
                                      mina=c(DOSE=0,TAU=0))
  
  
  
}
## simulation
{
  
  
  # model prediction
  plot_model_prediction(db11_dispersed,  model_num_points = 300, PI=TRUE,
                        y_lab = "Concentration")
  # model optimization
  set.seed(1013)
  output1 <- poped_optim(db11_dispersed, opt_xt =TRUE, opt_a = FALSE,
                         parallel=TRUE, dSeed = 1013)
  set.seed(1013) #setting the seed again here is required
  output1_adj <- poped_optim(db11_windowadj, opt_xt =TRUE, opt_a = FALSE,
                             parallel=FALSE, dSeed = 1013)
  
  summary(output1)
  summary(output1_adj)
  
  
  # model evaluation
  #db11_dispersed db11_skewed db12_less1 db12_less2 db12_1 db12_freq1 db12_freq2 
  #table 1
  (e12_D1dispersed <- evaluate_design(db11_dispersed))
  (e11_D2skewed <- evaluate_design(db11_skewed))
  (e13_D3optimized <- evaluate_design(db11_optimized))
  (e14_D4optimized_adj <- evaluate_design(db11_optimized_adj))
  efficiency(e11_D2skewed$ofv, e12_D1dispersed$ofv, db11_dispersed)
  efficiency(e13_D3optimized$ofv, e12_D1dispersed$ofv, db11_dispersed)
  efficiency(e14_D4optimized_adj$ofv, e12_D1dispersed$ofv, db11_dispersed)
  
  
  
  #table2
  (less2_D41 <- evaluate_design(db12_less2))
  (less1_D42 <- evaluate_design(db12_less1))
  (freq1_D43 <- evaluate_design(db12_freq1))
  (freq2_D44 <- evaluate_design(db12_freq2))
  efficiency(less2_D41$ofv, e14_D4optimized_adj$ofv, db11_optimized_adj)
  efficiency(less1_D42$ofv, e14_D4optimized_adj$ofv, db11_optimized_adj)
  efficiency(freq1_D43$ofv, e14_D4optimized_adj$ofv, db11_optimized_adj)
  efficiency(freq2_D44$ofv, e14_D4optimized_adj$ofv, db11_optimized_adj)
  
  
  
}




### example 2 ###
## set-up for example 2
{
  #############excerpt from Shotwell(2016)#############START
  pk_basic_solution <-
    function(k_10, k_12, k_21, v_1, k_R, c_0=c(0,0)) {
      K <- matrix(c(-(k_10+k_12), k_21, k_12, -k_21),
                  2,2, byrow=TRUE)
      T_K <- -(k_10 + k_12 + k_21)
      D_K <- (k_10+k_12)*k_21 - k_12*k_21
      lambda <- T_K/2 + c(-1,1)*(T_K^2/4 - D_K)^(1/2)
      N_21 <- (k_21+lambda[2])/k_21
      N_22 <- (k_21+lambda[1])/k_21
      N <- matrix(c(-1,N_21,-1,N_22),2,2)
      gamma <- c(k_R/k_10/v_1, k_12*k_R/k_10/k_21/v_1)
      r <- c(( N_22*(c_0[1]-(k_R/v_1/k_10)) +
                 (c_0[2]-(k_R/v_1/k_10)*k_12/k_21)),
             (-N_21*(c_0[1]-(k_R/v_1/k_10)) -
                (c_0[2]-(k_R/v_1/k_10)*k_12/k_21)))
      r <- r * k_21/(diff(lambda))
      c_1 <- function(t)
        gamma[1] + r[1]*N[1,1]*exp(lambda[1]*t) +
        r[2]*N[1,2]*exp(lambda[2]*t)
      c_2 <- function(t)
        gamma[2] + r[1]*N[2,1]*exp(lambda[1]*t) +
        r[2]*N[2,2]*exp(lambda[2]*t)
      return(list(c_1=c_1, c_2=c_2))
    }
  #############excerpt from Shotwell(2016)#############END
  
  # 1) model function - using the above solution, we define the model function
  
  #for the specific model presented in Shotwell(2016)
  shotwell_function <- function(xt, k_R, k_10, k_12, k_21, v_1,
                                c_0=c(0,0), type= "central") {
    ###excerpt from Shotwell(2016)START
    # find solution for first hour: 3000 mg/hr IV infusion
    sol0 <- pk_basic_solution(k_10, k_12, k_21, v_1,
                              k_R, c_0=c(0,0))
    # find solution for next two hours: regular elimination
    sol1 <- pk_basic_solution(k_10, k_12, k_21, v_1,
                              k_R=0, c_0=c(sol0$c_1(1),sol0$c_2(1)))
    # find solution for next two hours: dialytic clearance 15L/h
    sol2 <- pk_basic_solution(k_10+15/v_1, k_12, k_21, v_1,
                              k_R=0, c_0=c(sol1$c_1(2),sol1$c_2(2)))
    # find solution for next three hours: regular elimination
    sol3 <- pk_basic_solution(k_10, k_12, k_21, v_1,
                              k_R=0, c_0=c(sol2$c_1(2),sol2$c_2(2)))
    ###excerpt from Shotwell(2016)END
    central_y1 = {sol0$c_1(xt)*(xt >= 0 & xt < 1) +
        sol1$c_1(xt-1)*(xt >= 1 & xt < 3) + sol2$c_1(xt-3)*(xt >= 3 & xt < 5) +
        sol3$c_1(xt-5)*(xt >= 5)}
    peripheral_y2 = {sol0$c_2(xt)*(xt >= 0 & xt < 1) +
        sol1$c_2(xt-1)*(xt >= 1 & xt < 3) + sol2$c_2(xt-3)*(xt >= 3 & xt < 5) +
        sol3$c_2(xt-5)*(xt >= 5)}
    if(type == "central") {y <- central_y1}
    else if(type == "peripheral") {y <- peripheral_y2}
    return(y)
  } 
  
  
  
  myff_central <- function(model_switch,xt,parameters,poped.db){
    with(as.list(parameters),{
      y=xt
      y=shotwell_function(xt=xt, k_R=k_R, k_10=k_10, k_12=k_12, k_21=k_21,
                          v_1=v_1, c_0=c_0, type="central")
      return(list( y=y,poped.db=poped.db))
    })
  }
  
  myff_peripheral <- function(model_switch,xt,parameters,poped.db){
    with(as.list(parameters),{
      y=xt
      y=shotwell_function(xt=xt, k_R=k_R, k_10=k_10, k_12=k_12, k_21=k_21,
                          v_1=v_1, c_0=c_0, type="peripheral")
      return(list( y=y,poped.db=poped.db))
    })
  }
  
  # 2) parameter definition function
  mysfg <- function(x,a,bpop,b,bocc){
    parameters=c( k_10=bpop[1]*exp(b[1]),
                  k_12=bpop[2]*exp(b[2]),
                  k_21=bpop[3]*exp(b[3]),
                  v_1=bpop[4]*exp(b[4]),
                  k_R=a)
    return(parameters) 
  }
  
  # 3) residual function
  feps <- function(model_switch,xt,parameters,epsi,poped.db){
    returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
    y <- returnArgs[[1]]
    poped.db <- returnArgs[[2]]
    
    y = y*(1+epsi[,1])+epsi[,2]
    
    return(list( y= y,poped.db =poped.db )) 
  }
  
  
  # 4) design space
  db3_central <- create.poped.database(ff_fun="myff_central",
                                       fg_fun="mysfg",
                                       fError_fun="feps",
                                       bpop=c(k_10=0.35, k_12=3.50, k_21=1.50, v_1=10), 
                                       notfixed_bpop=c(1,1,1,1),
                                       d=c(k_10=log((0.25^2) + 1), k_12=log((0.25^2) + 1),
                                           k_21=log((0.25^2) + 1), v_1=log((0.25^2) + 1)),
                                       sigma=c(0.04,5e-6),
                                       notfixed_sigma=c(1,0),
                                       groupsize=100,
                                       xt=c(2,4,6,8),
                                       minxt=c(0,0,0,0),
                                       maxxt=c(8,8,8,8),
                                       bUseGrouped_xt=1,
                                       a=list(c(a=3000)),
                                       maxa=c(DOSE=4500),
                                       mina=c(DOSE=1500))
  
  db3_central_alt <- create.poped.database(ff_fun="myff_central",
                                           fg_fun="mysfg",
                                           fError_fun="feps",
                                           bpop=c(k_10=0.35, k_12=3.50, k_21=1.50, v_1=10), 
                                           notfixed_bpop=c(1,1,1,1),
                                           d=c(k_10=log((0.25^2) + 1), k_12=log((0.25^2) + 1),
                                               k_21=log((0.25^2) + 1), v_1=log((0.25^2) + 1)),
                                           sigma=c(0.04,5e-6),
                                           notfixed_sigma=c(1,0),
                                           groupsize=100,
                                           xt=c(0.5,1,6,8),
                                           minxt=c(0,0,0,0),
                                           maxxt=c(8,8,8,8),
                                           bUseGrouped_xt=1,
                                           a=list(c(a=3000)),
                                           maxa=c(DOSE=4500),
                                           mina=c(DOSE=1500))
  
  db3_peripheral <- create.poped.database(ff_fun="myff_peripheral",
                                          fg_fun="mysfg",
                                          fError_fun="feps",
                                          bpop=c(k_10=0.35, k_12=3.50, k_21=1.50, v_1=10), 
                                          notfixed_bpop=c(1,1,1,1),
                                          d=c(k_10=log((0.25^2) + 1), k_12=log((0.25^2) + 1),
                                              k_21=log((0.25^2) + 1), v_1=log((0.25^2) + 1)),
                                          sigma=c(0.04,5e-6),
                                          notfixed_sigma=c(1,0),
                                          groupsize=100,
                                          xt=c(2,4,6,8),
                                          minxt=c(0,0,0,0),
                                          maxxt=c(8,8,8,8),
                                          bUseGrouped_xt=1,
                                          a=list(c(a=3000)),
                                          maxa=c(DOSE=4500),
                                          mina=c(DOSE=1500))
  
  db3_peripheral_alt <- create.poped.database(ff_fun="myff_peripheral",
                                              fg_fun="mysfg",
                                              fError_fun="feps",
                                              bpop=c(k_10=0.35, k_12=3.50, k_21=1.50, v_1=10), 
                                              notfixed_bpop=c(1,1,1,1),
                                              d=c(k_10=log((0.25^2) + 1), k_12=log((0.25^2) + 1),
                                                  k_21=log((0.25^2) + 1), v_1=log((0.25^2) + 1)),
                                              sigma=c(0.04,5e-6),
                                              notfixed_sigma=c(1,0),
                                              groupsize=100,
                                              xt=c(0.5,1,6,8),
                                              minxt=c(0,0,0,0),
                                              maxxt=c(8,8,8,8),
                                              bUseGrouped_xt=1,
                                              a=list(c(a=3000)),
                                              maxa=c(DOSE=4500),
                                              mina=c(DOSE=1500))
  
}


## simulation
{
  # model prediction
  plot_model_prediction(db3_central, IPRED=TRUE, IPRED.lines = FALSE, DV=T, PI=TRUE,
                        model_num_points = 300, y_lab = "Concentration")
  plot_model_prediction(db3_peripheral, IPRED=FALSE, IPRED.lines = FALSE, DV=T,  PI=TRUE,
                        model_num_points = 300, y_lab = "Concentration")
  # model evaluation
  (e3_central <- evaluate_design(db3_central))
  (e3_central_alt  <- evaluate_design(db3_central_alt))
  (e3_peripheral <- evaluate_design(db3_peripheral))
  (e3_peripheral_alt  <- evaluate_design(db3_peripheral_alt))
  
  efficiency(e3_central_alt$ofv, e3_central$ofv, db3_central)
  efficiency(e3_peripheral_alt$ofv, e3_peripheral$ofv, db3_peripheral)
  
  
  # model optimization
  #central 
  set.seed(1013)
  output3_central <- poped_optim(db3_central, opt_xt =TRUE, opt_a = FALSE,
                                 parallel=TRUE, dSeed=1013)
  
  #database function for optimized design
  db3_central_opt <- create.poped.database(ff_fun="myff_central",
                                           fg_fun="mysfg",
                                           fError_fun="feps",
                                           bpop=c(k_10=0.35, k_12=3.50, k_21=1.50, v_1=10), 
                                           notfixed_bpop=c(1,1,1,1),
                                           d=c(k_10=log((0.25^2) + 1), k_12=log((0.25^2) + 1),
                                               k_21=log((0.25^2) + 1), v_1=log((0.25^2) + 1)),
                                           sigma=c(0.04,5e-6),
                                           notfixed_sigma=c(1,0),
                                           groupsize=100,
                                           xt=c(0.0008256,  1.125,  3.463,   8),
                                           minxt=c(0,0,0,0),
                                           maxxt=c(8,8,8,8),
                                           bUseGrouped_xt=1,
                                           a=list(c(a=3000)),
                                           maxa=c(DOSE=4500),
                                           mina=c(DOSE=1500))
  e3_central_opt <- evaluate_design(db3_central_opt)
  summary(output3_central); efficiency(e3_central_opt$ofv, e3_central$ofv, db3_central)
  
  # central discrete
  db3_central_dis <- create.poped.database(db3_central,discrete_xt = list(0:8))
  set.seed(1013)
  output3_central_dis <- poped_optim(db3_central_dis, opt_xt =TRUE, opt_a = FALSE, 
                                     parallel=TRUE, dSeed=1013)
  #database function for optimized design
  db3_central_dis_opt <- create.poped.database(ff_fun="myff_central",
                                               fg_fun="mysfg",
                                               fError_fun="feps",
                                               bpop=c(k_10=0.35, k_12=3.50, k_21=1.50, v_1=10), 
                                               notfixed_bpop=c(1,1,1,1),
                                               d=c(k_10=log((0.25^2) + 1), k_12=log((0.25^2) + 1),
                                                   k_21=log((0.25^2) + 1), v_1=log((0.25^2) + 1)),
                                               sigma=c(0.04,5e-6),
                                               notfixed_sigma=c(1,0),
                                               groupsize=100,
                                               xt=c(1,2,4,8),
                                               minxt=c(0,0,0,0),
                                               maxxt=c(8,8,8,8),
                                               bUseGrouped_xt=1,
                                               a=list(c(a=3000)),
                                               maxa=c(DOSE=4500),
                                               mina=c(DOSE=1500))
  e3_central_dis_opt <- evaluate_design(db3_central_dis_opt)
  summary(output3_central_dis); efficiency(e3_central_dis_opt$ofv, e3_central$ofv, db3_central)
  
  
  
  #peripheral
  set.seed(1013)
  output3_peripheral <- poped_optim(db3_peripheral, opt_xt =TRUE, opt_a = FALSE,
                                    parallel=TRUE, dSeed=1013)
  #database function for optimized design
  db3_peripheral_opt <- create.poped.database(ff_fun="myff_central",
                                              fg_fun="mysfg",
                                              fError_fun="feps",
                                              bpop=c(k_10=0.35, k_12=3.50, k_21=1.50, v_1=10), 
                                              notfixed_bpop=c(1,1,1,1),
                                              d=c(k_10=log((0.25^2) + 1), k_12=log((0.25^2) + 1),
                                                  k_21=log((0.25^2) + 1), v_1=log((0.25^2) + 1)),
                                              sigma=c(0.04,5e-6),
                                              notfixed_sigma=c(1,0),
                                              groupsize=100,
                                              xt=c(0.0144,  3.041,  5.133,      8),
                                              minxt=c(0,0,0,0),
                                              maxxt=c(8,8,8,8),
                                              bUseGrouped_xt=1,
                                              a=list(c(a=3000)),
                                              maxa=c(DOSE=4500),
                                              mina=c(DOSE=1500))
  e3_peripheral_opt <- evaluate_design(db3_peripheral_opt)
  summary(output3_peripheral); efficiency(e3_peripheral_opt$ofv, e3_peripheral$ofv, db3_peripheral)
  
  # peripheral discrete
  db3_peripheral_dis <- create.poped.database(db3_peripheral,discrete_xt = list(0:8))
  set.seed(1013)
  output3_peripheral_dis <- poped_optim(db3_peripheral_dis, opt_xt =TRUE, opt_a = FALSE,
                                        parallel=TRUE, dSeed=1013)
  #database function for optimized design
  db3_peripheral_dis_opt <- create.poped.database(ff_fun="myff_central",
                                                  fg_fun="mysfg",
                                                  fError_fun="feps",
                                                  bpop=c(k_10=0.35, k_12=3.50, k_21=1.50, v_1=10), 
                                                  notfixed_bpop=c(1,1,1,1),
                                                  d=c(k_10=log((0.25^2) + 1), k_12=log((0.25^2) + 1),
                                                      k_21=log((0.25^2) + 1), v_1=log((0.25^2) + 1)),
                                                  sigma=c(0.04,5e-6),
                                                  notfixed_sigma=c(1,0),
                                                  groupsize=100,
                                                  xt=c(1,3,5,8),
                                                  minxt=c(0,0,0,0),
                                                  maxxt=c(8,8,8,8),
                                                  bUseGrouped_xt=1,
                                                  a=list(c(a=3000)),
                                                  maxa=c(DOSE=4500),
                                                  mina=c(DOSE=1500))
  e3_peripheral_dis_opt <- evaluate_design(db3_peripheral_dis_opt)
  summary(output3_peripheral_dis); efficiency(e3_peripheral_dis_opt$ofv, e3_peripheral$ofv, db3_peripheral)
  
  
  
  
  #Table 5
  round(evaluate_design(db3_central)$rse)
  round(evaluate_design(db3_central_alt)$rse)
  round(evaluate_design(db3_central_opt)$rse)
  round(evaluate_design(db3_central_dis_opt)$rse)
  
  round(evaluate_design(db3_peripheral)$rse)
  round(evaluate_design(db3_peripheral_alt)$rse)
  round(evaluate_design(db3_peripheral_opt)$rse)
  round(evaluate_design(db3_peripheral_dis_opt)$rse)
  
}
