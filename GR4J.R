# Modelo GR4J en R
# Autor: Jorge Ivan Ramirez Tamayo jiramirezt@unal.edu.co

# Instrucciones de uso 
# Insertar los datos base del modelo 
df <- read.csv("datos.csv", header=T)
# Insertar los parametros del modelo
x1 <- 1.19346441072412
x2 <- -7.2618083411362
x3 <- 7.88516466971487
x4 <- 0.309427096741274

# Area e la cuenca
area <- 101.63


#Modelo
GR4J <- function(df,x1,x2,x3,x4){
  
  # Parametros modificados
  x1 <- exp(x1)
  x2 <- sinh(x2)
  x3 <- exp(x3)
  x4 <- exp(x4) + 0.5
  
  # Desarrollo de los hidrogramas unitarios sinteticos
  HUS <- data.frame(matrix(nrow = 4, ncol = 20))
  colnames(HUS) <- c(1:20) 
  
  SS1 <- c()
  SS2 <- c()
  HU1 <- c()
  HU2 <- c() 
  
  # Llenado del HUS
  for (i in 1:ncol(HUS)){
    # SS1  
    if(as.numeric(colnames(HUS)[i])<x4){
      val <- as.numeric(colnames(HUS)[i])
      ss1 <- (val/x4)^(2.5)
      
    }else {
      ss1 <- 1 
    }
    SS1[i] <- ss1
    
    
    # SS2
    if(as.numeric(colnames(HUS)[i])<x4){
      val <- as.numeric(colnames(HUS)[i])
      ss2 <- (1/2)*(val/x4)^(2.5)
    }
    else if((as.numeric(colnames(HUS)[i]) > x4) & (as.numeric(colnames(HUS)[i]) < 2*x4)){
      val <- as.numeric(colnames(HUS)[i])
      ss2 <- 1-(1/2)*(2-(val/x4))^(2.5)
    } else{
      ss2 <- 1 
    } 
    SS2[i] <- ss2
  }
  
  HU1[1] <- SS1[1]
  HU2[1] <- SS2[1]
  
  for (i in 2:ncol(HUS)){
    HU1[i] <- SS1[i] - SS1[i-1]
    HU2[i] <- SS2[i] - SS2[i-1]
  }
  
  SS1 <- SS1[1:10]
  HU1 <- HU1[1:10] 
  
  lista_HUS <- list(SS1,SS2,HU1,HU2)
  
  for (i in 1:nrow(HUS)){
    HUS[i,] <- lista_HUS[[i]]
  }
  
  HUS[1,][11:20] <-NA
  HUS[3,][11:20] <-NA
  
  # Valores iniciales 
  S0_x1 <- 0.6
  R0_x3 <- 0.7
  
  vector_sox1 <- c(S0_x1)
  vector_rox1 <- c(R0_x3) 
  lista_HU1 <- list()
  lista_HU2 <- list()
  Pn_acum <- c()
  En_acum <- c()
  Ps_acum <- c()
  Es_acum <- c()
  sx1_f <- c()
  rx3_f <- c()
  Qsim <- c() 
  
  for (k in 1:nrow(df)){
    # Calculo de S/X1
    if (k == 1){
      S_x1 <- S0_x1
    } else {
      S_x1 <- vector_sox1[k]
    }
    
    # Calculo de Pn
    if (df$P[k] >= df$ETP[k]){
      Pn <- df$P[k] - df$ETP[k]
    } else {
      Pn <- 0 
    }
    Pn_acum[k]<- Pn
    
    # Calculo de En
    if (df$P[k] <= df$ETP[k]){
      En <- df$ETP[k] - df$P[k]
    } else {
      En <- 0 
    }
    En_acum[k] <- En
    
    # Calculo de Ps 
    if (Pn >= 0){
      Ps <- (x1*(1-(S_x1)^2)*tanh(Pn/x1))/(1+S_x1*tanh(Pn/x1)) 
    } else{
      Ps <- 0  
    }
    Ps_acum[k] <- Ps
    # Calculo de Es
    if (En >= 0){
      Es <- (S_x1*x1*(2-S_x1)*tanh(En/x1))/(1+(1-S_x1)*tanh(En/x1))  
    } else{
      Es <- 0 
    }
    Es_acum[k] <- Es
    # Calculo de S/X1
    S_x1 <- S_x1 + (Ps-Es)/x1
    
    # Calculo de Perc
    Perc <- S_x1*x1*(1-(1+((4/9)*S_x1)^4)^(-0.25))  
    
    # Calculo de S/X1
    S_x1 <- S_x1 - Perc/x1
    vector_sox1[k+1] <- S_x1 
    sx1_f[k] <- S_x1
    # Calculo del Pr
    Pr <- Perc + (Pn - Ps)
    
    # Calculo del R/x3
    if (k == 1){
      R_x3 <- R0_x3
    } else {
      R_x3 <- vector_rox1[k]
    }
    
    
    # Calculo del F
    F <- x2*(R_x3)^(7/2)
    
    # Calculo de variables asociadas al HU1 y HU2
    v_HU1 <- c()
    v_HU2 <- c()
    
    #HU1
    if (k == 1){
      for (i in 1:length(HU1)){
        # HU1
        vhu1 <- Pr *0.9* HU1[i]
        v_HU1[i] <- vhu1
      }
    } else {
      for (i in 1:length(HU1)){
        # HU1
        if (i < length(HU1)){
          vhu1 <- lista_HU1[[k-1]][i+1]+Pr*0.9* HU1[i]
          v_HU1[i] <- vhu1
        }else{
          # HU1
          vhu1 <- Pr *0.9* HU1[i]
          v_HU1[i] <- vhu1
        }
      }
    } 
    
    #HU2
    if (k == 1){
      for (i in 1:length(HU2)){
        # HU2
        vhu2 <- Pr *0.1* HU2[i]
        v_HU2[i] <- vhu2
      }
    } else {
      for (i in 1:length(HU2)){
        if (i < length(HU2)){
          # HU2
          vhu2 <- lista_HU2[[k-1]][i+1]+Pr *0.1* HU2[i]
          v_HU2[i] <- vhu2
        }else{
          # HU2
          vhu2 <- Pr *0.1* HU2[i]
          v_HU2[i] <- vhu2
        }
      }
    } 
    
    lista_HU1[[k]] <- v_HU1  
    lista_HU2[[k]] <- v_HU2  
    
    # Calculo de R/X3
    R_x3 <- max(0,(R_x3+(v_HU1[1]+F)/x3))
    
    
    # Calculo del QR
    Qr <- R_x3*x3*(1-(1+R_x3^4)^(-0.25)) 
    
    # Calculo de R/X3
    R_x3 <- R_x3 - Qr/x3
    vector_rox1[k+1] <- R_x3
    rx3_f[k] <- R_x3 
    
    # Calcular el QD
    Qd <- max(0,(F+v_HU2[1])) 
    
    # Calcular Q
    Q <- Qr + Qd
    Qsim[k] <- Q 
  }
  return(Qsim)
}

# Se calcula el caudal simulado
Qsim <- GR4J(df,x1,x2,x3,x4)
Qsim <- Qsim/86.4*area 
  
df$Qsim <- Qsim

write.csv(df,"salida.csv",row.names = F)