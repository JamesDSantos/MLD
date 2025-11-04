#=====Models

require(magic)
require(mvtnorm)
require(Matrix)

#-----F e G para os modelos polinomiais (de crescimento)---------
modPoly <- function(order,m0=rep(0,order), C0 = diag(1,order),delta=.99){
  if (order < 1) stop("A ordem deve ser pelo menos 1.")
  Ft <- matrix(c(1, rep(0, order - 1)), ncol = 1)
  
  Gt <- diag(1, order) # Cria uma matriz identidade
  if (order > 1) {
    # Preenche a primeira super-diagonal com 1s
    indices <- seq.int(from = order + 1, by = order + 1, length.out = order - 1)
    Gt[indices] <- 1
  }
  
  q = order 
  mod <- list( Ft = Ft, Gt = Gt, m = m0, C = C0, q = q, delta = delta)
  class(mod) = c("dglm")   
  return(mod)
}

#-----Seasonal model - form free ---------
modSaz = function(p,m0=rep(0,p-1),C0=diag(1,p-1),W=diag(1,p-1) ,delta = .8,V=1){
  W0 = W
  Ft = matrix(c(1,rep(0,p-2)) , ncol = 1)
  Gt = rbind( rep( -1,p-1 ) , cbind( diag(1,p-2),0 ) )
  q = length(m0)
  mod = list(Ft =Ft,Gt=Gt,m=m0,C=C0,W=W0,q = q, delta=delta,arma=F,V=V)
  class(mod) = c("dglm")
  return(mod)
}

#-----F e G para os modelos sazonais (forma trigonometrica)---------
modTrig = function(p, harm = 1, delta = .8){
  require(Matrix)
  
  max.h <- floor(p/2) # Número máximo de harmônicos
  
  #--- Suas excelentes mensagens de erro (mantidas) ---
  if(any(harm > max.h) || length(unique(harm)) != length(harm)){
    stop("Verifique os harmônicos solicitados. Eles não podem se repetir ou exceder p/2.")
  }
  #----------------------------------------------------
  
  # A matriz de rotação está perfeita
  J2 <- function(h) {
    omega <- h * 2 * pi / p
    matrix(c(cos(omega), sin(omega), -sin(omega), cos(omega)), 2, byrow = TRUE)
  }
  
  Ft_list <- list()
  Gt_list <- list()
  
  # Loop único que trata todos os casos
  for (h in harm) {
    # Verifica se é o caso especial de Nyquist
    if (h == p/2 && p %% 2 == 0) {
      Ft_list[[length(Ft_list) + 1]] <- 1
      Gt_list[[length(Gt_list) + 1]] <- -1 # Correção principal: deve ser -1
    } else {
      # Caso normal para todos os outros harmônicos
      Ft_list[[length(Ft_list) + 1]] <- c(1, 0)
      Gt_list[[length(Gt_list) + 1]] <- J2(h)
    }
  }
  
  # Monta as matrizes finais a partir das listas
  Ft <- matrix(unlist(Ft_list), ncol = 1)
  Gt <- as.matrix(bdiag(Gt_list))
  
  # O resto do seu código está ótimo
  q <- dim(Gt)[1]
  m0 <- rep(0, q)
  C0 <- diag(1, q)
  W0 <- C0
  mod <- list(Ft = Ft, Gt = Gt, m = m0, C = C0, q = q, delta = delta)
  class(mod) <- c("dglm")
  
  return(mod)
}

#------Superposition theorem----------------------------
"+.dglm" <- function(mod1,mod2){
  
  delta = c(mod1$delta,mod2$delta)
  q = c(mod1$q,mod2$q)
  
  Ft = rbind(mod1$Ft,mod2$Ft)
  Gt <- as.matrix( bdiag(mod1$Gt, mod2$Gt))
  C <- as.matrix( bdiag(mod1$C, mod2$C))
  delta = c(mod1$delta,mod2$delta)
  m = c( mod1$m, mod2$m)
  mod = list( Ft=Ft , Gt=Gt , m=m , C=C,q = q,delta=delta)
  class(mod) = c("dglm")
  return(mod)
}


#=======Discount
W_discount = function(Pt,mod){
  
  delta = mod$delta
  q = mod$q
  Wt    = array(0, c( sum(q), sum(q)))
  k <- 0
  for(i in 1:length(q)){
    Wt[k+1:q[i], k+1:q[i]] <- Pt[k+1:q[i], k+1:q[i]]*( 1/delta[i] - 1)
    k <- cumsum(q)[i]
  }
  Wt <- .5*( Wt + t(Wt))
  return(Wt)
}
#=======Filter
dglmFilter_normal = function(yt,mod, m0 = NULL, C0 = NULL, nt=1,st=.01, Vt= 'unknown', Wt = 'discount', outliers = NULL){
  n <- length(yt)
  Ft <- mod$Ft
  Gt <- mod$Gt
  q <- dim(Gt)[1]
  
  
  # priori at time t
  at <- array(NA_real_, c(n,q))
  Rt <- array(NA_real_, c(n,q,q))
  
  # forecast for yt
  ft <- array(NA_real_, n)
  Qt <- array(NA_real_, n)
  
  # posterior at time t
  mt <- array( 0, c(n+1,q))
  Ct <- array( NA_real_, c(n+1,q,q))
  
  # initial information
  if(is.null(m0)){
    mt[1,][1] <- yt[1] 
  } else {
    mt[1,] <- m0
  }
  if(is.null(C0)){
    Ct[1,,] <- diag(1,q) 
  } else {
    Ct[1,,] <- C0
  }
  
  if(is.null(outliers)){
    outliers <- rep(1,n)
  } else{
    one_vector <- rep(1,n)
    one_vector[outliers] <- 0
    outliers <- one_vector
  }
  
  if( Vt == 'unknown'){ Vt <- 1} 
  # filtro de kalman sem desconto
  if( Wt != 'discount'){
    for(t in 1:n){
      # atualiza a priori
      at[t,]   <- Gt %*% mt[t,]
      Rt[t,,]  <- Gt %*% Ct[t,,] %*% t(Gt) + Wt
      
      # previsão para o tempo t
      ft[t] <- t(Ft) %*% at[t,]
      Qt[t] <- t(Ft) %*% Rt[t,,] %*% Ft + Vt
      
      # atualiza a posteriori
      At <- Rt[t,,] %*% Ft %*% solve(Qt[t])
      mt[ t + 1, ] <- at[t,] + outliers[t] *(At %*% ( yt[t]- ft[t]))
      Ct[ t+1,,]   <- Rt[t,,] - outliers[t]*(At %*% Qt[t] %*% t(At))
      
      
      nt[t+1] <- nt[t] + 1
      st[t+1] <- (nt[t]/nt[t+1])*st[t] + ( yt[t] - ft[t])^2 / (nt[t+1]*Qt[t]) 
    }
  }
  # filtro de kalman com desconto
  if( Wt == 'discount'){
    for(t in 1:n){
      # atualiza a priori
      at[t,]   <- Gt %*% mt[t,]
      Pt <- Gt %*% Ct[t,,] %*% t(Gt)
      Rt[t,,]  <- Pt + W_discount(Pt, mod)
      
      # previsão para o tempo t
      ft[t] <- t(Ft) %*% at[t,]
      Qt[t] <- t(Ft) %*% Rt[t,,] %*% Ft + Vt
      
      # atualiza a posteriori
      At <- Rt[t,,] %*% Ft %*% solve(Qt[t])
      mt[ t + 1, ] <- at[t,] + outliers[t]*(At %*% ( yt[t]- ft[t]))
      Ct[ t+1,,]   <- Rt[t,,] - outliers[t]*(At %*% Qt[t] %*% t(At))
      
      
      nt[t+1] <- nt[t] + 1
      st[t+1] <- (nt[t]/nt[t+1])*st[t] + ( yt[t] - ft[t])^2 / (nt[t+1]*Qt[t]) 
    }
  }
  if( Vt == 'unknown'){
    Rt <- sapply( 1:n, function(t) st[t]*Rt[t,,])
    Qt <- sapply( 1:n, function(t) st[t]*Qt[t])
    Ct <- sapply( 1:(n+1), function(t) st[t]*Ct[t,,])
  }
  
  resp <- list( yt = yt, mt = mt, Ct = Ct, at = at,
                Rt = Rt, ft = ft, Qt = Qt, Ft = Ft, Gt = Gt)
  class(resp) <- 'filter_dglm'
  return(resp)
  
}
