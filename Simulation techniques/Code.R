params <- list(m = 4, x0 = 13, eta =4, lambda0 = 30, alpha = 0.5, v = 2, h0 = 300)

densite_non_normalisee <- function(x, params)
{
   1/ (params$m + abs(x - params$x0)^params$eta)
}

u <- 0.5 + atan(-params$x0/sqrt(params$m))/pi

coef_inverse <- (params$m/(params$m+1))*(sqrt(params$m))/pi

rintensite <- function(n,params)
{
  L <- c()
  while(length(L)<n){
    
    x <- runif(n-length(L),u,1)
    
    x <- qcauchy(x,location = params$x0, scale=sqrt(params$m))
    
    y <- runif(length(x), min=0, max=1)
    
    D <- x[coef_inverse*densite_non_normalisee(x,params)>=y*dcauchy(x,params$x0,sqrt(params$m))]
    
    L <- append(L,D)
    
  }
  
  return(L[1:n])  
}

lambda <- function(t) {
  params$lambda0 * (1 + params$alpha * sin(4 * pi * t))
}

t_pluies <- function(T,intensite)
{
  X <- rexp(floor(params$lambda0), intensite)
  Y <- cumsum(X)
  if (Y[floor(params$lambda0)] > T)
  {
    Tn <- Y[Y < T]
  }
  else
  {
    Tn <- Y
    TN <- Tn[length(Tn)]
    T0 <- TN
    X <- rexp(1, intensite)
    TN <- T0 + X
    while (TN < T)
    {
      Tn <- c(Tn, TN)
      T0 <- TN
      X <- rexp(1, intensite)
      TN <- T0 + X
      
    }
  }
  U <- runif(length(Tn), 0, 1)
  X <- Tn[intensite * U <= lambda(Tn)]
  return(X)
}

simul_h_pluies <- function(Tp, params) {
  I <- rintensite(length(Tp), params)
  H <- I[1]
  for (i in 1:(length(Tp) - 1))
  {
    if (H>params$h0){return(1)}
    H <- H * exp(params$v * (Tp[i] - Tp[i+1])) + I[i + 1]
    
  }
  if (H>params$h0){return(1)}
  return(0)
}

petoile<-function(n.simul,params)
{
  TH <- 0
  for(i in 1:n.simul)
  {
    TH<-TH+simul_h_pluies(t_pluies(1,params$lambda0 * (1 + params$alpha)),params)
  }
  result <- list(p = TH / n.simul, demi.largeur = 1.96 * sqrt((TH / n.simul) * (1 - TH / n.simul) / n.simul))
  return(result)
}




#petoile.petite, on multiplie le nombre moyen de pluies par k.
k=0.011

t_pluies_frequentes <- function(T,intensite)
{
  X <- rexp(floor(params$lambda0), intensite)
  Y <- cumsum(X)
  if (Y[floor(params$lambda0)] > T)
  {
    Tn <- Y[Y < T]
  }
  else
  {
    Tn <- Y
    TN <- Tn[length(Tn)]
    T0 <- TN
    X <- rexp(1, intensite)
    TN <- T0 + X
    while (TN < T)
    {
      Tn <- c(Tn, TN)
      T0 <- TN
      X <- rexp(1, intensite)
      TN <- T0 + X
      
    }
  }
  U <- runif(length(Tn), 0, 1)
  X <- Tn[intensite * U <= (1+k)*lambda(Tn)]
  return(X)
}

simul_h_grande_pluies <- function(Tp, params) {
  I <- rintensite(length(Tp), params)
  H <- I[1]
  for (i in 1:(length(Tp) - 1))
  {
    if (H>params$h0){return(exp(k*params$lambda0)*(1/(1+k))^length(Tp))}
    
    H <- H * exp(params$v * (Tp[i] - Tp[i+1])) + I[i+1]
    
  }
  if (H>params$h0){return(exp(k*params$lambda0)*(1/(1+k))^length(Tp))}
  return(0)
}

petoile.petite<-function(n.simul,params)
{
  TH <- 0
  for(i in 1:n.simul)
  {
    TH<-TH+simul_h_grande_pluies(t_pluies_frequentes(1,(1+k)*params$lambda0 * (1 + params$alpha)),params)
  }
  result <- list(p = TH / n.simul, demi.largeur = 1.96 * sqrt((TH / n.simul) * (1 - TH / n.simul) / n.simul))
  return(result)
}

#Utilisation C++
petoile.c <- function(n.simul,params)
{
  library(Rcpp)
  cppFunction("
int hauteurRcpp(NumericVector a, NumericVector b, double paramshauteur, double paramsv) {
  
  double H = b[0];
  
  for (int i = 0; i < a.size() - 1; ++i) {
    if (H > paramshauteur) {
      return 1;
    }
    H = H * exp(paramsv * (a[i] - a[i + 1])) + b[i + 1];
  }
  
  if (H > paramshauteur) {
    return 1;
  }
  
  return 0;
}")
  cppFunction("
	NumericVector t_pluies_rcpp(double paramslambda0, double intensite, double T, double paramsalpha)
      {
		RNGScope scope;
            NumericVector n(2);
            NumericVector X(floor(paramslambda0)), U(floor(paramslambda0 * 10));
            n = rpois(2,intensite * T);
            int lenL = 0;
            NumericVector Tn(2 * n[1]);
            Tn=runif(n[1], 0, T);
            U=runif(n[1], 0, T);
            for(int i = 0; i < n[1]; i++)
            {
            	if(intensite * U[i] <= paramslambda0 * (1 + paramsalpha * sin(4 * M_PI * Tn[i])))
            	{
                  	lenL = lenL + 1;
            	}
            }
            NumericVector L(lenL);
            int lenA = 0;
            for(int i = 0; i < n[1]; i++)
            {
            	if(intensite * U[i] <= paramslambda0 * (1 + paramsalpha * sin(4 * M_PI * Tn[i])))
            	{
                  	L[lenA] = Tn[i];
	                  lenA = lenA+1;
                  }
            }
            std::sort(L.begin(), L.end());
            return(L);
	}
	")
  
  nb=0
  for (i in 1:n.simul){
    a=t_pluies_rcpp(params$lambda0,params$lambda0*(1+params$alpha),1,params$alpha)
    b=rintensite(length(a),params)
    nb=nb+hauteurRcpp(a, b, params$h0, params$v)
  }
  result <- list(p = nb / n.simul, demi.largeur = 1.96 * sqrt((nb / n.simul) * (1 - nb / n.simul) / n.simul))
  return(result)
}

