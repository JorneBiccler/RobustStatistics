library(survival)

trimmedCox = function(t,
                  delta,
                  x,
                  alpha = 0.1,
                  poscaso = NULL,
                  D = NULL,
                  kmax = NULL,
                  rseqmax = NULL,
                  n.multistart = NULL) {
  # t: times
  # delta: event indicators
  # x: matrix of covariates
  # poscaso: optional initial solution
  # alpha: trimming level
  # D, kmax, rseqmax, n.multistart: tuning parameters for the fitting algorithm
  
  n = length(t)
  
  if (is.null(kmax)) {
    kmax = 1000
  }
  if (!is.matrix(x)) {
    x = matrix(x)
  }
  if (is.null(rseqmax)) {
    rseqmax = 20
  }
  if (is.null(D)) {
    D = 0.1
  }
  if (is.null(poscaso)) {
    poscaso = sample(n, ceiling(n * (1 - alpha)))
  }
  if (is.null(n.multistart)) {
    n.multistart = 10
  }
  
  max.lik = -Inf
  
  for (starts in 1:n.multistart) {
    if (starts > 1) {
      poscaso = sample(n, ceiling(n * (1 - alpha)))
    }
    
    t.sub = t[poscaso]
    delta.sub = delta[poscaso]
    x.sub = x[poscaso, ]
    
    lt = length(poscaso)
    pos2 = poscaso
    stime = coxph(Surv(t.sub, delta.sub) ~ ., data = data.frame(x.sub))
    
    beta = stime$coefficients
    
    lik = stime$loglik[2]
    likold = -Inf
    k = 1
    bestlik = lik
    b2 = beta
    p2 = pos2
    rseq = 0
    
    while (k < kmax && rseq < rseqmax) {
      if (lik - likold != 0) {
        rseq = 0
      }
      if (lik - likold == 0) {
        rseq = rseq + 1
      }
      
      likold = lik
      
      for (position in 1:lt) {
        Tk = log(k + 1, 2) / (lt * D)
        k = k + 1
        
        candidato = sample((1:n)[-pos2], 1)
        
        aa = coxph(Surv(t[c(pos2[-position], candidato)], delta[c(pos2[-position], candidato)]) ~
                     x[c(pos2[-position], candidato), ])
        
        lik.cand = aa$loglik[2]
        
        pii = min(exp(Tk * (lik.cand - lik)), 1)
        
        if (rbinom(1, 1, pii) == 1)  {
          pos2 = c(pos2[-position], candidato)
          likold = lik
          lik = lik.cand
          beta = aa$coefficients
          if (lik > bestlik) {
            bestlik = lik
            b2 = beta
            p2 = pos2
          }
        }
        
        
      }
      
    }
    
    lik = bestlik
    beta = b2
    pos2 = p2
    
    if (lik > max.lik)
    {
      max.lik = lik
      max.pos2 = pos2
      max.beta = beta
    }
    
  }
  
  return(list(
    lik = max.lik,
    outliers = (1:n)[-max.pos2],
    units = max.pos2,
    beta = max.beta,
    trim = alpha
  ))
}
