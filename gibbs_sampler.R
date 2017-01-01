
#example of using algorithm for gibbs sampler with metropolis hastings random walk
# two unknown parameters. theta is sampled using an function that samples from a gamme distribution, which is the full conditional distribution.
#  nu is sampled using a metroplis-hastings random walk step to approximate a random sampling distribution for nu.
set.seed(1212)
data<-read.table("SET YOUR DATA PATH")


#gibbs sampler function
gibbs5 <- function(n.sims, x, burnin, thin) 
{
 
  #initialize parameters of the gibbs sampler
  n<-length(x)
  theta.draws <- c()   # initialize vector that will store draws from the full conditional of rho
  nu.draws<- c() # initialize vector that will store draws from theh full conditional of rho
  nu.cur<-1
  
  theta.update <- function(nu,x) { # updates theta using the full conditional distribution of theta

      shape = n*nu+2
      rate = 2 + sum(x)
      return(rgamma(1,shape = shape,  rate = rate))   #sample from the full conditional distribution of theta parameter.

    }
  }
  
  nu.update <- function(nu.cur, theta,x) { # updates nu using the MH-RW
    mh.rw(nu.cur,theta,x)
  }

  for (i in 1:n.sims) {  # simulates and calls update functions to simulate theta and nu
    theta.cur <- theta.update(nu.cur,x)
    nu.cur <- nu.update(nu.cur, theta.cur, x)
    if (i > burnin & (i - burnin)%%thin == 0) {  # applys burn-in and thining to the simulated data
      theta.draws[(i - burnin)/thin] <- theta.cur
      nu.draws[(i - burnin)/thin] <- nu.cur
    }
  }
  
  sims <- cbind(theta.draws, nu.draws)
  post.mean <- cbind(mean(sims[,1]),mean(sims[,2]))
  return (sims)
}


#metropolis hastings random walk function
mh.rw<-function(nu,theta,x){
  
  n<-length(x)
  log.f.nu <- function(nu){ #log of full conditional of nu
    (nu*n)*log(theta) -nu + 2*log(nu) - n*log(gamma(nu)) + (nu-1)*sum(log(x)) 
  }
  jac<-function(nu){
    exp(nu) # jacobian
  }
  
  nu.pro <-  rnorm(1, mean=log(nu), sd=1)
  
  p.pro <- log.f.nu(jac(nu.pro)) + nu.pro
  p.cur <- log.f.nu(nu) + log(nu)
  
  alpha <- p.pro-p.cur   
  
  if(log(runif(1)) < alpha) #acceptance rule
  {
    nu <- jac(nu.pro)
  }
  
  return(nu)
}

s<-gibbs5(1000,data,100,1)
cbind(mean(s[,1]),mean(s[,2]))
#          [theta]     [nu]     post mean
#   [1,] 1.034744 2.569902
quantile(s[,1],c(.025,.975))
#     2.5%     97.5%  for theta
#  0.6441589 1.5877641 

quantile(s[,2],c(.025,.975))
#     2.5%    97.5% 
#  1.697733 3.855781 



# (b)
