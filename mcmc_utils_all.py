import numpy as np
import temp_functions_all as tfa

#PLOTTING STUFF
import matplotlib.pyplot as plt
import scipy.stats as stats
'''
## This file contains utility functions for manipulating mcmc samps,
## and calculating/plotting predictions

library(IDPmisc)
library('rjags')

## first some functions to put the mcmc samples into a properly
## labeled data frame, depending on the functional form

# Added by Daniel W, for the linear, columns = inter, slope, sigma
make.linear.samps<-function(coda.samps, nchains=2, samp.lims=c(151, 5000), sig=TRUE){
  n.inter<-slope<-sigma<-NULL
  l1<-samp.lims[1]
  l2<-samp.lims[2]
  for(i in 1:nchains){
    n.inter<-c(n.inter, coda.samps[[i]][l1:l2,1])
    slope<-c(slope, coda.samps[[i]][l1:l2,2])
    if(sig) sigma<-c(sigma, coda.samps[[i]][l1:l2,3])
  }
  if(sig){
    samps<-data.frame(matrix(c(n.inter, slope, sigma), ncol=3, byrow=FALSE))
    names(samps)<-c("n.inter", "slope", "sigma")
  }
  else{
    samps<-data.frame(matrix(c(n.inter, slope), ncol=2, byrow=FALSE))
    names(samps)<-c("n.inter", "slope", "sigma")
    
  }
  return(samps)
}

make.nlinear.samps<-function(coda.samps, nchains=2, samp.lims=c(151, 5000), sig=TRUE){
  inter<-n.slope<-sigma<-NULL
  l1<-samp.lims[1]
  l2<-samp.lims[2]
  for(i in 1:nchains){
    inter<-c(inter, coda.samps[[i]][l1:l2,1])
    n.slope<-c(n.slope, coda.samps[[i]][l1:l2,2])
    if(sig) sigma<-c(sigma, coda.samps[[i]][l1:l2,3])
  }
  if(sig){
    samps<-data.frame(matrix(c(inter, n.slope, sigma), ncol=3, byrow=FALSE))
    names(samps)<-c("inter", "n.slope", "sigma")
  }
  else{
    samps<-data.frame(matrix(c(inter, n.slope), ncol=2, byrow=FALSE))
    names(samps)<-c("inter", "n.slope", "sigma")
    
  }
  return(samps)
}


##for the briere bit, columns = T0, Tm, c, sigma
make.briere.samps<-function(coda.samps, nchains=2, samp.lims=c(151, 5000), sig=TRUE){
  T0<-Tm<-cc<-sigma<-NULL
  l1<-samp.lims[1]
  l2<-samp.lims[2]
  for(i in 1:nchains){
    T0<-c(T0, coda.samps[[i]][l1:l2,1])
    Tm<-c(Tm, coda.samps[[i]][l1:l2,2])
    cc<-c(cc, coda.samps[[i]][l1:l2,3])
    if(sig) sigma<-c(sigma, coda.samps[[i]][l1:l2,4])
  }
  if(sig){
    samps<-data.frame(matrix(c(T0, Tm, cc, sigma), ncol=4, byrow=FALSE))
    names(samps)<-c("T0", "Tm", "c", "sigma")
  }
  else{
    samps<-data.frame(matrix(c(T0, Tm, cc), ncol=3, byrow=FALSE))
    names(samps)<-c("T0", "Tm", "c")
  }

  return(samps)
}

##for the quadratic, columns = inter, quad, sigma, slope
make.quad.samps<-function(coda.samps, nchains=2, samp.lims=c(151, 5000), sig=TRUE){
  T0<-Tm<-qd<-sigma<-NULL
  l1<-samp.lims[1]
  l2<-samp.lims[2]
  for(i in 1:nchains){
    T0<-c(T0, coda.samps[[i]][l1:l2,1])
    Tm<-c(Tm, coda.samps[[i]][l1:l2,2])
    qd<-c(qd, coda.samps[[i]][l1:l2,3])
    if(sig) sigma<-c(sigma, coda.samps[[i]][l1:l2,4])
  }
  if(sig){
    samps<-data.frame(matrix(c(T0, Tm, qd, sigma), ncol=4, byrow=FALSE))
  names(samps)<-c("T0", "Tm", "qd", "sigma")
  }
  else{
    samps<-data.frame(matrix(c(T0, Tm, qd), ncol=3, byrow=FALSE))
    names(samps)<-c("T0", "Tm", "qd")

  }
  return(samps)
}

make.pos.quad.samps<-function(coda.samps, nchains=2, samp.lims=c(151, 5000)){
  inter<-n.slope<-qd<-tau<-NULL
  l1<-samp.lims[1]
  l2<-samp.lims[2]
  for(i in 1:nchains){
    inter<-c(inter, coda.samps[[i]][l1:l2,1])
    n.slope<-c(n.slope, coda.samps[[i]][l1:l2,2])
    qd<-c(qd, coda.samps[[i]][l1:l2,3])
    tau<-c(tau, coda.samps[[i]][l1:l2,4])
  }
  samps<-data.frame(matrix(c(inter, n.slope, qd, tau), ncol=4, byrow=FALSE))
  names(samps)<-c("inter", "n.slope", "qd", "tau")
  return(samps)
}


plot.hists<-function(samps, my.par=c(2,2), n.hists=4, priors=NA){
  par(mfrow=my.par, bty="n")
  for(i in 1:n.hists){
    hist(samps[,i], main=names(samps)[i], xlab="", freq=FALSE)
    if(!is.na(priors)){
      x<-c(seq(0, 1, by=0.00001), seq(1, 4000, by=0.1))
      h<-priors$hyper[,i]
      if(priors$fun[i]=="beta") lines(x, dbeta(x, shape1=as.numeric(h[1]), shape2=as.numeric(h[2])), col=2)
      if(priors$fun[i]=="gamma") lines(x, dgamma(x, shape=as.numeric(h[1]), rate=as.numeric(h[2])), col=2)
      if(priors$fun[i]=="uniform") lines(x, dunif(x, min=as.numeric(h[1]), max=as.numeric(h[2])), col=2)
      if(priors$fun[i]=="exp") lines(x, dexp(x, rate=as.numeric(h[1])), col=2)
      if(priors$fun[i]=="min.gamma") lines(h[3]+x, dgamma(x, shape=as.numeric(h[1]), rate=as.numeric(h[2])), col=2)
      if(priors$fun[i]=="normal") lines(x, dnorm(x, mean=as.numeric(h[1]), sd=as.numeric(h[2])), col=2)
      
    }
  }
}
    

'''



def make_sims_temp_resp(sim, samps, Temps, thin_factor):
  #thin it out by the thin factor

  samps= samps[:][0::thin_factor]

  #create the output data
  data_sim= np.zeros((len(Temps),len(samps)))

  for i in range(0,len(samps)-1):

    if(sim=="briere"): 
      c = samps['c'][i]
      Tm = samps['Tm'][i]
      T0 = samps['T0'][i]
      #select temperatures within a specific range


      data_sim[:,i] = tfa.briere(Temps,c,Tm,T0)


    if(sim=="briere_trunc"): 
      c = samps['c'][i]
      Tm = samps['Tm'][i]
      T0 = samps['T0'][i]
      #select temperatures within a specific range


      data_sim[:,i] = tfa.briere_trunc(Temps,c,Tm,T0)
   

    if(sim=="quad_trunc"): 
      c = samps['c'][i]
      Tm = samps['Tm'][i]
      T0 = samps['T0'][i]
      #select temperatures within a specific range


      data_sim[:,i] = tfa.quad_trunc(Temps,c,Tm,T0)
     
    if(sim=="quad"): 
      c = samps['c'][i]
      Tm = samps['Tm'][i]
      T0 = samps['T0'][i]
      #select temperatures within a specific range


      data_sim[:,i] = tfa.quad_2(Temps,c,Tm,T0)
 
  return data_sim



'''
make.sims.temp.resp<-function(sim="quad", samps, Temps, thinned, p.name="PDR", trunc.num=0){

  out<-data.sim<-NULL
  out<-list()
  data.sim<-matrix(NA, nrow=length(Temps), ncol=length(thinned))
  for(i in 1:length(thinned)){

    if(sim=="briere"){
      c<-as.numeric(samps$c[thinned[i]])
      Tm<-as.numeric(samps$Tm[thinned[i]])
      T0<-as.numeric(samps$T0[thinned[i]])
      w0<-which(Temps<=T0)
      wM<-which(Temps>=Tm)
      data.sim[,i]<-briere(Temps, c, Tm, T0)
      data.sim[c(w0,wM),i]<-0
    }

    if(sim=="briere.trunc"){
      c<-as.numeric(samps$c[thinned[i]])
      Tm<-as.numeric(samps$Tm[thinned[i]])
      T0<-as.numeric(samps$T0[thinned[i]])
      w0<-which(Temps<=T0)
      wM<-which(Temps>=Tm)
      data.sim[,i]<-briere.trunc(Temps, c, Tm, T0)
      data.sim[c(w0,wM),i]<-0
    }
    
    
    
    if(sim=="quad"){
      T0<-as.numeric(samps$T0[i])
      Tm<-as.numeric(samps$Tm[i])
      qd<-as.numeric(samps$qd[i])
      data.sim[,i]<-quad.2(Temps, T0=T0, Tm=Tm, qd=qd)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }


    if(sim=="quad.pos.trunc"){    # added by EAM on 8/31/15
      inter<-as.numeric(samps$inter[i])
      n.slope<- as.numeric(samps$n.slope[i])
      qd<-as.numeric(samps$qd[i])
      data.sim[,i]<-quad(Temps, inter=inter, n.slope=n.slope, qd=qd)
      w<-which(data.sim[,i]<trunc.num)
      data.sim[w,i]<-trunc.num
    }
    
      if(sim=="quad.trunc"){
      T0<-as.numeric(samps$T0[i])
      Tm<-as.numeric(samps$Tm[i])
      qd<-as.numeric(samps$qd[i])
      data.sim[,i]<-quad.trunc(Temps, T0=T0, Tm=Tm, qd=qd)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }
    
        if(sim=="linear"){
          n.inter<-as.numeric(samps$n.inter[i])
          slope<-as.numeric(samps$slope[i])
          data.sim[,i]<-linear(Temps, inter=-n.inter, slope=slope)
          w<-which(data.sim[,i]<0)
          data.sim[w,i]<-0
        }
    
    if(sim=="nlinear"){
      inter<-as.numeric(samps$inter[i])
      n.slope<-as.numeric(samps$n.slope[i])
      data.sim[,i]<-linear(Temps, inter=inter, slope=-n.slope)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }
    
    
  }
  out<-list(param = p.name,
            T   = Temps,
            fits = data.sim)
  return(out)
}

'''

'''




make.all.temp.resp<-function(a.p, n.p, thinned, Temps){
  out<-list()

  for(i in 1:n.p){
    out[[i]]<-make.sims.temp.resp(sim=all.params[[i]]$func,
                                  samps=all.params[[i]]$samps,
                                  Temps, thinned=thinned,
                                  p.name=all.params[[i]]$param)
  }

  return(out)

}
'''



def temp_sim_quants(sim_data,Temps):
  
  p= np.array([2.5,97.5])
  q= np.zeros((len(p), len(Temps)))

  for i in range(0, len(Temps)-1):
    q[:,i]=np.percentile(sim_data[i,:],p)

  return q 


'''
temp.sim.quants<-function(sim.data, l.temps, byCol=FALSE,
                     probs=c(0.025, 0.975)){

  q<-matrix(NA, nrow=length(probs), ncol=l.temps)
  if(byCol) for(i in 1:l.temps) q[,i]<-quantile(sim.data[,i], probs, na.rm=TRUE)
  else for(i in 1:l.temps) q[,i]<-quantile(sim.data[i,], probs, na.rm=TRUE)
  
  return(q)
}

'''






'''

make.all.r0<-function(){


}
'''






def add_sim_lines(Temps,sim_data, q, Y, T, figure_count): 

  figure_count= figure_count+1
  plt.figure(figure_count)
  #plot the quantiles
  plt.plot(Temps,q[0,:],'r--')
  plt.plot(Temps,q[1,:],'r--')

  #plot the mean
  plt.plot(Temps,np.mean(sim_data, axis=1),'k')

  #plot the data
  plt.plot(T, Y,'bo')

  plt.axis([0, 45, np.min(q[0,:])-(np.min(q[0,:])/100), np.max(q[1,:])+(np.min(q[0,:])/100) ])
  plt.xlabel('Temperature (C)')
 
  plt.show()
  return figure_count


'''
add.sim.lines<-function(t, sim.data=NULL, q=NULL, q2=NULL, mycol="blue", lwd=1){

  if(!is.null(sim.data)){
    ## plot the predicted mean
    if(!is.null(dim(sim.data)[1])) lines(t, rowMeans(sim.data, na.rm=TRUE),
                                         col=mycol, lwd=lwd)
    else lines(t, sim.data, col=mycol, lwd=lwd)
  }
  if(!is.null(q)){
    ## plot the predicted quantiles for the mean
    lines(t, q[1,], col=mycol, lty="dashed", lwd=(lwd+1))
    lines(t, q[2,], col=mycol, lty="dashed", lwd=(lwd+1))
  }
  if(!is.null(q2)){
    ## plot the predicted quantiles for the mean
    lines(t, q2[1,], col=mycol, lty="dotted", lwd=(lwd+2))
    lines(t, q2[2,], col=mycol, lty="dotted", lwd=(lwd+2))
  }
  
}
'''



#SOME EXTRA FUNCTIONS 

def gamma_fit(data, plot_boolean):
  
  gamma= stats.gamma
  y= data
  x = np.linspace(0, y.max(), 100)
    # fit
  param = gamma.fit(y,floc=0 )
  
  #aka do you want there to be plotting
  if(plot_boolean):
    pdf_fitted = gamma.pdf(x, *param)
    plt.plot(x,pdf_fitted, color='r')
    # plot the histogram
    plt.hist(y, normed=True, bins=30)

  else: 
    return np.array([param[0],param[2]]) #only export the scale and shape, not loc


'''
def create_2x2_histograms(data, figure_count):

  figure_count= figure_count+1
  plt.figure(figure_count)

  f, axarr = plt.subplots(2, 2)
  axarr[0, 0].hist(data['T0'][:], 50, normed=1, facecolor='g', alpha=0.75)
  axarr[0, 0].set_title('T0')
  axarr[0, 1].hist(data['Tm'][:], 50, normed=1, facecolor='g', alpha=0.75)
  axarr[0, 1].set_title('Tm')
  axarr[1, 0].hist(data['c'][:], 50, normed=1, facecolor='g', alpha=0.75)
  axarr[1, 0].set_title('c')
  axarr[1, 1].hist(data['tau'][:], 50, normed=1, facecolor='g', alpha=0.75)


  axarr[1, 1].set_title('tau')

  plt.show()

  return figure_count
'''


def create_2x2_histograms(data, figure_count):

  plot_boolean= True
  figure_count= figure_count+1
  plt.figure(figure_count)

  plt.subplot(221)
  gamma_fit(data['T0'][:],plot_boolean)

  plt.subplot(222)
  gamma_fit(data['Tm'][:],plot_boolean)

  plt.subplot(223)
  gamma_fit(data['c'][:],plot_boolean)
  
  plt.subplot(224)
  gamma_fit(data['tau'][:],plot_boolean)

  plt.show()

  return figure_count






