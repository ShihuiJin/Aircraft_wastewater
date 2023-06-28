#load case data and stool distribution
load('cases')
load('stool.dist')
load('modelled.cases')
load('urban.pop')

library(dplyr)
#functions for calculating stool probabilities
{
  #stool frequency
  stool.prob.cal=function(flight.time,n.stool,gamma=1,delta.before=0,delta.after=0){
    #gamma (<=1): relative probability that one would use toilet on board
    gamma=max(0,min(gamma,1))
    #turn hour of day to minutes
    flight.time=flight.time+c(delta.after,-delta.before) #account for time when toilet is not allowed to use
    t.temp=round(flight.time*60)
    if(max(t.temp)<=24*60){
      t.range=min(t.temp):max(t.temp)
    }else{
      t.range=c(1:(max(t.temp)-24*60),min(t.temp):(24*60))
    }
    do.call('sum',lapply(1:length(n.stool), function(i)
      min(n.stool.dist.agg[i], gamma*sum(n.stool[[i]][t.range])) ))
    # as.numeric(diff(range(flight.time))>(24/i/gamma))*n.stool.dist.agg[i]+
    #   gamma*as.numeric(diff(range(flight.time))<=(24/i/gamma))*sum(n.stool[[i]][t.range]) ))
    
  }
  #stool frequency: equally distributed during one day
  stool.prob.cal.ave=function(flight.time,n.stool,gamma=1,delta.before=0,delta.after=0){
    gamma=max(0,min(gamma,1))
    flight.time=flight.time+c(delta.after,-delta.before)
    delta=diff(round(range(flight.time)*60)/60)
    do.call('sum',lapply(1:length(n.stool), function(i) min(delta*i/24*gamma, 1)*n.stool.dist.agg[i] ))
  }
}
#functions for calculating the probabilities
{
  #d is the start day in a week, default: Monday
  flight.schedule=function(flight.out,n.t,d=1){
    n.t1=n.t+6
    flight.chart.week=flight.out[,c('Mon','Tue','Wed','Thu','Fri','Sat','Sun')]
    flight.chart.week[is.na(flight.chart.week)]=0
    if(n.t1%%7==0){
      flight.chart=as.matrix(do.call('cbind',lapply(1:(n.t1/7), function(i) flight.chart.week)))
    }else{
      flight.chart=as.matrix(cbind(as.matrix(do.call('cbind',lapply(1:floor(n.t1/7), function(i) flight.chart.week))),
                                   flight.chart.week[,1:(n.t1%%7)]))
    }
    flight.chart[,1:n.t+d-1]
  }
  
  #with known start day of the week
  prob.positive.global=function(airport,tau=1,gamma=1, p.test=0.5,infectious.por=infectious.p,d=1,cut.time=c(0,0)){
    n.t=length(infectious.por)
    idx=which(flight.out$airport==airport)
    duration=flight.out$duration[idx[1]]
    flight.chart=flight.schedule(flight.out,n.t,d)
    1-apply(as.matrix(do.call('rbind',lapply(idx, function(i){
      v=flight.out$n.passengers[i];take.off=flight.out$departure[i];
      if(is.na(take.off)){
        1-tau*(1-dbinom(0,v,stool.prob.cal.ave(c(0,duration),n.stool,gamma,cut.time[1],cut.time[2])*p.test*infectious.por))*flight.chart[i,] 
      }else{
        1-tau*(1-dbinom(0,v,stool.prob.cal(c(0,duration)+take.off,n.stool,gamma,cut.time[1],cut.time[2])*p.test*infectious.por))*flight.chart[i,] 
      }
    }))), 2,prod)
  }
  
  #aggregating the probabilities, for all airports version
  #return: cumulative & new detection probabilities
  prob.mutliple.aiports=function(airport.list,tau,gamma,p.test,infectious.p,d,cut.time){
    #daily prob.
    temp1=1-matrixStats::colProds(as.matrix(do.call('rbind',lapply(unique(flight.out$airport), function(i) 
      1-prob.positive.global(i,tau,gamma,p.test,infectious.p,d,cut.time) ))), na.rm=T) 
    #prob. of detection by day t
    temp2=1-cumprod(1-temp1)
    #prob. of first detection on day t
    temp3=temp1*c(1,1-temp2[1:length(temp2)-1])
    rbind(temp1, temp2, temp3)
  }
  
  prob.by.epicenter=function(epicenter,tau,gamma,p.test,infectious.p,cut.time=c(0,0)){
    if(epicenter=='Wuhan'){
      #known starting day: Monday
      prob.mutliple.aiports(unique(flight.out$airport), tau,gamma,p.test,infectious.p,1,cut.time)
    }else{
      #unknown starting day of the outbreak, can be any day of the week
      temp.prob=as.matrix(do.call('rbind',lapply(1:7, function(d) 
        prob.mutliple.aiports(unique(flight.out$airport), 
                              tau,gamma,p.test,infectious.p,d,cut.time))))
      rbind(colMeans(temp.prob[1+0:6*3,]),
            colMeans(temp.prob[2+0:6*3,]),
            colMeans(temp.prob[1:7*3,]))
    }
  }
}

#calculating the probabilities
#baseline scenario: testing all planes, no change in travel/defecation behavior, testing sensitivity 0.5, no ban from toilet use
baseline=list()
for(c in 1:5)
{
  epicenter=urban.pop$city[c]; pop=urban.pop$pop.m[c]*1e6
  flight.out=as.data.frame(read.csv(paste0(address.pre,'data/',epicenter,'.travel.csv')))
  #airport.list=unique(flight.out$airport)
  if(epicenter!='Wuhan'){
    flight.out=flight.out%>%left_join(aircraft.load,by='aircraft')
    infectious.p=infectious/pop #percentage of infections in the total population
  }else{
    infectious.p=infectious[1:45]/pop #until Jan 22, 2020 only
  }
  n.t=length(infectious.p)
  baseline[[c]]=prob.by.epicenter(epicenter, 1,1, 0.5,infectious.p,c(0,0))
}