
#----------------------- G-function -----------------------#
Gfunc=function(x,r){
  if(r==0){return(x)}
  else{return(log(1+r*x)/r)}
}

#----------------------- Derivative of G-function -----------------------#
Gdfunc=function(x,r){
  if(r==0){return(1)}
  else{return(1/(1+r*x))}
}

#----------------------- WVIC -----------------------#
s_cal_adj=function(x,jumppts,cumhaz,regcoef,r){
  Ltime=x[1]; Rtime=x[2]
  expterm=exp(t(x[-c(1,2)])%*%regcoef)
  cumhaz_L=ifelse(Ltime==0,0,cumhaz[tail(which(Ltime>=jumppts),1)])
  cumhaz_R=ifelse(Rtime==Inf,Inf,cumhaz[tail(which(Rtime>=jumppts),1)])
  chLexp=cumhaz_L*expterm
  chRexp=cumhaz_R*expterm
  if(Ltime==0)
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp)/(1-exp(-Gfunc(chRexp,r))))
  else if(Rtime==Inf)
    return(-Gdfunc(chLexp,r)*chLexp)
  else
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp-exp(-Gfunc(chLexp,r))*Gdfunc(chLexp,r)*chLexp)/(exp(-Gfunc(chLexp,r))-exp(-Gfunc(chRexp,r))))
}

s_cal_noadj=function(x,jumppts,cumhaz,r){
  Ltime=x[1]; Rtime=x[2]
  expterm=1
  cumhaz_L=ifelse(Ltime==0,0,cumhaz[tail(which(Ltime>=jumppts),1)])
  cumhaz_R=ifelse(Rtime==Inf,Inf,cumhaz[tail(which(Rtime>=jumppts),1)])
  chLexp=cumhaz_L*expterm
  chRexp=cumhaz_R*expterm
  if(Ltime==0)
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp)/(1-exp(-Gfunc(chRexp,r))))
  else if(Rtime==Inf)
    return(-Gdfunc(chLexp,r)*chLexp)
  else
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp-exp(-Gfunc(chLexp,r))*Gdfunc(chLexp,r)*chLexp)/(exp(-Gfunc(chLexp,r))-exp(-Gfunc(chRexp,r))))
}


#----------------------- WVICLT -----------------------#
s_cal_adj_trunc=function(x,jumppts,cumhaz,regcoef,r){
  V0time=x[1]; Ltime=x[2]; Rtime=x[3]
  expterm=exp(t(x[-c(1,2,3)])%*%regcoef)
  cumhaz_V0toL=ifelse(Ltime==V0time,0,cumhaz[tail(which(Ltime>=jumppts),1)]-cumhaz[tail(which(V0time>=jumppts),1)])
  cumhaz_V0toR=ifelse(Rtime==Inf,Inf,cumhaz[tail(which(Rtime>=jumppts),1)]-cumhaz[tail(which(V0time>=jumppts),1)])
  chLexp=cumhaz_V0toL*expterm
  chRexp=cumhaz_V0toR*expterm
  if(Ltime==V0time)
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp)/(1-exp(-Gfunc(chRexp,r))))
  else if(Rtime==Inf)
    return(-Gdfunc(chLexp,r)*chLexp)
  else
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp-exp(-Gfunc(chLexp,r))*Gdfunc(chLexp,r)*chLexp)/(exp(-Gfunc(chLexp,r))-exp(-Gfunc(chRexp,r))))
}

s_cal_noadj_trunc=function(x,jumppts,cumhaz,r){
  V0time=x[1]; Ltime=x[2]; Rtime=x[3]
  expterm=1
  cumhaz_V0toL=ifelse(Ltime==V0time,0,cumhaz[tail(which(Ltime>=jumppts),1)]-cumhaz[tail(which(V0time>=jumppts),1)])
  cumhaz_V0toR=ifelse(Rtime==Inf,Inf,cumhaz[tail(which(Rtime>=jumppts),1)]-cumhaz[tail(which(V0time>=jumppts),1)])
  chLexp=cumhaz_V0toL*expterm
  chRexp=cumhaz_V0toR*expterm
  if(Ltime==V0time)
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp)/(1-exp(-Gfunc(chRexp,r))))
  else if(Rtime==Inf)
    return(-Gdfunc(chLexp,r)*chLexp)
  else
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp-exp(-Gfunc(chLexp,r))*Gdfunc(chLexp,r)*chLexp)/(exp(-Gfunc(chLexp,r))-exp(-Gfunc(chRexp,r))))
}


