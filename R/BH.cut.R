BH.cut<-function(pv, q)
{ 
  # the input 
    # pv the set of the two-sided p-values
    # q the desired FDR level
  # the output is a list with
    # the first element (st.pv) the sorted p-values
    # the second element (k) the number of hypothesis to be rejected
    # the third element (pk) the threshold for the p-values
    # the fourth element (reject) the set of rejected hypotheses
    # the fifth element (accept) the set of accepted hypotheses
  m=length(pv)
  st.pv<-sort(pv)   
  pvi<-st.pv/1:m
  if( sum(pvi<=(q/m))==0)
    {
      reject=0
      accept=0
      k=0
      pk=0
    }else{
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    reject<-which(pv<=pk)
    accept<-which(pv>pk)
  }
  y<-list(spv=st.pv, nr=k, threshold=pk, re=reject, ac=accept)
  return (y)
}
