# resolveambiguous function updated from Aleks Obradovic 2019 R toolkit that doesn't 
# remove adequately promiscuous clones

resolveambiguous<-function(data,c1,c2,ratio=5){
  normalize(data[,c(c1,c2)])
  c1indices=which(data[,c1]>0 & data[,c2]>0 & data[,c1]>ratio*data[,c2])
  c2indices=which(data[,c1]>0 & data[,c2]>0 & data[,c2]>ratio*data[,c1])
  # ambiindices=which(data[,c1]>0 & data[,c2]>0)
  # ambiindices=setdiff(setdiff(ambiindices,c1indices),c2indices) #this just means data that is > 0 in both columns but under ratio
  data[c1indices,c2]=0
  data[c2indices,c1]=0
  # data[ambiindices,c(c1,c2)]=0
  return(data)
}

normalize <- function(data) {
  nc = ncol(data)
  for (i in 1:nc) {
    data[,i] = data[,i ] / sum(data[,i])
  }
  return(data)
}
