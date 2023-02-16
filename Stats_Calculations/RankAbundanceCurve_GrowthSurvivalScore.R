#For Creating Rank Abundance Curves
#method = "final"
#numsim = 500

args <- commandArgs(trailingOnly = TRUE)
out=as.character(args[1])
method=as.character(args[2])
outdir =as.character(args[3])
method=as.character(args[4])
startsim=as.numeric(args[5])
endsim=as.numeric(args[6])

startsim=startsim+1
endsim=endsim
numsim= endsim-startsim+1

year=10
gap_number=1
j=1
simtype = ""
recruittype = ""
type = "all"

nsp = 328
library(stats)
source("sourcefunctions.R")
traitvals <- read.csv("TraitValues.csv") #what the trait values are for the species
ng_num = ""
count_ngs = j- 1

dist_testRA <- matrix(0, ncol=numsim,nrow=nsp)
dist_actRA <- matrix(0, ncol=numsim*100,nrow=nsp)
average_factor_value_act <-matrix(0, ncol=1,nrow=nsp)
average_factor_value_sim <-matrix(0, ncol=1,nrow=nsp)
countact <-matrix(0, ncol=1,nrow=nsp)
countsim <-matrix(0, ncol=1,nrow=nsp)
title = paste0(out,"/",method,ng_num,simtype,recruittype,"expectedin",type, "at",year,"yrs.txt")
totalcompositionactual <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")
densityactual = density(totalcompositionactual) #figure out density of the actual gap for subsampling

inds<- c(0)
##### rarefy samples######     
count = 0
for(s in 1:numsim){
  #rarefy 100 times to get accurate inference for each sample (averaging rarefied samples)
  for(y in 1:100){
    count = count + 1
    actcomp = totalcompositionactual
    title = paste0(out, "/",method,"simulation",s-1,type,"speciescountyr",year,".txt")
    if(file.exists(title)){
      totalcomposition <- read.delim(title,header=F)
    }else{
      cat("dne: ")
      cat(s-1)
      cat("\n")
    }
    
    #go to the column we want in the data space
    col_num = j + count_ngs
    testSAD<- totalcomposition[,col_num:(col_num+1)]
    colnames(testSAD)<-c("sp","n")
    
    #get the density for rarefying
    densitySAD <- density(testSAD)
    if(densitySAD==0){
      cat("empty: ")
      cat(s-1)
      cat("\n")
    }
    
    #RAREFY
    inds <-c(0)
    compgroup <-c(inds)
    if(isTRUE(densitySAD > densityactual)){
      print("DENSITY IS BIGGER")
      count = 10000000000000
      #get the difference in number of individuals
      #rebuild the SAD
      testSAD <- 1
    }else if(isTRUE(densitySAD < densityactual)){
      x = densityactual-densitySAD
      if(x != 0){
        #expand the composition vector into a group of individuals
        for(d in 1:dim(actcomp)[1]){
          if(actcomp$n[d] != 0){
            for(z in 1:actcomp$n[d]){
              compgroup <- rbind(compgroup, c(actcomp$sp[d]))
            }
          }
        }
        compgroup <- compgroup[-c(1),]
        rm =  sample(1:length(compgroup),x,replace = F)
        compgroup <- compgroup[-c(rm)]
        actcomp <- createSAD(compgroup, nsp)
      }
    }
    
    # #rank abundance curves
    testrankab = makerankab(testSAD)
    actrankab = makerankab(actcomp)
   
    print(count)
    dist_testRA[,s] = testrankab[,2] 
    #print(testrankab)
    dist_actRA[,count] = actrankab[,2] 
    #print(actrankab)

    for(ijk in 1:dim(average_factor_value_act)[1]){
      if(!is.na(traitvals$growthsurv[actrankab[ijk,1]+1]) && (actrankab[ijk,2] != 0)){
  	    countact[ijk,1] = countact[ijk,1] +1
  	average_factor_value_act[ijk,1] = traitvals$growthsurv[actrankab[ijk,1]+1] +  average_factor_value_act[ijk,1] 
      }
      if(!is.na(traitvals$growthsurv[testrankab[ijk,1]+1]) && (testrankab[ijk,2] != 0)){
      countsim[ijk,1] = countsim[ijk,1] +1
      average_factor_value_sim[ijk,1] = traitvals$growthsurv[testrankab[ijk,1]+1] +  average_factor_value_sim[ijk,1]
      }
    }
    #print(average_factor_value_act)
    #print(average_factor_value_sim)
  }
}
x= dist_actRA
y=dist_testRA
write.csv((x), paste0("plots/actualgaprankabsatyr",year, "_nolog.csv"), row.names = F)
write.csv((y), paste0("plots/simulatedgaprankabsatyr",year, "_nolog.csv"), row.names = F)
write.csv(log10(x), paste0("plots/actualgaprankabsatyr",year, ".csv"), row.names = F)
write.csv(log10(y), paste0("plots/simulatedgaprankabsatyr",year, ".csv"), row.names = F)
write.csv(average_factor_value_act[,1]/countact[,1], paste0("plots/actualgaprankabsFACTORatyr",year, ".csv"), row.names = F)
write.csv(average_factor_value_sim[,1]/countsim[,1], paste0("plots/simulatedgaprankabsFACTORatyr",year, ".csv"), row.names = F)
