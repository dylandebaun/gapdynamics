method = "final"
type= "all"
nsp = 328
library(stats)
source("sourcefunctions.R")
traitvals <- read.csv("wdfactorvalues.csv")
recruittype=""
count_ngs =0
args <- commandArgs(trailingOnly = TRUE)
out=as.character(args[1])
method=as.character(args[2])

#find the minimum density among all things across all years
min_density_old = 10000
for(j in 1:10){
  if(j == 1){
    simtype = ""
  }else{
    simtype = "nongap"
    count_ngs = j - 1
  }
  for(yr in 1:7){
    if(simtype == "nongap"){ ng_num = j-2 }else{ ng_num = ""}#get nongap number (from 0 to 9)
    year = (yr-1)*5
    
    title = paste0(out,"/",method,simtype,ng_num,recruittype,"expectedin",type, "at",year,"yrs.txt")
    totalcompositionactual <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")
    densityactual = density(totalcompositionactual) #figure out density of the actual gap for subsampling
    
    if(densityactual < min_density_old){
      min_density_old = densityactual
    }
  }
}

#rarefy to minimum sample and make caclulations
for(yr in 1:7){
  year = (yr-1)*5
  dist_act <- matrix(0, nrow=10, ncol=17)
  colnames(dist_act)<-c("rich","even","H","expH","wd","factor","SM","RGR","wdvar","factorvar","SMvar","RGRvar","wdcov","factorcov","SMcov","RGRcov","density")
  
  min_density = 10000
  for(j in 1:10){
    if(j == 1){
      simtype = ""
    }else{
      simtype = "nongap"
      count_ngs = j - 1
    }
    if(simtype == "nongap"){ ng_num = j-2 }else{ ng_num = ""}#get nongap number (from 0 to 9)
    year = (yr-1)*5
    
    title = paste0(out,"/",method,simtype,ng_num,recruittype,"expectedin",type, "at",year,"yrs.txt")
    totalcompositionactual <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")
    densityactual = density(totalcompositionactual) #figure out density of the actual gap for subsampling
    
    if(densityactual < min_density){min_density = densityactual}
  }
  
  for(j in 1:10){
    if(j == 1){
      simtype = ""
    }else{
      simtype = "nongap"
      count_ngs = j - 1
    }
    if(simtype == "nongap"){ ng_num = j-2 }else{ ng_num = ""}#get nongap number (from 0 to 9)
    year = (yr-1)*5
    
    title = paste0(out,"/",method,simtype,ng_num,recruittype,"expectedin",type, "at",year,"yrs.txt")
    totalcompositionactual <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")
    densityactual = density(totalcompositionactual) #figure out density of the actual gap for subsampling
    
    #rarefy it 100 times
    for(y in 1:100){
      actcomp = totalcompositionactual
      testSAD = actcomp
      inds <-c(0)
      compgroup <-c(inds)
      if(isTRUE(min_density < densityactual)){
        #get the difference in number of individuals
        x = densityactual-min_density
        if(y == 1){print(x)}
        if(x != 0){
          #expand the composition vector into a group of individuals
          for(d in 1:dim(testSAD)[1]){
            if(testSAD$n[d] != 0){
              for(z in 1:testSAD$n[d]){
                compgroup <- rbind(compgroup, c(testSAD$sp[d]))
              }
            }
          }
          compgroup <- compgroup[-c(1),]
          #randomly remove the individuals
          rm =  sample(1:length(compgroup),x,replace = F)
          compgroup <- compgroup[-c(rm)]
          #rebuild the SAD
          actcomp <- createSAD(compgroup, nsp)
        }
      }
      
      s=j
      dist_act[s,"rich"] = rich(actcomp) + dist_act[s,"rich"]
      dist_act[s,"even"] = E(actcomp) + dist_act[s,"even"]
      dist_act[s,"H"] = shannon(actcomp)+ dist_act[s,"H"]
      dist_act[s,"expH"] = exp(shannon(actcomp)) + dist_act[s,"expH"]
      dist_act[s,"wd"] = trait(traitvals$altwd, actcomp) + dist_act[s,"wd"]
      dist_act[s,"factor"] = trait(traitvals$growthsurv, actcomp)  + dist_act[s,"factor"]
      dist_act[s,"RGR"] = trait(traitvals$RGR, actcomp) + dist_act[s,"RGR"]
      dist_act[s,"SM"] = trait(traitvals$sm, actcomp)  + dist_act[s,"SM"]
      dist_act[s,"wdcov"] = traitcov(traitvals$altwd, actcomp) + dist_act[s,"wdcov"]
      dist_act[s,"factorcov"] = traitcov(traitvals$growthsurv, actcomp)  + dist_act[s,"factorcov"]
      dist_act[s,"RGRcov"] = traitcov(traitvals$RGR, actcomp) + dist_act[s,"RGRcov"]
      dist_act[s,"SMcov"] = traitcov(traitvals$sm, actcomp)  + dist_act[s,"SMcov"]
      dist_act[s,"wdvar"] = traitvar(traitvals$altwd, actcomp) + dist_act[s,"wdvar"]
      dist_act[s,"factorvar"] = traitvar(traitvals$growthsurv, actcomp)  + dist_act[s,"factorvar"]
      dist_act[s,"RGRvar"] = traitvar(traitvals$RGR, actcomp) + dist_act[s,"RGRvar"]
      dist_act[s,"SMvar"] = traitvar(traitvals$sm, actcomp)  + dist_act[s,"SMvar"]
      if(j == 1){
        dist_act[s,"density"] = min_density  + dist_act[s,"density"]
      }else{
        dist_act[s,"density"] = min_density_old  + dist_act[s,"density"]
      }
      title = paste0("results/trash",type,method,year,simtype,recruittype,"distresults_totalssnongap.csv")
      write.csv(dist_act,title,row.names = F)
    }
  }
  dist_act = dist_act/100
  title = paste0("results/",type,method,year,simtype,recruittype,"distresults_peryearssnongap.csv")
  write.csv(dist_act,title,row.names = F)
}
