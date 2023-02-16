#GAP VS NONGAP - TOTAL NONSUBSAMPLING
method = "final"
type= "all"
nsp = 328
library(stats)
source("sourcefunctions.R")
traitvals <- read.csv("TraitValues.csv")
recruittype=""
count_ngs =0

args <- commandArgs(trailingOnly = TRUE)
out=as.character(args[1])
method=as.character(args[2])
#find the minimum density among all things across all years
min_density = 10000
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
      
    title = paste0(out,"/",method,ng_num,simtype,recruittype,"expectedin",type, "at",year,"yrs.txt")
    totalcompositionactual <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")
    densityactual = density(totalcompositionactual) #figure out density of the actual gap for subsampling
    
    if(densityactual < min_density){
      min_density = densityactual
    }
  }
}

#rarefy to minimum sample and make caclulations
for(yr in 1:7){
  year = (yr-1)*5
dist_act <- matrix(0, ncol=37,nrow=1)
        colnames(dist_act)<-c("ks","ksunsub","pp","ppunsub","bc","bcunsub","rich","even","H","expH","wd","factor","SM","RGR","wdvar","factorvar","SMvar","RGRvar","wdcov","factorcov","SMcov","RGRcov","density","lma","lmavar","ldmc","ldmcvar","P_shade","P_shadevar","N_shade","N_shadevar","P_sun","P_sunvar","N_sun","N_sunvar","pca2","pca2var")
#  dist_act <- matrix(0, nrow=10, ncol=18)
 # colnames(dist_act)<-c("bc","rich","even","H","expH","wd","factor","SM","RGR","wdvar","factorvar","SMvar","RGRvar","wdcov","factorcov","SMcov","RGRcov","density")
  
  for(j in 1:10){
    if(j == 1){
      simtype = ""
    }else{
      simtype = "nongap"
      count_ngs = j - 1
    }
    if(simtype == "nongap"){ ng_num = j-2 }else{ ng_num = ""}#get nongap number (from 0 to 9)
    year = (yr-1)*5
    
    title = paste0(out, "/",method,ng_num,simtype,recruittype,"expectedin",type, "at",year,"yrs.txt")
    totalcompositionactual <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")
    densityactual = density(totalcompositionactual) #figure out density of the actual gap for subsampling
    
    title = paste0(out, "/",method,ng_num,simtype,recruittype,"expectedin",type, "at",0,"yrs.txt")
    yr0comp <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")

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
      dist_act[s,"bc"] = braycurtis(totalcompositionactual,yr0comp) + dist_act[s,"bc"]
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
      dist_act[s,"density"] = densityactual  + dist_act[s,"density"]
              dist_act[s,"lmavar"] = traitvar(traitvals$lma, actcomp)  + dist_act[s,"lmavar"]
        dist_act[s,"ldmcvar"] = traitvar(traitvals$ldmc, actcomp) + dist_act[s,"ldmcvar"]
        dist_act[s,"N_shadevar"] = traitvar(traitvals$N_shade, actcomp)  + dist_act[s,"N_shadevar"]
        dist_act[s,"P_shadevar"] = traitvar(traitvals$P_shade, actcomp)  + dist_act[s,"P_shadevar"]
        dist_act[s,"N_sunvar"] = traitvar(traitvals$N_sun, actcomp)  + dist_act[s,"N_sunvar"]
        dist_act[s,"P_sunvar"] = traitvar(traitvals$P_sun, actcomp)  + dist_act[s,"P_sunvar"]
        dist_act[s,"pca2var"] = traitvar(traitvals$pca2, actcomp)  + dist_act[s,"pca2var"]
        dist_act[s,"lma"] = trait(traitvals$lma, actcomp) + dist_act[s,"lma"]
        dist_act[s,"ldmc"] = trait(traitvals$ldmc, actcomp)  + dist_act[s,"ldmc"]
        dist_act[s,"N_shade"] = trait(traitvals$N_shade, actcomp) + dist_act[s,"N_shade"]
        dist_act[s,"P_shade"] = trait(traitvals$P_shade, actcomp)  + dist_act[s,"P_shade"]
              dist_act[s,"N_sun"] = trait(traitvals$N_sun, actcomp) + dist_act[s,"N_sun"]
        dist_act[s,"P_sun"] = trait(traitvals$P_sun, actcomp)  + dist_act[s,"P_sun"]
        dist_act[s,"pca2"] = trait(traitvals$pca2, actcomp) + dist_act[s,"pca2"]

      title = paste0("results/trash",type,method,year,simtype,recruittype,"distresults_totalssnongap.csv")
      write.csv(dist_act,title,row.names = F)
    }
  }
  dist_act = dist_act/100
  title = paste0("results/",type,method,year,simtype,recruittype,"distresults_totalssnongap.csv")
  write.csv(dist_act,title,row.names = F)
}

