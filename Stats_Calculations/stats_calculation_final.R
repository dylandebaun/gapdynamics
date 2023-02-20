#method = "final"
#numsim = 500
subsample_to_find_avg_distribution = FALSE
nsp = 328
count_ngs =0
library(stats)
source("sourcefunctions.R")
traitvals <- read.csv("TraitValues.csv") #what the trait values are for the species
type="all"

args <- commandArgs(trailingOnly = TRUE)
year=as.numeric(args[1])
gap_number=as.numeric(args[2])
outdir =as.character(args[3])
method=as.character(args[4])
startsim=as.numeric(args[5])
endsim=as.numeric(args[6])

startsim=startsim+1
endsim=endsim
numsim= endsim-startsim+1

for(j in gap_number:gap_number){ #j looks through each gap and nongap community
  if(j == 1){
    simtype = ""
  }else{
    simtype = "nongap"
    count_ngs = j- 1
  }
      if(year == 0){
      	if(simtype == "nongap"){ ng_num = j-2 }else{ ng_num = ""}#get nongap number (from 0 to 9)
        title = paste0(outdir,"/",method,ng_num,simtype,"expectedinall", "at",year,"yrs.txt")
        totalcompositionactual <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")
        densityactual = density(totalcompositionactual) #figure out density of the actual gap for subsampling
        dist_act <- matrix(0, ncol=40,nrow=numsim)
               colnames(dist_act)<-c("ks","ksunsub","pp","ppunsub","bc","bcunsub","rich","even","H","expH","wd","factor","SM","RGR","wdvar","factorvar","SMvar","RGRvar","wdcov","factorcov","SMcov","RGRcov","density","lma","lmavar","lma_unpub","lmavar_unpub","ldmc","ldmcvar","P_shade","P_shadevar","N_shade","N_shadevar","P_sun","P_sunvar","N_sun","N_sunvar","pca2","pca2var","wind")
        actcomp = totalcompositionactual

	      s=1
        dist_act[s,"rich"] = rich(actcomp) + dist_act[s,"rich"]
        dist_act[s,"even"] = E(actcomp) + dist_act[s,"even"]
        dist_act[s,"expH"] = exp(shannon(actcomp))+ dist_act[s,"expH"]
	      dist_act[s,"H"] = shannon(actcomp)+ dist_act[s,"H"]
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
        dist_act[s,"wind"] = trait(traitvals$wind, actcomp) + dist_act[s,"wind"]

        colnames(dist_act)<-c("KS_CDF_test","KS_CDF_test_unsamplesavg","Chi_Sq_Preston_Plot","Chi_Sq_Preston_Plot_unsamplesavg","Bray_Curtis_Dissimilarity","Bray_Curtis_Dissimilarity_unsamplesavg","Richness","Evenness","Shannons_H_Diversity","Effective_Species_expH","Wood_Density","Factor_Score","Seed_Mass","RGR","WoddDensity_Variance","Factor_variance","SMvariance","RGRvariance","WDcovar","Factorcovar","SMcovar","RGRcovar","density","LMA_published_DISC","LMA_published_DISC_variance","LMA_UNpublished_DISC_variance","LMA_UNpublished_DISC_variance","LDMC_unpublished","LDMC_unpublished_variance","Pshade","Pshadevar","Nshade","Nshadevar","Phosphorous_sun","Phosphorous_sun_var","Nitrogen_sun","Nitrogen_sun_variance","pca2","pca2var","WindDispersal")

        dist_sim = dist_act
        
        title = paste0("results/",method,year,ng_num,simtype,"distresultsactual.csv")
        write.csv(dist_act,title,row.names = F)
        title = paste0("results/",method,year,ng_num,simtype,"distresultssimulated.csv")
        write.csv(dist_sim,title,row.names = F)
      }else{
        if(simtype == "nongap"){ ng_num = j-2 }else{ ng_num = ""} #get nongap number (from 0 to 9)
        title = paste0(outdir,"/",method,ng_num,simtype,"expectedinall", "at",year,"yrs.txt")
        totalcompositionactual <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")
        densityactual = density(totalcompositionactual) #figure out density of the observed gap for subsampling
  
        #create lists to hold the compositions and cdfs
        actcomp1 <- list() #observed composition
        testSAD1 <- list() #expected composition
        actcdf1 <- list()
        testcdf1 <- list()
        
        #create preston plot,cdf, avreage composition data frames
        ppcount <-as.data.frame(matrix(c(c("1","2 to 3","4 to 7","8 to 15","16-31","32+"), rep(0,6)),nrow = 6, ncol =2))
        colnames(ppcount) <- c("category","n")
        ppcount$n <-c(rep(0,6))
        cdfavg<-as.data.frame(matrix(c(0:((nsp-1)), rep(0,nsp)),nrow = nsp, ncol =2))
        colnames(cdfavg) <-c("abundance","n")
        cdfavg$n <- c(rep(0,nsp))
        compavg <-as.data.frame(matrix(c(0:(nsp-1), rep(0,nsp)),nrow = nsp, ncol =2))
        colnames(compavg) <- c("sp","n")
        
        #create distance matrices (difference, simulated, actual)
        dist <- matrix(0, ncol=40,nrow=numsim)
                colnames(dist)<-c("ks","ksunsub","pp","ppunsub","bc","bcunsub","rich","even","H","expH","wd","factor","SM","RGR","wdvar","factorvar","SMvar","RGRvar","wdcov","factorcov","SMcov","RGRcov","density","lma","lmavar","lma_unpub","lmavar_unpub","ldmc","ldmcvar","P_shade","P_shadevar","N_shade","N_shadevar","P_sun","P_sunvar","N_sun","N_sunvar","pca2","pca2var","wind")
                dist_sim <- matrix(0, ncol=40,nrow=numsim)
                colnames(dist_sim)<-c("ks","ksunsub","pp","ppunsub","bc","bcunsub","rich","even","H","expH","wd","factor","SM","RGR","wdvar","factorvar","SMvar","RGRvar","wdcov","factorcov","SMcov","RGRcov","density","lma","lmavar","lma_unpub","lmavar_unpub","ldmc","ldmcvar","P_shade","P_shadevar","N_shade","N_shadevar","P_sun","P_sunvar","N_sun","N_sunvar","pca2","pca2var","wind")
                dist_act <- matrix(0, ncol=40,nrow=numsim)
                colnames(dist_act)<-c("ks","ksunsub","pp","ppunsub","bc","bcunsub","rich","even","H","expH","wd","factor","SM","RGR","wdvar","factorvar","SMvar","RGRvar","wdcov","factorcov","SMcov","RGRcov","density","lma","lmavar","lma_unpub","lmavar_unpub","ldmc","ldmcvar","P_shade","P_shadevar","N_shade","N_shadevar","P_sun","P_sunvar","N_sun","N_sunvar","pca2","pca2var","wind")

        if(subsample_to_find_avg_distribution == FALSE){
          ######make unsumbsampled distributions######
          #run through simulations
          for(s in 1:numsim){
            testcdf <-as.data.frame(matrix(c(0:(nsp-1), rep(0,nsp)),nrow = nsp, ncol =2))
            colnames(testcdf) <- c("abundance","n")
            title = paste0(outdir,"/",method,"simulation",s-1,"allspeciescountyr",year,".txt") #read in file for that simulation
            if(file.exists(title)){
              totalcomposition <- read.delim(title,header=F)
            }else{
              cat("dne: ")
              cat(s-1)
              cat("\n")
            }
            
            #make colnum
            col_num = j + count_ngs
            testSAD1<- totalcomposition[,(col_num):(col_num+1)]
            colnames(testSAD1)<-c("sp","n")
            
            if(s == 1){
              print(paste("Doing gap column:",col_num,"\n"))
              print(paste("Note that this should read 1,3,5,7,9,11,...\n"))
            }
            
            testcdf= abundCDF(testSAD1,testcdf)
            cdfavg <- pptally(testcdf,cdfavg)
            #pp
            ppcount <- pptally(prestonplot(testSAD1,'ng110cm'),ppcount)
            #for comp
            compavg <- pptally(testSAD1,compavg)
          }
          ppavg <-averagePP(ppcount,numsim)
          compavg <- averagePP(compavg,numsim)
          cdfavgprop <-propCDF(cdfavg)
          for(b in 1:dim(ppavg)[1]){
            if(ppavg$n[b]<1){
              ppavg$n[b] = 0
            }
          }
          ppavg_unsub <- ppavg
          compavg_unsub <- compavg
          cdfavgprop_unsub <- cdfavgprop
        }else{
          ######make sumbsampled CDFs######
          #find minimum density
          simulation_minimum_density = densityactual
          for(s in startsim:endsim){
            title = paste0(outdir,"/",method,"simulation",s-1,"allspeciescountyr",year,".txt") #read in file for that simulation
            if(file.exists(title)){
              totalcomposition <- read.delim(title,header=F)
            }else{
              cat("dne: ")
              cat(s-1)
              cat("\n")
            }
            
            #make colnum
            col_num = j + count_ngs
            testSAD<- totalcomposition[,(col_num):(col_num+1)]
            colnames(testSAD)<-c("sp","n")
            density1 = density(testSAD)
            if(density1 < simulation_minimum_density){simulation_minimum_density = density1}
          }
          
          #caclulate avg cdf/pp
          for(s in startsim:endsim){
            testcdf <-as.data.frame(matrix(c(0:(nsp-1), rep(0,nsp)),nrow = nsp, ncol =2))
            colnames(testcdf) <- c("abundance","n")
            title = paste0(outdir,"/",method,"simulation",s-1,"allspeciescountyr",year,".txt") #read in file for that simulation
            if(file.exists(title)){
              totalcomposition <- read.delim(title,header=F)
            }else{
              cat("dne: ")
              cat(s-1)
              cat("\n")
            }
            
            #make colnum
            col_num = j + count_ngs
            print(col_num)
            print(j)
            print(count_ngs)
            testSAD <- totalcomposition[,(col_num):(col_num+1)]
            colnames(testSAD)<-c("sp","n")
            density1 = density(testSAD)
            
            if(s == 1){
              print(paste("Doing gap column:",col_num,"\n"))
              print(paste("Note that this should read 1,3,5,7,9,11,...\n"))
            }
          }
          
          #RAREFY
          inds <-c(0)
          compgroup <-c(inds)
          #get the difference in number of individuals
          x = density1-simulation_minimum_density
          print(x)
          print(density1)
          print(simulation_minimum_density)
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
            testSAD <- createSAD(compgroup, nsp)
          }
          
          testcdf= abundCDF(testSAD,testcdf)
          cdfavg <- pptally(testcdf,cdfavg)
          #pp
          ppcount <- pptally(prestonplot(testSAD,'ng110cm'),ppcount)
          #for comp
          compavg <- pptally(testSAD,compavg)
          print("success")
        
          #compute the average disttibutions
          ppavg <-averagePP(ppcount,numsim)
          compavg <- averagePP(compavg,numsim)
          cdfavgprop <-propCDF(cdfavg)
          for(b in 1:dim(ppavg)[1]){
            if(ppavg$n[b]<1){
              ppavg$n[b] = 0
            }
          }
        }
        
        ##### rarefy samples######     
        for(s in startsim:endsim){
          #rarefy 100 times to get accurate inference for each sample (averaging rarefied samples)
          for(y in 1:100){
            actcomp = totalcompositionactual
            title = paste0(outdir,"/",method,"simulation",s-1,"allspeciescountyr",year,".txt")
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
              #get the difference in number of individuals
              x = densitySAD-densityactual
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
                testSAD <- createSAD(compgroup, nsp)
                cat("\n")
                cat("SAD>actual")
                cat("\n")
              }
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
            #transfer both to a anew vector
            actcomp1[[s]] = actcomp
            testSAD1[[s]] = testSAD
            
            actualCDF <-as.data.frame(matrix(c(0:(nsp-1), rep(0,nsp)),nrow = nsp, ncol =2))
            colnames(actualCDF) <- c("abundance","n")
            testCDF <-as.data.frame(matrix(c(0:(nsp-1), rep(0,nsp)),nrow = nsp, ncol =2))
            colnames(testCDF) <- c("abundance","n")
            
            #calculate the cdfs (proportions)
            actcdf = abundCDF(actcomp, actualCDF)
            propactcdf = propCDF(actcdf)
            testCDF = abundCDF(testSAD, testCDF)
            proptestcdf = propCDF(testCDF)
            actcdf1[[s]] = actcdf
            testcdf1[[s]] = testCDF
            
            #manual calc of KS stat
            ksact = propactcdf$n[1]-cdfavgprop$n[1]
            for(u in 2:dim(propactcdf)[1]){
              if(abs(propactcdf$n[u]-cdfavgprop$n[u]) > ksact){
                ksact = abs(propactcdf$n[u]-cdfavgprop$n[u])
              }
            }
            kssim = proptestcdf$n[1]-cdfavgprop$n[1]
            for(u in 2:dim(proptestcdf)[1]){
              if(abs(proptestcdf$n[u]-cdfavgprop$n[u]) > kssim){
                kssim = abs(proptestcdf$n[u]-cdfavgprop$n[u])
              }
            }
            
            #manual calc of KS stat UNSUB
            ksactunsub = propactcdf$n[1]-cdfavgprop_unsub$n[1]
            for(u in 2:dim(propactcdf)[1]){
              if(abs(propactcdf$n[u]-cdfavgprop_unsub$n[u]) > ksactunsub){
                ksactunsub = abs(propactcdf$n[u]-cdfavgprop_unsub$n[u])
              }
            }
            kssimunsub = proptestcdf$n[1]-cdfavgprop_unsub$n[1]
            for(u in 2:dim(proptestcdf)[1]){
              if(abs(proptestcdf$n[u]-cdfavgprop_unsub$n[u]) > kssimunsub){
                kssimunsub = abs(proptestcdf$n[u]-cdfavgprop_unsub$n[u])
              }
            }
            dist_act[s,"ksunsub"] = ksactunsub + dist_act[s,"ksunsub"]
            dist_sim[s,"ksunsub"] = kssimunsub + dist_sim[s,"ksunsub"]
            dist[s,"ksunsub"] = ksactunsub-kssimunsub + dist[s,"ksunsub"]
            
            dist_act[s,"ks"] = ksact + dist_act[s,"ks"]
            dist_sim[s,"ks"] = kssim + dist_sim[s,"ks"]
            dist[s,"ks"] = ksact-kssim + dist[s,"ks"]
            
            #FOR PP
            actpp = prestonplot(actcomp,"ng110cm")
            testpp = prestonplot(testSAD,"ng110cm")
            dist_act[s,"pp"] = chi(actpp,ppavg) + dist_act[s,"pp"]
            dist_sim[s,"pp"] = chi(testpp,ppavg) + dist_sim[s,"pp"]
            dist[s,"pp"] = chi(actpp,ppavg) - chi(testpp,ppavg) + dist[s,"pp"]
            
            dist_act[s,"ppunsub"] = chi(actpp,ppavg_unsub) + dist_act[s,"ppunsub"]
            dist_sim[s,"ppunsub"] = chi(testpp,ppavg_unsub) + dist_sim[s,"ppunsub"]
            dist[s,"ppunsub"] = chi(actpp,ppavg_unsub) - chi(testpp,ppavg_unsub) + dist[s,"ppunsub"]
            
            #FOR BC
            dist_act[s,"bc"] = braycurtis(actcomp,compavg) + dist_act[s,"bc"]
            dist_sim[s,"bc"] = braycurtis(testSAD,compavg) + dist_sim[s,"bc"]
            dist[s,"bc"] = braycurtis(actcomp,compavg) - braycurtis(testSAD,compavg) + dist[s,"bc"]
            
            dist_act[s,"bcunsub"] = braycurtis(actcomp,compavg_unsub) + dist_act[s,"bcunsub"]
            dist_sim[s,"bcunsub"] = braycurtis(testSAD,compavg_unsub) + dist_sim[s,"bcunsub"]
            dist[s,"bcunsub"] = braycurtis(actcomp,compavg_unsub) - braycurtis(testSAD,compavg_unsub) + dist[s,"bcunsub"]
            
            
            #for richness
            dist_act[s,"rich"] = rich(actcomp) + dist_act[s,"rich"]
            dist_sim[s,"rich"] = rich(testSAD) + dist_sim[s,"rich"]
            dist[s,"rich"] = rich(actcomp) - rich(testSAD) + dist[s,"rich"]
            
            #for evenness
            dist_act[s,"even"] = E(actcomp) + dist_act[s,"even"]
            dist_sim[s,"even"] = E(testSAD) + dist_sim[s,"even"]
            dist[s,"even"] = E(actcomp) - E(testSAD) + dist[s,"even"]
            dist_act[s,"expH"] = exp(shannon(actcomp))+ dist_act[s,"expH"]
            dist_sim[s,"expH"] = exp(shannon(testSAD)) + dist_sim[s,"expH"]
            dist[s,"expH"] = exp(shannon(actcomp)) - exp(shannon(testSAD)) + dist[s,"expH"]
            #for diversity
            dist_act[s,"H"] = shannon(actcomp)+ dist_act[s,"H"]
            dist_sim[s,"H"] = shannon(testSAD) + dist_sim[s,"H"]
            dist[s,"H"] = shannon(actcomp) - shannon(testSAD) + dist[s,"H"]
            
            #for wd
            dist_act[s,"wd"] = trait(traitvals$altwd, actcomp) + dist_act[s,"wd"]
            dist_sim[s,"wd"] =trait(traitvals$altwd, testSAD) + dist_sim[s,"wd"]
            dist[s,"wd"] = trait(traitvals$altwd, actcomp) - trait(traitvals$altwd, testSAD) + dist[s,"wd"]
            
            #for factor
            dist_act[s,"factor"] = trait(traitvals$growthsurv, actcomp)  + dist_act[s,"factor"]
            dist_sim[s,"factor"] = trait(traitvals$growthsurv, testSAD) + dist_sim[s,"factor"]
            dist[s,"factor"] = trait(traitvals$growthsurv, actcomp) - trait(traitvals$growthsurv, testSAD) + dist[s,"factor"]
            #("ks","pp","bc","rich","even","H","wd","factor","SM","RGR","wdvar","factorvar","SMvar","RGRvar","wdcov","factorcov","SMcov","RGRcov","density")
            
            dist[s,"SM"] = trait(traitvals$sm, actcomp) - trait(traitvals$sm, testSAD) + dist[s,"SM"]
            #for factor
            dist[s,"RGR"] = trait(traitvals$RGR, actcomp) - trait(traitvals$RGR, testSAD) + dist[s,"RGR"]
            
            dist[s,"wdvar"] = traitvar(traitvals$altwd, actcomp) - traitvar(traitvals$altwd, testSAD) + dist[s,"wdvar"]
            #for factor
            dist[s,"factorvar"] = traitvar(traitvals$growthsurv, actcomp) - traitvar(traitvals$growthsurv, testSAD) + dist[s,"factorvar"]
            
            dist[s,"SMvar"] = traitvar(traitvals$sm, actcomp) - traitvar(traitvals$sm, testSAD) + dist[s,"SMvar"]
            #for factor
            dist[s,"RGRvar"] = traitvar(traitvals$RGR, actcomp) - traitvar(traitvals$RGR, testSAD) + dist[s,"RGRvar"]
            
            dist[s,"wdcov"] = traitcov(traitvals$altwd, actcomp) - traitcov(traitvals$altwd, testSAD) + dist[s,"wdcov"]
            #for factor
            dist[s,"factorcov"] = traitcov(traitvals$growthsurv, actcomp) - traitcov(traitvals$growthsurv, testSAD) + dist[s,"factorcov"]
            
            dist[s,"SMcov"] = traitcov(traitvals$sm, actcomp) - traitcov(traitvals$sm, testSAD) + dist[s,"SMcov"]
            #for factor
            dist[s,"RGRcov"] = traitcov(traitvals$RGR, actcomp) - traitcov(traitvals$RGR, testSAD) + dist[s,"RGRcov"]
            
            dist_act[s,"RGR"] = trait(traitvals$RGR, actcomp) + dist_act[s,"RGR"]
            dist_sim[s,"RGR"] =trait(traitvals$RGR, testSAD) + dist_sim[s,"RGR"]
            #for factor
            dist_act[s,"SM"] = trait(traitvals$sm, actcomp)  + dist_act[s,"SM"]
            dist_sim[s,"SM"] = trait(traitvals$sm, testSAD) + dist_sim[s,"SM"]
            
            dist_act[s,"wdcov"] = traitcov(traitvals$altwd, actcomp) + dist_act[s,"wdcov"]
            dist_sim[s,"wdcov"] =traitcov(traitvals$altwd, testSAD) + dist_sim[s,"wdcov"]
            #for factor
            dist_act[s,"factorcov"] = traitcov(traitvals$growthsurv, actcomp)  + dist_act[s,"factorcov"]
            dist_sim[s,"factorcov"] = traitcov(traitvals$growthsurv, testSAD) + dist_sim[s,"factorcov"]
            
            dist_act[s,"RGRcov"] = traitcov(traitvals$RGR, actcomp) + dist_act[s,"RGRcov"]
            dist_sim[s,"RGRcov"] =traitcov(traitvals$RGR, testSAD) + dist_sim[s,"RGRcov"]
            #for factor
            dist_act[s,"SMcov"] = traitcov(traitvals$sm, actcomp)  + dist_act[s,"SMcov"]
            dist_sim[s,"SMcov"] = traitcov(traitvals$sm, testSAD) + dist_sim[s,"SMcov"]
            
            dist_act[s,"wdvar"] = traitvar(traitvals$altwd, actcomp) + dist_act[s,"wdvar"]
            dist_sim[s,"wdvar"] =traitvar(traitvals$altwd, testSAD) + dist_sim[s,"wdvar"]
            #for factor
            dist_act[s,"factorvar"] = traitvar(traitvals$growthsurv, actcomp)  + dist_act[s,"factorvar"]
            dist_sim[s,"factorvar"] = traitvar(traitvals$growthsurv, testSAD) + dist_sim[s,"factorvar"]
            
            dist_act[s,"RGRvar"] = traitvar(traitvals$RGR, actcomp) + dist_act[s,"RGRvar"]
            dist_sim[s,"RGRvar"] =traitvar(traitvals$RGR, testSAD) + dist_sim[s,"RGRvar"]
            #for factor
            dist_act[s,"SMvar"] = traitvar(traitvals$sm, actcomp)  + dist_act[s,"SMvar"]
            dist_sim[s,"SMvar"] = traitvar(traitvals$sm, testSAD) + dist_sim[s,"SMvar"]
            
            dist_act[s,"density"] = densityactual  + dist_act[s,"density"]
            dist_sim[s,"density"] =  densitySAD + dist_sim[s,"density"]
            dist[s,"density"] = densityactual - densitySAD + dist[s,"density"]
            
            #dist for var
            dist[s,"lmavar"] = traitvar(traitvals$lma, actcomp) - traitvar(traitvals$lma, testSAD) + dist[s,"lmavar"]
            dist[s,"ldmcvar"] = traitvar(traitvals$ldmc, actcomp) - traitvar(traitvals$ldmc, testSAD) + dist[s,"ldmcvar"]
            dist[s,"N_shadevar"] = traitvar(traitvals$N_shade, actcomp) - traitvar(traitvals$N_shade, testSAD) + dist[s,"N_shadevar"]
            dist[s,"P_shadevar"] = traitvar(traitvals$P_shade, actcomp) - traitvar(traitvals$P_shade, testSAD) + dist[s,"P_shadevar"]
            dist[s,"N_sunvar"] = traitvar(traitvals$N_sun, actcomp) - traitvar(traitvals$N_sun, testSAD) + dist[s,"N_sunvar"]
            dist[s,"P_sunvar"] = traitvar(traitvals$P_sun, actcomp) - traitvar(traitvals$P_sun, testSAD) + dist[s,"P_sunvar"]
            dist[s,"pca2var"] = traitvar(traitvals$pca2, actcomp) - traitvar(traitvals$pca2, testSAD) + dist[s,"pca2var"]
            
            dist_act[s,"lmavar"] = traitvar(traitvals$lma, actcomp)  + dist_act[s,"lmavar"]
            dist_sim[s,"lmavar"] = traitvar(traitvals$lma, testSAD) + dist_sim[s,"lmavar"]
            
            dist_act[s,"ldmcvar"] = traitvar(traitvals$ldmc, actcomp) + dist_act[s,"ldmcvar"]
            dist_sim[s,"ldmcvar"] =traitvar(traitvals$ldmc, testSAD) + dist_sim[s,"ldmcvar"]
            
            dist_act[s,"N_shadevar"] = traitvar(traitvals$N_shade, actcomp)  + dist_act[s,"N_shadevar"]
            dist_sim[s,"N_shadevar"] = traitvar(traitvals$N_shade, testSAD) + dist_sim[s,"N_shadevar"]
            
            dist_act[s,"P_shadevar"] = traitvar(traitvals$P_shade, actcomp)  + dist_act[s,"P_shadevar"]
            dist_sim[s,"P_shadevar"] = traitvar(traitvals$P_shade, testSAD) + dist_sim[s,"P_shadevar"]
            
            dist_act[s,"N_sunvar"] = traitvar(traitvals$N_sun, actcomp)  + dist_act[s,"N_sunvar"]
            dist_sim[s,"N_sunvar"] = traitvar(traitvals$N_sun, testSAD) + dist_sim[s,"N_sunvar"]
            
            dist_act[s,"P_sunvar"] = traitvar(traitvals$P_sun, actcomp)  + dist_act[s,"P_sunvar"]
            dist_sim[s,"P_sunvar"] = traitvar(traitvals$P_sun, testSAD) + dist_sim[s,"P_sunvar"]
            
            dist_act[s,"pca2var"] = traitvar(traitvals$pca2, actcomp)  + dist_act[s,"pca2var"]
            dist_sim[s,"pca2var"] = traitvar(traitvals$pca2, testSAD) + dist_sim[s,"pca2var"]
            
            dist_act[s,"lma"] = trait(traitvals$lma, actcomp) + dist_act[s,"lma"]
            dist_sim[s,"lma"] =trait(traitvals$lma, testSAD) + dist_sim[s,"lma"]
            dist[s,"lma"] = trait(traitvals$lma, actcomp) - trait(traitvals$lma, testSAD) + dist[s,"lma"]
            
            dist_act[s,"ldmc"] = trait(traitvals$ldmc, actcomp)  + dist_act[s,"ldmc"]
            dist_sim[s,"ldmc"] = trait(traitvals$ldmc, testSAD) + dist_sim[s,"ldmc"]
            dist[s,"ldmc"] = trait(traitvals$ldmc, actcomp) - trait(traitvals$ldmc, testSAD) + dist[s,"ldmc"]
            
            dist_act[s,"N_sun"] = trait(traitvals$N_sun, actcomp) + dist_act[s,"N_sun"]
            dist_sim[s,"N_sun"] =trait(traitvals$N_sun, testSAD) + dist_sim[s,"N_sun"]
            dist[s,"N_sun"] = trait(traitvals$N_sun, actcomp) - trait(traitvals$N_sun, testSAD) + dist[s,"N_sun"]
            
            dist_act[s,"N_shade"] = trait(traitvals$N_shade, actcomp) + dist_act[s,"N_shade"]
            dist_sim[s,"N_shade"] =trait(traitvals$N_shade, testSAD) + dist_sim[s,"N_shade"]
            dist[s,"N_shade"] = trait(traitvals$N_shade, actcomp) - trait(traitvals$N_shade, testSAD) + dist[s,"N_shade"]
            
            dist_act[s,"P_sun"] = trait(traitvals$P_sun, actcomp)  + dist_act[s,"P_sun"]
            dist_sim[s,"P_sun"] = trait(traitvals$P_sun, testSAD) + dist_sim[s,"P_sun"]
            dist[s,"P_sun"] = trait(traitvals$P_sun, actcomp) - trait(traitvals$P_sun, testSAD) + dist[s,"P_sun"]
            
            dist_act[s,"P_shade"] = trait(traitvals$P_shade, actcomp)  + dist_act[s,"P_shade"]
            dist_sim[s,"P_shade"] = trait(traitvals$P_shade, testSAD) + dist_sim[s,"P_shade"]
            dist[s,"P_shade"] = trait(traitvals$P_shade, actcomp) - trait(traitvals$P_shade, testSAD) + dist[s,"P_shade"]
            
            dist_act[s,"pca2"] = trait(traitvals$pca2, actcomp) + dist_act[s,"pca2"]
            dist_sim[s,"pca2"] =trait(traitvals$pca2, testSAD) + dist_sim[s,"pca2"]
            dist[s,"pca2"] = trait(traitvals$pca2, actcomp) - trait(traitvals$pca2, testSAD) + dist[s,"pca2"]
            
            dist_act[s,"wind"] = trait(traitvals$wind, actcomp) + dist_act[s,"wind"]
            dist_sim[s,"wind"] =trait(traitvals$wind, testSAD) + dist_sim[s,"wind"]
            dist[s,"wind"] = trait(traitvals$wind, actcomp) - trait(traitvals$wind, testSAD) + dist[s,"wind"]
            
            
            title = paste0("results/trash",method,year,simtype,ng_num,"distresultsactual.csv")
            write.csv(dist_act,title,row.names = F)
            title = paste0("results/trash",method,year,simtype,ng_num,"distresultssimulated.csv")
            write.csv(dist_sim,title,row.names = F)
            title = paste0("results/trash",method,year,simtype,ng_num,"distresults.csv")
            write.csv(dist,title,row.names = F)
            
          }
          cat("\n")
          cat("done with simulation")
          cat("\n")
        } 
        #average across rarefications
        dist = dist/100
        dist_act = dist_act/100
        dist_sim = dist_sim/100
        colnames(dist)<-c("KS_CDF_test","KS_CDF_test_unsamplesavg","Chi_Sq_Preston_Plot","Chi_Sq_Preston_Plot_unsamplesavg","Bray_Curtis_Dissimilarity","Bray_Curtis_Dissimilarity_unsamplesavg","Richness","Evenness","Shannons_H_Diversity","Effective_Species_expH","Wood_Density","Factor_Score","Seed_Mass","RGR","WoddDensity_Variance","Factor_variance","SMvariance","RGRvariance","WDcovar","Factorcovar","SMcovar","RGRcovar","density","LMA_published_DISC","LMA_published_DISC_variance","LMA_UNpublished_DISC_variance","LMA_UNpublished_DISC_variance","LDMC_unpublished","LDMC_unpublished_variance","Pshade","Pshadevar","Nshade","Nshadevar","Phosphorous_sun","Phosphorous_sun_var","Nitrogen_sun","Nitrogen_sun_variance","pca2","pca2var","WindDispersal")
        colnames(dist_sim)<-c("KS_CDF_test","KS_CDF_test_unsamplesavg","Chi_Sq_Preston_Plot","Chi_Sq_Preston_Plot_unsamplesavg","Bray_Curtis_Dissimilarity","Bray_Curtis_Dissimilarity_unsamplesavg","Richness","Evenness","Shannons_H_Diversity","Effective_Species_expH","Wood_Density","Factor_Score","Seed_Mass","RGR","WoddDensity_Variance","Factor_variance","SMvariance","RGRvariance","WDcovar","Factorcovar","SMcovar","RGRcovar","density","LMA_published_DISC","LMA_published_DISC_variance","LMA_UNpublished_DISC_variance","LMA_UNpublished_DISC_variance","LDMC_unpublished","LDMC_unpublished_variance","Pshade","Pshadevar","Nshade","Nshadevar","Phosphorous_sun","Phosphorous_sun_var","Nitrogen_sun","Nitrogen_sun_variance","pca2","pca2var","WindDispersal")
        colnames(dist_act)<-c("KS_CDF_test","KS_CDF_test_unsamplesavg","Chi_Sq_Preston_Plot","Chi_Sq_Preston_Plot_unsamplesavg","Bray_Curtis_Dissimilarity","Bray_Curtis_Dissimilarity_unsamplesavg","Richness","Evenness","Shannons_H_Diversity","Effective_Species_expH","Wood_Density","Factor_Score","Seed_Mass","RGR","WoddDensity_Variance","Factor_variance","SMvariance","RGRvariance","WDcovar","Factorcovar","SMcovar","RGRcovar","density","LMA_published_DISC","LMA_published_DISC_variance","LMA_UNpublished_DISC_variance","LMA_UNpublished_DISC_variance","LDMC_unpublished","LDMC_unpublished_variance","Pshade","Pshadevar","Nshade","Nshadevar","Phosphorous_sun","Phosphorous_sun_var","Nitrogen_sun","Nitrogen_sun_variance","pca2","pca2var","WindDispersal")
        
        
        title = paste0("results/all",method,year,simtype,ng_num,"distresultsactual.csv")
        write.csv(dist_act,title,row.names = F)
        title = paste0("results/all",method,year,simtype,ng_num,"distresultssimulated.csv")
        write.csv(dist_sim,title,row.names = F)
        title = paste0("results/all",method,year,simtype,ng_num,"distresults.csv")
        write.csv(dist,title,row.names = F)
      }
}

