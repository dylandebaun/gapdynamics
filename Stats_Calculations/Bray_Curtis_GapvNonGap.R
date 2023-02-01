method = "final"
numsim = 500
type="all"
simtype=""
recruittype=""
args <- commandArgs(trailingOnly = TRUE)
year=as.numeric(args[1])
gap_number=as.numeric(args[2])
out=as.character(args[3])
method=as.character(args[4])

nsp = 328
numsim=500
library(stats)
source("sourcefunctions.R")
ng_num = ""
title = paste0(out, "/",method,simtype,ng_num,recruittype,"expectedin",type, "at",year,"yrs.txt")
totalcompositionactual <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")
densityactual = density(totalcompositionactual) #figure out density of the actual gap for subsampling

#create distance matrices
dist <- matrix(0, ncol=3,nrow=numsim)
colnames(dist)<-c("bc_act","bc_sim","bc_diff")

for(j in gap_number:gap_number){
	simtype="nongap"
	ng_num = j-1
  title = paste0("final_actual/","final",simtype,j-1,recruittype,"expectedin",type, "at",year,"yrs.txt")
  totalcompositionactual_ng <- read.delim(title, header =  F, col.names =  c("sp", "n"),sep="")
  densityactual_ng = density(totalcompositionactual_ng) #figure out density of the actual gap for subsampling
  
##### rarefy samples to match simulation######   
  #nongaps
  for(s in 1:numsim){
    #rarefy 100 times to get accurate inference for each sample (averaging rarefied samples)
    for(y in 1:100){
      #nongap
      actcomp_ng = totalcompositionactual_ng
      title = paste0("final10cm/",method,"simulation",s-1,type,"speciescountyr",year,".txt")
      if(file.exists(title)){
        totalcomposition_ng <- read.delim(title,header=F)
      }else{
        cat("dne: ")
        cat(s-1)
        cat("\n")
      }
      #go to the column we want in the data space
      print(j)
      col_num=0
      print(col_num)
      if(j == 1){
	      col_num = 3
	      testSAD_ng<- totalcomposition_ng[,3:4]
      }else if(j == 2){
              col_num = 5
	      testSAD_ng<- totalcomposition_ng[,5:6]
      }else if(j == 3){
              testSAD_ng<- totalcomposition_ng[,7:8]
      }else if(j == 4){
              col_num = 9
	      testSAD_ng<- totalcomposition_ng[,9:10]
      }else if(j == 5){
              col_num = 11
	      testSAD_ng<- totalcomposition_ng[,11:12]
      }else if(j == 6){
              col_num = 13
	      testSAD_ng<- totalcomposition_ng[,13:14]
      }else if(j == 7){
              col_num = 15
	      testSAD_ng<- totalcomposition_ng[,15:16]
      }else if(j == 8){
              col_num = 17
	      testSAD_ng<- totalcomposition_ng[,17:18]
      }else if(j == 9){
              col_num = 19
	      testSAD_ng<- totalcomposition_ng[,19:20]
      }
      print(col_num)
      #testSAD_ng<- totalcomposition_ng[,col_num:(col_num+1)]
      #cat("fail1")
      colnames(testSAD_ng)<-c("sp","n")
      
      #get the density for rarefying
      densitySAD_ng <- density(testSAD_ng)
      print(densitySAD_ng)
      if(densitySAD_ng==0){
        cat("empty: ")
        cat(s-1)
        cat("\n")
      }
      
      #RAREFY
      if(year != 0){
        inds <-c(0)
        compgroup <-c(inds)
        if(isTRUE(densitySAD_ng > densityactual_ng)){
          #get the difference in number of individuals
          x = densitySAD_ng-densityactual_ng
          if(x != 0){
            #expand the composition vector into a group of individuals
            for(d in 1:dim(testSAD_ng)[1]){
              if(testSAD_ng$n[d] != 0){
                for(z in 1:testSAD_ng$n[d]){
                  compgroup <- rbind(compgroup, c(testSAD_ng$sp[d]))
                }
              }
            }
            compgroup <- compgroup[-c(1),]
            #randomly remove the individuals
            rm =  sample(1:length(compgroup),x,replace = F)
            compgroup <- compgroup[-c(rm)]
            #rebuild the SAD
            testSAD_ng <- createSAD(compgroup, nsp)
	    print(testSAD_ng)
            cat("\n")
            cat("SAD>actual")
            cat("\n")
          }
        }else if(isTRUE(densitySAD_ng < densityactual_ng)){
          x = densityactual_ng-densitySAD_ng
          if(x != 0){
            #expand the composition vector into a group of individuals
            for(d in 1:dim(actcomp_ng)[1]){
              if(actcomp_ng$n[d] != 0){
                for(z in 1:actcomp_ng$n[d]){
                  compgroup <- rbind(compgroup, c(actcomp_ng$sp[d]))
                }
              }
            }
            compgroup <- compgroup[-c(1),]
            rm =  sample(1:length(compgroup),x,replace = F)
            compgroup <- compgroup[-c(rm)]
            actcomp_ng <- createSAD(compgroup, nsp)
          print("actual > SAD")
	  }
        }
      }
      #gap
      actcomp = totalcompositionactual
      title = paste0("final10cm/",method,"simulation",s-1,type,"speciescountyr",year,".txt")
      if(file.exists(title)){
        totalcomposition <- read.delim(title,header=F)
      }else{
        cat("dne: ")
        cat(s-1)
        cat("\n")
      }
      #go to the column we want in the data space
     # col_num =  1
      testSAD<- totalcomposition[,1:(1+1)]
      #cat("fail1")
      colnames(testSAD)<-c("sp","n")
      
      #get the density for rarefying
      densitySAD <- density(testSAD)
      if(densitySAD==0){
        cat("empty: ")
        cat(s-1)
        cat("\n")
      }
      
      #RAREFY
      if(year != 0){
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
      }
      
      
      #FOR BC
      dist[s,"bc_act"] = braycurtis(actcomp,actcomp_ng) +  dist[s,"bc_act"]
      dist[s,"bc_sim"] = braycurtis(testSAD,testSAD_ng) +dist[s,"bc_sim"]
      #dist[s,"bc_diff"] = dist[s,"bc_act"] - dist[s,"bc_sim"] + dist[s,"bc_diff"]
      dist[s,"bc_diff"] =  braycurtis(actcomp,actcomp_ng) - braycurtis(testSAD,testSAD_ng) + dist[s,"bc_diff"]
      title = paste0("results/trash",type,method,year,simtype,recruittype,ng_num,"distresults_gngbc.csv")
      write.csv(dist,title,row.names = F)
    }
    cat("\n")
    cat("done with simulation")
    cat("\n")
  } 
#average across rarefications
dist = dist/100

title = paste0("results/",type,method,year,simtype,recruittype,ng_num,"distresults_gngbc.csv")
write.csv(dist,title,row.names = F)
}

