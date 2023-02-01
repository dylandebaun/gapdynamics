#CHANGED COV TO OUTPUT IQR AND VAR TO OUTPUT SD
makerankab <- function(SAD){
  sadcount <- matrix(c(1:328, rep(0,328)),nrow=328, ncol = 2)
maxtomin <- SAD[order(SAD[,2],decreasing=T),]
total <- sum(SAD[,2])
sadcount[,2] = maxtomin[,2]
sadcount[,2] = sadcount[,2]/total
sadcount[,1] = maxtomin[,1]
#sadcount[,2] = log10(sadcount[,2])
#sadcount[is.infinite(sadcount[,2]),2] = 0
	return(sadcount)
}

orderSAD <- function(SAD){
  for (k in 1:dim(SAD)[1]){
    for (l in 1:dim(SAD)[1]){
      if (isTRUE(SAD$n[l] < SAD$n[k])) {
        tmp = SAD$n[k];
        tmp2 = SAD$sp[k];
        SAD$n[k] = SAD$n[l];
        SAD$n[l] = tmp;
        SAD$sp[k] = SAD$sp[l];
        SAD$sp[l] = tmp2;
      }
    }
  }
  return(SAD)
}

remove0sSAD <- function(SAD){
  for(p in 1:100){
    for(j in 1:dim(SAD)[1]){
      if(isTRUE(SAD[j, 3]==0)){
        SAD<-SAD[-c(j),]
      }
    }
  }
  return(SAD)
}

prestonplotx <- function(SAD){
  abundancecategories <-c("1","2 to 3","4 to 7","8 to 15","16-31","32-63",
                          "64-127","128-255","256-511", ">511")
  n <- c(rep(0,10))
  prestonplot <- data.frame(abundancecategories,n)
  for(i in 1:dim(SAD)[1]){
    if(isTRUE(SAD$n[i] <= 1) && isTRUE(SAD$n[i] > 0)){
      prestonplot$n[1] = prestonplot$n[1] +1
    } else if(isTRUE(SAD$n[i] >= 2) && isTRUE(SAD$n[i] <= 3)){
      prestonplot$n[2] = prestonplot$n[2] +1
    }else if(isTRUE(SAD$n[i] >= 4) && isTRUE(SAD$n[i] <= 7)){
      prestonplot$n[3] = prestonplot$n[3] +1
    }else if(isTRUE(SAD$n[i] >= 8) && isTRUE(SAD$n[i] <= 15)){
      prestonplot$n[4] = prestonplot$n[4] +1
    }else if(isTRUE(SAD$n[i] >= 16) && isTRUE(SAD$n[i] <= 31)){
      prestonplot$n[5] = prestonplot$n[5] +1
    }else if(isTRUE(SAD$n[i] >= 32) && isTRUE(SAD$n[i] <= 63)){
      prestonplot$n[6] = prestonplot$n[6] +1
    }else if(isTRUE(SAD$n[i] >= 64) && isTRUE(SAD$n[i] <= 127)){
      prestonplot$n[7] = prestonplot$n[7] +1
    }else if(isTRUE(SAD$n[i] >= 128) && isTRUE(SAD$n[i] <= 255)){
      prestonplot$n[8] = prestonplot$n[8] +1
    }else if(isTRUE(SAD$n[i] >= 256) && isTRUE(SAD$n[i] <= 511)){
      prestonplot$n[9] = prestonplot$n[9] +1
    }else if(isTRUE(SAD$n[i] >= 512)){
      prestonplot$n[10] = prestonplot$n[10] +1
    }
  }
  colnames(prestonplot)<- c("category", "n")
  return(prestonplot)
}
prestonplot <- function(SAD,method){
  if(method == "10cm85"|| method == "ng110cm"|| method == "ng210cm"){
  abundancecategories <-c("1","2 to 3","4 to 7","8 to 15","16-31","32+")
  n <- c(rep(0,6))
  prestonplot <- data.frame(abundancecategories,n)
  for(i in 1:dim(SAD)[1]){
    if(isTRUE(SAD$n[i] <= 1) && isTRUE(SAD$n[i] > 0)){
      prestonplot$n[1] = prestonplot$n[1] +1
    } else if(isTRUE(SAD$n[i] >= 2) && isTRUE(SAD$n[i] <= 3)){
      prestonplot$n[2] = prestonplot$n[2] +1
    }else if(isTRUE(SAD$n[i] >= 4) && isTRUE(SAD$n[i] <= 7)){
      prestonplot$n[3] = prestonplot$n[3] +1
    }else if(isTRUE(SAD$n[i] >= 8) && isTRUE(SAD$n[i] <= 15)){
      prestonplot$n[4] = prestonplot$n[4] +1
    }else if(isTRUE(SAD$n[i] >= 16) && isTRUE(SAD$n[i] <= 31)){
      prestonplot$n[5] = prestonplot$n[5] +1
    }else if(isTRUE(SAD$n[i] >= 32)){
      prestonplot$n[6] = prestonplot$n[6] +1
    }
  }
  }else if(method == "rf85"|| method == "ng1rf"|| method == "ng2rf"){
    abundancecategories <-c("1","2 to 3","4 to 7","8 to 15","16-31","32-63",
                            "64-127","128+")
    n <- c(rep(0,8))
    prestonplot <- data.frame(abundancecategories,n)
    for(i in 1:dim(SAD)[1]){
      if(isTRUE(SAD$n[i] <= 1) && isTRUE(SAD$n[i] > 0)){
        prestonplot$n[1] = prestonplot$n[1] +1
      } else if(isTRUE(SAD$n[i] >= 2) && isTRUE(SAD$n[i] <= 3)){
        prestonplot$n[2] = prestonplot$n[2] +1
      }else if(isTRUE(SAD$n[i] >= 4) && isTRUE(SAD$n[i] <= 7)){
        prestonplot$n[3] = prestonplot$n[3] +1
      }else if(isTRUE(SAD$n[i] >= 8) && isTRUE(SAD$n[i] <= 15)){
        prestonplot$n[4] = prestonplot$n[4] +1
      }else if(isTRUE(SAD$n[i] >= 16) && isTRUE(SAD$n[i] <= 31)){
        prestonplot$n[5] = prestonplot$n[5] +1
      }else if(isTRUE(SAD$n[i] >= 32) && isTRUE(SAD$n[i] <= 63)){
        prestonplot$n[6] = prestonplot$n[6] +1
      }else if(isTRUE(SAD$n[i] >= 64) && isTRUE(SAD$n[i] <= 127)){
        prestonplot$n[7] = prestonplot$n[7] +1
      }else if(isTRUE(SAD$n[i] >= 128)){
        prestonplot$n[8] = prestonplot$n[8] +1
      }
    }
  }
  colnames(prestonplot)<- c("category", "n")
  return(prestonplot)
}

sortgaps <-function(gaplocations){
  for (k in 1:dim(gaplocations)[1]){
    for (l in 1:dim(gaplocations)[1]){
      if (isTRUE((gaplocations$size[l]) < (gaplocations$size[k]))) {
        tmp = gaplocations$gap[k]
        tmp2 = gaplocations$size[k]
        tmp3 = gaplocations$gapnum[k]
        gaplocations$gap[k] = gaplocations$gap[l]
        gaplocations$gap[l] = tmp
        gaplocations$size[k] = gaplocations$size[l]
        gaplocations$size[l] = tmp2
        gaplocations$gapnum[k] = gaplocations$gapnum[l]
        gaplocations$gapnum[l] = tmp3
      }
    }
  }
  return(gaplocations)
 }
# sortgaps <-function(gaplocations){
#   for (k in 1:dim(gaplocations)[1]){
#     for (l in 1:dim(gaplocations)[1]){
#       if (isTRUE((gaplocations$gapsizex[l]*gaplocations$gapsizey[l]) < (gaplocations$gapsizex[k]*gaplocations$gapsizey[k]))) {
#         tmp = gaplocations$x[k]
#         tmp2 = gaplocations$y[k]
#         tmp3=gaplocations$gapsizex[k]
#         tmp4 = gaplocations$gapsizey[k]
#         tmp5 = gaplocations$gapnum[k]
#         gaplocations$x[k] = gaplocations$x[l]
#         gaplocations$x[l] = tmp
#         gaplocations$y[k] = gaplocations$y[l]
#         gaplocations$y[l] = tmp2
#         gaplocations$gapsizex[k] = gaplocations$gapsizex[l]
#         gaplocations$gapsizex[l] = tmp3
#         gaplocations$gapsizey[k] =  gaplocations$gapsizey[l]
#         gaplocations$gapsizey[l] = tmp4
#         gaplocations$gapnum[k] =  gaplocations$gapnum[l]
#         gaplocations$gapnum[l] = tmp5
#         
#       }
#     }
#   }
#   return(gaplocations)
# }

orderbyotherSAD <- function(averageSAD,actualSAD){
    for(i in 1:dim(actualSAD)[1]){
      for(j in 1:dim(averageSAD)[1]){
        if(isTRUE(actualSAD$sp[i] == averageSAD$sp[j])){
          temp = averageSAD$sp[i]
          averageSAD$sp[i] = averageSAD$sp[j]
          averageSAD$sp[j] = temp
          temp1 = averageSAD$n[i]
          averageSAD$n[i] = averageSAD$n[j]
          averageSAD$n[j] = temp1
        }
      }
    }
  return(averageSAD)
}

pptally <- function(ppsim, ppcount){
  for(i in 1:dim(ppsim)[1]){
    ppcount$n[i] = ppsim$n[i] + ppcount$n[i]
  }
  return(ppcount)
}

averagePP <- function(ppcount,numsim){
  for(i in 1:dim(ppcount)[1]){
    ppcount$n[i] = ppcount$n[i]/numsim
  }
  return(ppcount)
}

propPP <- function(avgpp){
  sum=0
  for(i in 1:dim(avgpp)[1]){
    sum= avgpp$n[i] + sum
  }
  for(i in 1:dim(avgpp)[1]){
    avgpp$n[i] = avgpp$n[i]/sum
  }
  return(avgpp)
}

orderbysp <- function(SAD){
  for (k in 1:dim(SAD)[1]){
    for (l in 1:dim(SAD)[1]){
      if (isTRUE(SAD$sp[l] > SAD$sp[k])) {
        tmp = SAD$n[k];
        tmp2 = SAD$sp[k];
        SAD$n[k] = SAD$n[l];
        SAD$n[l] = tmp;
        SAD$sp[k] = SAD$sp[l];
        SAD$sp[l] = tmp2;
      }
    }
  }
  return(SAD)
}

abundCDF <- function(SAD, sadcount){
  max = max(SAD$n, na.rm = TRUE)
  max = max+1
  if(isFALSE(max == 0))
  for(i in 2:(dim(sadcount)[1])){
      count = 0
      if(isTRUE(i <= max)){
        for(j in 1:dim(SAD)[1]){
          if(isTRUE(SAD$n[j] == sadcount$abundance[i])){
              count = count + 1
          }
        }
        sadcount$n[i] = sadcount$n[i-1] + count 
      }else{
      sadcount$n[i] = sadcount$n[max]
      }
   # if(sadcount$abundance[i] >= 201){
     #   cat("above 200 abund!\n")
   # }
  }
  return(sadcount)
}

propCDF <- function(cdf){
  total = cdf$n[dim(cdf)[1]]
  if(isFALSE(total == 0)){
  for(i in 1: dim(cdf)[1]){
    cdf$n[i] = cdf$n[i]/total
  }
  }
  return(cdf)
}


shannon <- function(SAD){
  H = 0
  sumofsp = 0
  for(k in 1:dim(SAD)[1]){
    sumofsp = sumofsp + SAD$n[k]
  }
  for(k in 1:dim(SAD)[1]){
    if(isFALSE(SAD$n[k] == 0)){
      H = H + ((SAD$n[k]/sumofsp)*log(SAD$n[k]/sumofsp))
    }
  }
  return(-H)
}

E <- function(SAD){
  H = 0
  sumofsp = 0
  for(k in 1:dim(SAD)[1]){
    sumofsp = sumofsp + SAD$n[k]
  }
  for(k in 1:dim(SAD)[1]){
    if(isFALSE(SAD$n[k] == 0)){
      H = H + ((SAD$n[k]/sumofsp)*log(SAD$n[k]/sumofsp))
    }
  }
  evenness = -H/(log(sumofsp))
  return(evenness)
}
  

msetest <- function(SADx, SADy){
  mse = 0
  for(k in 1:dim(SADy)[1]){
    if(isFALSE(SADx$n[k] == 0)){
      mse = (((SADy$n[k]- SADx$n[k])^2)/SADx$n[k]) + mse
    }
  }
  #mse = mse/dim(SADy)[1]
  return(mse)
}

chi <-function(ppobs, ppexp){
  chisq = 0
  for(j in 1:dim(ppobs)[1]){
    if(isFALSE(ppexp$n[j] == 0)){
      chisq = ((ppobs$n[j]-ppexp$n[j])*(ppobs$n[j]-ppexp$n[j])/ppexp$n[j]) + chisq
      #sum of (observed-expected)^2/expected
    }
  }
  return(chisq)
}
# chi <-function(ppobs, ppexp){
#   chisq = 0
#   for(j in 1:dim(ppobs)[1]){
#     obs = ppobs$n[j] + 1
#     exp = ppexp$n[j] + 1
#       chisq = ((obs-exp)*(obs-exp)/exp) + chisq
#   }
#   return(chisq)
# }

trait<- function(traitvals, testSAD){
  avg = 0
  count = 0
  for(i in 1:dim(testSAD)[1]){
    if(isFALSE(is.na(traitvals[i])) && isTRUE(is.numeric(traitvals[i]))){
    avg = (testSAD$n[i] * traitvals[i]) + avg
    count = testSAD$n[i] +count
    }
  }
  avg = avg/count
  var = 0
  for(i in 1:dim(testSAD)[1]){
        if(isFALSE(is.na(traitvals[i])) && isTRUE(is.numeric(traitvals[i]))){
      var = var + (((testSAD$n[i] * traitvals[i]) - avg)*((testSAD$n[i] * traitvals[i]) - avg))
    }
  }
  var = var/(count-1)
  return(avg)
}

traitvar<- function(traitvals, testSAD){
  avg = 0
  count = 0
  for(i in 1:dim(testSAD)[1]){
    if(isFALSE(is.na(traitvals[i])) && isTRUE(is.numeric(traitvals[i]))){
    avg = (testSAD$n[i] * traitvals[i]) + avg
    count = testSAD$n[i] +count
    }
  }
  avg = avg/count
    var = 0
  for(i in 1:dim(testSAD)[1]){
        if(isFALSE(is.na(traitvals[i])) && isTRUE(is.numeric(traitvals[i]))){
        for(p in 1:testSAD$n[i]){
            var = var + (((traitvals[i]) - avg)*((traitvals[i]) - avg))
          }
        }
      }
  var = var/(count-1)
  cov = sqrt(var)/avg

  return((var))
}
traitcov<- function(traitvals, testSAD){
    avg = 0
    count = 0
    traits <- c(0)
    for(i in 1:dim(testSAD)[1]){
 if(isFALSE(is.na(traitvals[i])) && isTRUE(is.numeric(traitvals[i]))){
            avg = (testSAD$n[i] * traitvals[i]) + avg
            count = testSAD$n[i] +count
	    traits <- c(traits,rep(traitvals[i],testSAD$n[i]))
        }
    }
    avg = avg/count
    var = 0
      var = 0
  for(i in 1:dim(testSAD)[1]){
     if(isFALSE(is.na(traitvals[i])) && isTRUE(is.numeric(traitvals[i]))){
        for(p in 1:testSAD$n[i]){
            var = var + (((traitvals[i]) - avg)*((traitvals[i]) - avg))
          }
        }
      }
  var = var/(count-1)
    cov = sqrt(var)/avg
	traits <-traits[-c(1)]

    return(cov)
     #return(IQR(traits))	
}

rich <- function(testSAD){
  richness = 0
  for(i in 1:dim(testSAD)[1]){
    if(isFALSE(testSAD$n[i] < 0.5)){
      richness = richness + 1
    }
  }
  return(richness)
}

braycurtis <- function(SAD1,SAD2){
  numerator = 0
  numindiv1 = 0
  numindiv2 = 0
  for(i in 1:dim(SAD1)[1]){
    numindiv1 = numindiv1 + SAD1$n[i]
    numindiv2 = numindiv2 + SAD2$n[i]
    if(isTRUE(SAD1$n[i] < SAD2$n[i])){
      numerator = numerator + SAD1$n[i]
    }else{
      numerator = numerator + SAD2$n[i]
    }
  }
  bc = 1 - (2*numerator/(numindiv1+numindiv2))
  return(bc)
}

cramer <- function(chivalue, numindiv, nsp){
  df = nsp-1
  V = sqrt(chivalue/(numindiv*df))
  return(V)
}

density <- function(SAD){
  count = 0
  for(i in 1:dim(SAD)[1]){
    count= count + SAD$n[i]
  }
  return(count)
}

presenceabscence <- function(SAD){
  for(i in 1:dim(SAD)[1]){
    if(SAD$n[i] >0){
      SAD$n[i] = 1
    }
  }
  return(SAD)
}

createSAD <-function(sample, nsp){
  sad<- as.data.frame(matrix(c(0:(nsp-1),rep(0,nsp)),nrow = nsp, ncol =2))
  colnames(sad)<- c("sp","n")
  for(i in 1:dim(sad)[1]){
    for(j in 1:length(sample)){
      if(isTRUE(sad$sp[i] == sample[j])){
        sad$n[i] = sad$n[i] +1
      }
    }
  }
  return(sad)
}


rowsalive <- function(year, bci,method){
  specieskey <-read.csv("~/Desktop/bci project/specieskeyrf.csv")
  rownames(specieskey)<-specieskey$sp
  tree <- c(0)
  ids <- c(tree)
  for(i in 1:dim(bci)[1]){
    #if(isTRUE(bci$dbh[i] > specieskey[as.character(bci$sp[i]),"rf"]) && isTRUE(bci$status[i] == "A")){
    if(isTRUE(bci$dbh[i] > 99) && isTRUE(bci$status[i] == "A")){
      ids<- rbind(ids, c(i))
    }
  }
  ids<-ids[-c(1),]
  title = paste0("~/Desktop/bci project/",method, year, "yralivetrees.csv")
  write.csv(ids,title, row.names = F)
}
