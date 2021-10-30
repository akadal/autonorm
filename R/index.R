library("GA")

require("GA")
require("doParallel")
require("sqldf")
options(sqldf.driver = "SQLite")

autonorm <- function(anRawdata,
                     reporting = T,
                     parallel = T, 
                     caching = F,
                     sampling = F,
                     w = rep(1, 4),
                     popsize = 100, 
                     maxiter = 1000,
                     pmutation = 0.2) {
  
  options("scipen"=10000, "digit"=22)
  
  # preparing data
  colnames(anRawdata) <- NULL
  
  # defining hvs
  hvs <<- ncol(anRawdata)
  
  # controling hvs
  if (hvs < 3) {
    stop ("Number of column must be 3 least.")
  }
  
  # defining m
  m <- ceiling(log2(hvs))
  if (m == 1) { m <- 2 }
  
  # defining tab
  tab <- (2^(m))
  
  # defining chsize
  chsize <- m * hvs
  
  # sampling
  if (is.numeric(sampling)) {
    
    rowCount <- NROW(anRawdata)
    sampleSize <- round(rowCount * sampling)
    anRawdata <- anRawdata[sample(rowCount, sampleSize),]
    
  }
  
  # caching
  if(caching) { anCaching() }
  
  # parallel
  if (parallel) {
    corecount <- detectCores()-1
  } else {
    corecount <- 1
  }
  
  cl <- makeCluster(corecount, type='PSOCK')
  registerDoParallel(cl)
  
  # GA !
  gaResult <- ga(
    type="binary",
    nBits=chsize,
    fitness=anFitnessFunction,
    anRawdata=anRawdata, 
    caching=caching,
    hvs = hvs,
    m = m,
    tab = tab,
    w = w,
    maxiter=maxiter,
    popSize=popsize,
    pmutation = pmutation,
    elitism = 2, 
    mutation = anMutation,
    parallel = corecount,
    monitor = anMonitor
  )
  
  stopCluster(cl)
  
  returnSet <- list()
  
  returnSet$gaResult <- gaResult
  returnSet$settings$ncol <- NCOL(anRawdata)
  returnSet$settings$nrow <- NROW(anRawdata)
  returnSet$settings$corecount <- corecount
  returnSet$settings$popsize <- popsize
  returnSet$settings$maxiter <- maxiter
  returnSet$settings$caching <- caching
  returnSet$settings$pmutation <- pmutation
  returnSet$settings$sampling <- sampling
  returnSet$weights <- w
    
  
  
  chm <- matrix(gaResult@solution[1,], ncol = m, byrow = TRUE)
  v <- apply(chm, 1, binary2decimal) + 1
  tables <- unique(v)
  
  
  returnSet$result$v <- length(v)
  returnSet$result$chm <- length(chm)
  returnSet$result$entities <- length(tables)
  
  # preparing the report
  report <- cat("\n\n### AUTONORM REPORT ###\n
Settings:
Dataset size - col:",NCOL(anRawdata),"row:",NROW(anRawdata),"
Used core count:",corecount,"
Population size:",popsize,"
Maximum iteration:",maxiter,"
Caching mode:",caching,"
Probability of mutation:",pmutation,"
Sampling mode:",sampling,"
Weights:",w,"\n
The inputted dataset has been analyzed using its ",length(v),"variables.
Autonorm suggests to create",length(tables),"entities.\n
Suggested database design:\n")
  
  for (i in 1:length(tables)) {
    
    report <- cat(report,"Entitiy",i,"has",length(which(tables[i] == v)),"variable(s):",which(tables[i] == v),"\n")
    returnSet$result$design$entity[[i]] <- which(tables[i] == v)
    
  }
  
  report <- cat(report,"* Each variable id shows column order in the given dataset.\n")
  
  if (length(tables) > 1) {
    
    report <- cat(report,"\nRelations suggestions:\n")
    
    for (i in 1:(length(tables) - 1)) {
      for (j in (i+1):length(tables)) {
        temprel <- anRelation(anRawdata, which(tables[i] == v), which(tables[j] == v))
        report <- cat(report,"[Entity ",i," and Entity ",j,"] ")
        report <- cat(report,"Relation type:",temprel['reltype']," | Success: ",round(1 / as.numeric(temprel['successrate'])),"\n")
        returnSet$result$relations$rel[[paste0(i,"-",j)]]$entity1 <- i
        returnSet$result$relations$rel[[paste0(i,"-",j)]]$entity2 <- j
        returnSet$result$relations$rel[[paste0(i,"-",j)]]$type <- temprel['reltype']
        returnSet$result$relations$rel[[paste0(i,"-",j)]]$success <- as.numeric(temprel['successrate'])
      }
    }
    
  }
  
  report <- cat(report,"\nIn order to cite AUTONORM in your paper, please check citation(\"autonorm\")\n\n")
  
  # end of the preparing the report
  
  if (reporting) { cat(report) }
  
  # if (caching) { anDestroyCache() }
  
  return (returnSet)
  
}

anFitnessFunction <- function (ch, anRawdata, caching, hvs, m, tab, w) {
  
  # caching
  if (caching) {
    
    fitval <- anCheckCache(anGetCache(),binary2decimal(ch))
    if(is.numeric(fitval)) {
      return(-fitval)
    }
    
  }
  
  # ch to table id vector
  chm <- matrix(ch, ncol = m, byrow = TRUE)
  v <- apply(chm, 1, binary2decimal) + 1
  
  # creating table vector
  tables <- unique(v)
  
  # stacking tables
  for (i in tables) {
    temptable <- unique(anRawdata[(which (v == i))])
    assign(paste("table",i,sep=""), temptable)
  }
  
  # objective 1
  totalcell <- 0
  for (i in tables) {
    temptable <- get(paste("table",i,sep=""))
    totalcell <- totalcell + ((ncol(temptable)) * nrow(temptable))
  }
  obj1 <- w[1] *  (totalcell / (nrow(anRawdata)*ncol(anRawdata)))
  
  # objective 2
  dupprates <- c()
  for (i in 1:length(v)) {
    mytable <- get(paste("table", v[i], sep="", collapse = ""))
    thecol <- mytable[,paste("V", i, sep="", collapse = "")]
    allrows <- length(thecol)
    distinctrows <- length(unique(thecol))
    myvalue <- ((allrows - distinctrows) / allrows)
    dupprates <- append(dupprates,myvalue)
  }
  obj2 <- w[2] * (sum(dupprates,na.rm = TRUE))
  
  # objective 3
  obj3 <- w[3] * (length(tables) / tab)
  
  # objective 4
  numericaltable <- 0
  for (i in tables) {
    if (is.numeric(as.vector(as.matrix(anRawdata[which(v == i)])))) {
      numericaltable <- numericaltable + 1
    }
  }
  obj4 <- w[4] * ((numericaltable / length(unique(v))) ^ (1/length(unique(v))))

  # summing weighted objectives
  fitnessvalue <- sum(    obj1 + 
                            obj2 + 
                            obj3 + 
                            obj4,
                          na.rm = TRUE  )
  
  if (caching) { anCaching(c(as.character(binary2decimal(ch)), as.character(fitnessvalue))) }
  
  return (-fitnessvalue)
}

anMutation <- function (object, parent, ...)  {
  mutate <- parent <- as.vector(object@population[parent, ])
  mymnumber <- hvs
  mycoin <- round(runif(1,1,5))
  if (mycoin == 1) {
    n <- length(parent)
    j <- sample(1:n, size = 1)
    mutate[j] <- abs(mutate[j] - 1)
  } else if (mycoin == 2) {
    chmatrix <- matrix(mutate,ncol = mymnumber,byrow = TRUE)
    chmatrix[round(runif(1,1,nrow(chmatrix))),] <- round(runif(mymnumber,0,1))
    mutate <- as.vector(t(chmatrix))
  } else if (mycoin == 3) {
    chmatrix <- matrix(mutate,ncol = mymnumber,byrow = TRUE)
    chmatrix2 <- chmatrix
    gene1 <- round(runif(1,1,nrow(chmatrix)))
    gene2 <- round(runif(1,1,nrow(chmatrix)))
    chmatrix[gene1,] <- chmatrix2[gene2,]
    chmatrix[gene2,] <- chmatrix2[gene1,]
    mutate <- as.vector(t(chmatrix))
  } else if (mycoin == 4) {
    chmatrix <- matrix(mutate,ncol = mymnumber,byrow = TRUE)
    chmatrix[round(runif(1,1,nrow(chmatrix))),] <- chmatrix[round(runif(1,1,nrow(chmatrix))),]
    mutate <- as.vector(t(chmatrix))
  } else if (mycoin == 5) {
    chmatrix <- matrix(mutate,ncol = mymnumber,byrow = TRUE)
    genepool <- unique(chmatrix)
    for (i in 1:nrow(chmatrix)) {
      chmatrix[i,] <- genepool[round(runif(1,1,nrow(genepool))),]
    }
    mutate <- as.vector(t(chmatrix))
  }
  return(mutate)
}

anMonitor <- function (object, digits = getOption("digits"), ...) {
  fitness <- na.exclude(object@fitness)
  sumryStat <- c(mean(fitness), max(fitness))
  sumryStat <- format(sumryStat, digits = digits)
  cat(paste("\rAUTONORM | Completed: ", round((object@iter / object@maxiter)*100) ,"% | Minimum Fitness Value = ", (-1 * as.numeric(sumryStat[2])),sep = ""))
  flush.console()
}

anRelation <- function(rawdata,set1,set2) {
  
  set1table <- rawdata[set1]
  set2table <- rawdata[set2]
  
  # soyut tablo yaratılıyor
  abstracttable <- as.data.frame(matrix(c(apply(set1table,1,paste,collapse=";"),apply(set2table,1,paste,collapse=";")),ncol = 2))
  
  # satır sayısı belirleniyor
  satirsay <- NROW(abstracttable)
  
  # geçiş sayıları belirleniyor
  gecis1 <- sqldf("SELECT count(distinct(V1)) FROM abstracttable GROUP BY V2 HAVING count(distinct(V1)) > 1")
  
  if (NROW(gecis1) > 0) {
    gecis1num <- NROW(gecis1)
  } else {
    gecis1num <- 0
  }
  
  gecis2 <- sqldf("SELECT count(distinct(V2)) FROM abstracttable GROUP BY V1 HAVING count(distinct(V2)) > 1")
  
  if (NROW(gecis2) > 0) {
    gecis2num <- NROW(gecis2)
  } else {
    gecis2num <- 0
  }
  
  # ilişki türü belirleniyor
  if (gecis1num == 0 && gecis2num == 0) {
    iliski <- "1-1"
  } else if (gecis1num > 0 && gecis2num == 0) {
    iliski <- "1-m"
  } else if (gecis1num == 0 && gecis2num > 0) {
    iliski <- "m-1"
  } else if (gecis1num > 0 && gecis2num > 0) {
    iliski <- "m-m"
  }
  
  # başarı puanı hesaplanıyor
  bp <- (satirsay/(satirsay + gecis1num^2))*(satirsay/(satirsay + gecis2num^2))
  
  # dönüş seti hazırlanıyor
  result <- c()
  result["transition1"] <- gecis1num
  result["transition2"] <- gecis2num
  result["reltype"] <- iliski
  result["successrate"] <- bp
  
  return (result)
}

anCaching <- function(data = NULL) {
  
  if (is.null(data)) {
  
    cat(c("key","val"),"\n",file = ".caching.tmp",append = F,sep = ";")
    
  } else {
    
    cat(data,"\n",file = ".caching.tmp",append = T,sep = ";")
    
  }
  
}

anGetCache <- function() {
  df <- read.csv(".caching.tmp", sep = ";")
  return(df)
}

anCheckCache <- function(df, key) {
  indices <- which(df$key == key)
  if (length(indices) == 0) {
    return(F)
  } else {
    return(as.numeric(df[[indices[1],"val"]]))
  }
}

anDestroyCache <- function() {
  file.remove(".caching.tmp")
}
