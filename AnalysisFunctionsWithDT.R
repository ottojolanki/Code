require(data.table)
require(dplyr)

customWrite <- function(DF, filename){
  write.table(DF, filename, sep = "\t", row.names = FALSE)
}
##function for cleaning hsbc(after getting the correct subjects in the hsbc)

cleanHsbc <- function(hsbcDF){
  ##get the hsbc right
  hsbcDF$vhash <- NULL
  hsbcDF$jhash <- NULL
  hsbcDF$chain <- NULL
  hsbcDF <- as.data.frame(t(hsbcDF))
  names(hsbcDF) <- c("age", "dbid", "gender", "positive_602", "positive_0301", "beforeafter", "casecontrol", "cell_type", "rs1483979")
  hsbcDF <- hsbcDF[-1,]
  hsbcDF$age <- as.numeric(as.character(hsbcDF$age))
  hsbcDF$gender <- factor(as.character(hsbcDF$gender))
  hsbcDF$positive_602 <- factor(as.character(hsbcDF$positive_602))
  hsbcDF$positive_0301 <- factor(as.character(hsbcDF$positive_0301))
  hsbcDF$casecontrol <- factor(as.character(hsbcDF$casecontrol))
  hsbcDF$cell_type <- factor(as.character(hsbcDF$cell_type))
  hsbcDF$rs1483979 <- factor(as.character(hsbcDF$rs1483979))
  return(hsbcDF)
}

##function for identifying two types of alpha-J24s
identifyJ24 <- function(DF = data.frame()){
  J24.df <- subset(DF,jhash=="TRAJ24*01")
  not.J24.df <- subset(DF, jhash!="TRAJ24*01")
  #get ones that match the "SWGK[FL][RQ]$" pattern
  J24.df <- J24.df[grepl("SWGK[FL][RQ]$", J24.df$covariate),]
  J24.df$J24type <- sapply(J24.df$covariate, function(x){ifelse(grepl("L[RQ]$", x), "01", "02")}, USE.NAMES = FALSE, simplify = TRUE)
  J24.df$jhash <- paste(J24.df$jhash, J24.df$J24type, sep = ":")
  J24.df$J24type <- NULL
  return(rbind(J24.df, not.J24.df))
}

##analysis functions

normalize<-function(mat=matrix()){
  col.sums<-colSums(mat)
  if(prod(col.sums)){
    t(t(mat)/col.sums)*min(col.sums)
  }
  else{
    stop("some columns sum to zero")
  }
}
removeRareClones<-function(m,N){
  m[which(rowSums(apply(m[,-1],c(1,2),function(x)x!=0))>N),]
}

wrapNormalize <- function(DF){
  mat <- as.matrix(sapply(DF[, -1], as.numeric))
  mat <- normalize(mat)
  return(cbind(name = DF[, 1], as.data.frame(mat)))
}

##does not order the result 
doLogmodel <- function(DF, hsbc = NULL){  
  name <- DF$name
  DF$name <- NULL
  bigmodel <- lm(log(t(DF)+1) ~ hsbc$age + hsbc$gender + hsbc$positive_0301 + hsbc$casecontrol + hsbc$rs1483979)
  bigsummary <- summary(bigmodel)
  coef     <- lapply(bigsummary, function(x)x$coefficients[,4])
  eff      <- lapply(bigsummary,function(x)x$coefficients[,1])
  coefDF   <- as.data.frame(coef)
  effDF    <- as.data.frame(eff)
  minpvals <- apply(coefDF,2,min)
  coefDF   <- rbind(coefDF,minpvals)
  names(coefDF) <- name 
  names(effDF) <- name 
  coefDF <- rbind(coefDF,effDF)
  return(coefDF)
}

##simple log(x+1) transformed model with just case control status as covariate

doSimpleLogmodel <- function(DF, hsbc = NULL){  
  name <- DF$name
  DF$name <- NULL
  bigmodel <- lm(log(t(DF)+1) ~ hsbc$casecontrol)
  bigsummary <- summary(bigmodel)
  coef     <- lapply(bigsummary, function(x)x$coefficients[,4])
  eff      <- lapply(bigsummary,function(x)x$coefficients[,1])
  coefDF   <- as.data.frame(coef)
  effDF    <- as.data.frame(eff)
  minpvals <- apply(coefDF,2,min)
  coefDF   <- rbind(coefDF,minpvals)
  names(coefDF) <- name 
  names(effDF) <- name 
  coefDF <- rbind(coefDF,effDF)
  return(coefDF)
}

##calculate standard error of a vector
std <- function(x) sd(x)/sqrt(length(x))

##function to perform wilcox test for case/control groups according to a character vector. By default perform an unpaired wilcoxon test.
doWilcoxtest<-function(DF = NULL, casecontrol = NULL, paired = FALSE  ){
  name       <- NA
  testp      <- NA
  meancase <- NA
  SEMcase  <- NA
  SEMctrl   <- NA
  meanctrl  <- NA
  
  tDF   <- as.data.frame(t(DF[, -1]))
  nimet <- as.character(DF[, 1])
  
  for (i in 1:ncol(tDF)) {
    testp[i]    <- wilcox.test(tDF[, i] ~ casecontrol, paired = paired)$p.value
    name[i]     <-nimet[i]
    meancase[i] <-mean(tDF[which(casecontrol == "case") ,i])
    meanctrl[i] <-mean(tDF[which(casecontrol == "control") ,i])
    SEMcase[i]  <-std(tDF[which(casecontrol == "case") ,i])
    SEMctrl[i]  <-std(tDF[which(casecontrol == "control") ,i])
    
  }
  
  bound <- as.data.frame(cbind(name, testp, meancase, SEMcase, meanctrl, SEMctrl))
  
  return(bound)
}

##countNonZeroClones function:
##count the number of cases and controls for which the clone count is non-zero
##input: dataframe with names of clones on the first column, and numeric data in the rest of the columns
##       character vector indicating cases ("case") and controls ("control") 
countNonZeroClones <- function(DF, casecontrol = character()){
  case    <- which(casecontrol == "case")
  control <- which(casecontrol == "control")
  case.nonzero    <- rowSums(apply(DF[, -1][, case], c(1,2), function(x){x != 0}))
  control.nonzero <- rowSums(apply(DF[, -1][, control], c(1,2), function(x){x != 0}))
  names(case.nonzero)    <- NULL
  names(control.nonzero) <- NULL
  return(data.frame(casenonzeroN = case.nonzero, ctrlnonzeroN = control.nonzero))
}

##calculateCaseControlMeans function:
##calculate mean counts of cases and controls indicated by character vector in dataframe with names of clones in the first column and numeric
##data in the rest of the columns.
##outputs a two column dataframe with columns case.mean and column.mean containing the numerical mean count values
calculateCaseControlMeans <- function(DF, casecontrol = character()){
  case.indices    <- which(casecontrol == "case")
  control.indices <- which(casecontrol == "control")
  case.n    <- length(case.indices)
  control.n <- length(control.indices)
  case.mean    <- rowSums(DF[, -1][, case.indices]) / case.n
  control.mean <- rowSums(DF[, -1][, control.indices]) / control.n
  return(data.frame(case.mean = case.mean, control.mean = control.mean))
}


##buildContingencyTables
##input: numerical vectors of same length x, and y, indicating the number of present condition, 
##numbers totalX and totalY indicating the total amounts of subjects
##(cases vs controls for example).
##builds a list (of length length(x)=length(y)) of tables ready-for-use for chisq.test function in {base}

buildContingencyTables <- function(x, y, totalX, totalY){
  if ((length(x) != length(y)) | (length(x)*length(y) == 0)){
    stop("the input vectors have to be of equal, nonzero, length")
  }
  result = list()
  for (i in 1:length(x)){
    xpresent <- x[[i]]
    xabsent  <- totalX - x[[i]]
    ypresent <- y[[i]]
    yabsent  <- totalY - y[[i]]
    t <- as.table(rbind(c(xpresent, xabsent), c(ypresent, yabsent)))
    dimnames(t) <- list(cactrl = c("case", "control"), presence = c("present", "absent"))
    result[[i]] <- t
  }
  return(result)
}