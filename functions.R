#Handy funcitons used for all converting of samplings for FISHECU


#fix name of column, print error if there are more than 1 column with pattern specified
fixColumnName <- function(x, xpattern, new.name){
  colidx <- grep(pattern = xpattern, x = names(x)); 
  if(! new.name %in% names(x)){
    if(length(colidx) > 1) print(paste("ERROR:more columns starting with pattern -", xpattern));
    if(length(colidx) == 0){x[,new.name] <- NA}else{names(x)[colidx] <- new.name};
  }
  return(x)  
}

#fix species
fixSpecies <- function(x, species.typos, species.db){
  x$Species <- sub("\\\u009d$", "", x$Species)
  #load map of species - first correct typos
  x <- merge(x, species.typos, by  = "Species", all.x = T)
  x[!is.na(correctsp), Species := correctsp]
  sort(unique(x$Species))
  species.unq <- sort(unique(x$Species))
  x <- merge(x, species.db[, .(czechcodename, czechname)], by.x = c("Species"), by.y= c("czechcodename"), all.x = T)
  #CHECKs
  #get species which are not in species table
  print(paste("Following species were not found in database:" , paste(unique(x[is.na(czechname)]$Species), collapse = ","))) # '0' OK
  x
}

#function to get for each input time if it was night or day
#input:
# x - vector of type POSIXct
# lat - numeric vector of lenght one
# lon - numeric vector of lenght one
# sunsetsunrise <- suncalc::getSunlightTimes(date = as.Date(c("2017-05-02", "2017-05-05", "2017-05-06")), lat = 35, lon = 32, keep = c("sunrise", "sunset"))[,c(4,5)]
# sunrise.sunset.t <- data.table(sunrise.set(49+(34/60)+((42/60)/60), 15+(15/60)+((21/60)/60), "2017-01-01" , num.days = 365, timezone = "UTC"))
# sunrise.sunset.t[, rize.sec := as.numeric(difftime(sunrise, round.POSIXt(x = sunrise, "days"), units = "secs"))]
# sunrise.sunset.t[, set.sec := as.numeric(difftime(round.POSIXt(x = sunset, "days"), sunset, units = "secs"))]
# 
# thresholds <- sunrise.sunset.t[, .(first.date = min(as.Date(sunrise)), minrise = min(rize.sec), maxrise = max(rize.sec), minset = min(set.sec), maxset = max(set.sec)), by = .(month = month(sunrise))] 
# thresholds[, data.table(period = rep(c("eveningshift", "night", "morningshift"), each = 4), 
#                         x = c(minset, minset, maxset, maxset,
#                               maxset, maxset, minrise, minrise, 
#                               minrise, minrise, maxrise, maxrise),
#                         y = rep(c(-Inf, Inf, Inf, -Inf), 3)), by = .(first.date)]
# 



getNightDay <- function(x, lat = 49.5765639, lon = 14.6637706, xtz = "UTC", returnSunrSuns = F, day = "day", night = "night"){
  if(length(x) == 0) stop("length of input is 0")
  if(all(is.na(x))) return(as.character(x))
  
  library(StreamMetabolism)
  library(data.table)
  from <- min(x, na.rm = T)
  to <- max(x, na.rm = T)
  sunrise.sunset <- as.data.table(sunrise.set(lat, lon, from , num.days = difftime(time1 = to, time2 = from, units = "days", tz = xtz)+2, timezone = xtz))
  sunrise.sunset$Date <- as.Date(sunrise.sunset$sunrise, tz = xtz)
  res <- merge(data.table(x, xorder = 1:length(x), Date = as.Date(x)), sunrise.sunset, by = "Date", all.x = T)
  res[, Diel.period := as.character(ifelse(x > sunrise & x < sunset, "day", "night"))]
  setkey(res, xorder)
  if(returnSunrSuns){
    return(res[, .(Diel.period, sunrise, sunset)])
  }else{
    return(res$Diel.period)  
  }
}



#recursive function to get samplings close to each other
getOvelappingSamplings<- function(samplingids, dates, days.distance = 10){
  library(data.table)
  if(length(samplingids) == 0 | length(samplingids) == 0) {print("ERROR: vector of length 0 as input"); return(NA)}
  x <- data.table(samplingid = samplingids, Date = dates)
  x$processed <- F
  x$group <- -1
  group.i <- 1
  for(i in 1:nrow(x)){
    if(x$processed[i] == F){
      oversamplings.idxs <- getOverSamplingsRec(i, x, days.distance)  
      x[oversamplings.idxs, processed := T]
      x[oversamplings.idxs, group := group.i]
      group.i <- group.i + 1
    }  
  }
  return(x[, .(group)])
}

getOverSamplingsRec <- function(idx, aa, days.distance){
  aa$processed[idx] <- T
  idxs <- unique(which(aa$processed == F & aa$Date %between% c(aa[idx]$Date-days.distance, aa[idx]$Date+days.distance)))
  if(length(idxs) == 0) return(idx)
  aa$processed[idxs] <- T
  return(unique(c(idx, unlist( sapply(unique(idxs), function(idx, aa, days.distance) getOverSamplingsRec(idx, aa, days.distance), days.distance = days.distance, aa = aa, simplify = T), recursive = TRUE))))
}

#Fix chzech letters - substitute special czech letters
#@input: character vector to be cleaned
#@output: cleaned character vector
fixCzechLetters <- function(x){
  #Vowels
  x <- gsub(pattern = "\u00E1", replacement = "a", x )#Lowercase long a (á)
  x <- gsub(pattern = "\u00C1", replacement = "A", x )#Uppercase long a (Á)
  x <- gsub(pattern = "\u00E9", replacement = "e", x )#Lowercase long e (é)
  x <- gsub(pattern = "\u00C9", replacement = "E", x )#Uppercase long e (É)
  x <- gsub(pattern = "\u011B", replacement = "e", x )#Lowercase soft e (ě)
  x <- gsub(pattern = "\u011A", replacement = "E", x )#Uppercase soft e (Ě)
  x <- gsub(pattern = "\u00ED", replacement = "i", x )#Lowercase long i (í)
  x <- gsub(pattern = "\u00CD", replacement = "I", x )#Uppercase long i (Í)
  x <- gsub(pattern = "\u00F3", replacement = "o", x )#Lowercase long o (ó)
  x <- gsub(pattern = "\u00D3", replacement = "O", x )#Uppercase long o (Ó)
  x <- gsub(pattern = "\u00FA", replacement = "u", x )#Lowercase long u (ú)
  x <- gsub(pattern = "\u00DA", replacement = "U", x )#Uppercase long u (Ú)
  x <- gsub(pattern = "\u016F", replacement = "u", x )#Lowercase ring u (ů)
  x <- gsub(pattern = "\u016E", replacement = "U", x )#Uppercase ring u (Ů)
  #Consonants
  x <- gsub(pattern = "\u010D", replacement = "c", x )#Lowercase soft c (č)
  x <- gsub(pattern = "\u010C", replacement = "C", x )#Uppercase soft c (Č)
  x <- gsub(pattern = "\u010F", replacement = "d", x )#Lowercase soft d (ď)
  x <- gsub(pattern = "\u010E", replacement = "D", x )#Uppercase soft d (Ď)
  x <- gsub(pattern = "\u0148", replacement = "n", x )#Lowercase soft n (ň)
  x <- gsub(pattern = "\u0147", replacement = "N", x )#Uppercase soft n (Ň)
  x <- gsub(pattern = "\u0159", replacement = "r", x )#Lowercase soft r (ř)
  x <- gsub(pattern = "\u0158", replacement = "R", x )#Uppercase soft r (Ř)
  x <- gsub(pattern = "\u0161", replacement = "s", x )#Lowercase soft s (š)
  x <- gsub(pattern = "\u0160", replacement = "S", x )#Uppercase soft s (Š)
  x <- gsub(pattern = "\u0165", replacement = "t", x )#Lowercase soft t (ť)
  x <- gsub(pattern = "\u0164", replacement = "T", x )#Uppercase soft t (Ť)
  x <- gsub(pattern = "\u00FD", replacement = "y", x )#Lowercase long y (ý)
  x <- gsub(pattern = "\u00DD", replacement = "Y", x )#Uppercase long y (Ý)
  x <- gsub(pattern = "\u017E", replacement = "z", x )#Lowercase soft z (ž)
  x <- gsub(pattern = "\u017D", replacement = "Z", x )#Uppercase soft z (Ž)
  return(x)
}

#Function to detect outliers in LW relationship  
detectOuts <- function(W, L, rlm.thresh){
  library(MASS)
  if(length(W) >1 & length(L) > 1){
    m.rlm <- rlm(log10(W) ~ log10(L), maxit = 200)
    return(m.rlm$w < rlm.thresh)
  }else{return(F)}
}

#general function computing value per unit of effort
#input - dt of samplings, catch, vector of catch split factors, vector of sampling split factors, samplingid column, value.var is name of column containing count/weight/whatever

getVPUE <- function(samplings, catch, split.factors.catch, 
                    split.factors.samplings, id.colname = "samplingid", effort.colname = 
                      "Effort", returnlist = F, value.var){
  #check input
  library(data.table)
  if(!is.character(split.factors.catch) | 
     !is.character(split.factors.samplings)){stop("Split.factors argument is 
                                                  expected to be character vector")  }
  if(!(is.data.table(samplings) & is.data.table(catch))){stop("Catch 
                                                              and sampling has to be data.table object")}
  if(length(id.colname) != 1){stop("Wrong length of id.colname 
                                   argument: expected length == 1")}
  if(!value.var %in% names(catch)){stop(paste("The value.var:", 
                                              value.var, " column is missing in catch argument"))}
  if(!id.colname %in% names(samplings)){stop(paste("The id.colname:", 
                                                   id.colname, " column is missing in samplings argument"))}
  if(!id.colname %in% names(catch)){stop(paste("The id.colname:", 
                                               id.colname, " column is missing in catch argument"))}
  if(!effort.colname %in%  names(samplings)){stop(paste("The 
                                                        effort.colname:", effort.colname, " column is missing in samplings 
                                                        argument"))}
  
  #check if there are some fish
  if(nrow(catch) == 0){warning("There was no fish caught (this might 
                               mean that sampling ids are not matching!)")}
  
  #dataframes should
  
  #TODO: check if all split factors are present in tables
  #add samplingid column for easy reference
  catch$samplingid <- catch[, id.colname, with = F]
  samplings$samplingid <- samplings[, id.colname, with = F]
  samplings$Effort <- samplings[, effort.colname, with = F]
  catch$value.var <- catch[, value.var, with = F]
  
  #subset catch only to catch from sampling table
  catch <- catch[samplingid %in% samplings$samplingid]
  
  #split into one vector
  split.factors <- c(split.factors.catch, split.factors.samplings)
  
  #keep only columns needed for computation
  catch <- catch[, c("samplingid", "value.var", split.factors.catch), 
                 with = F]
  samplings <- samplings[, c("samplingid","Effort", 
                             split.factors.samplings), with = F]
  
  #add number of samplings in factor group
  samplings[, sampcount := .N, by = split.factors.samplings]
  
  #if there are no catch factors
  if(length(split.factors.catch) > 0){
    catch.uniques <- do.call(what = CJ, args =  c(catch[, 
                                                        split.factors.catch, with = F], list(unique = T)))
    #add all species
    samplings.sp <- CJ.tables(samplings, catch.uniques)
  }else{ samplings.sp <- samplings}
  
  #if all fish dont have weight - return error
  catch.sum <- catch[, .(value.var.sum = sum(value.var), fishcount = 
                           .N), by = c("samplingid", split.factors.catch)]
  catch.sampl <- merge(catch.sum, samplings.sp, all = T, by = 
                         names(samplings.sp)[names(samplings.sp) %in% names(catch.sum)])
  catch.sampl[is.na(value.var.sum), value.var.sum :=  0]
  catch.sampl[is.na(fishcount), fishcount :=  0]
  catch.sampl[, VPUE := value.var.sum / Effort]
  
  #populate sampling table by factor values
  catch.agg <- catch.sampl[, .(VPUE.mean = mean(VPUE), VPUE.se = plotrix::std.error(VPUE), 
                               VPUE.sum = sum(value.var.sum),  sampcount = sampcount[1], 
                               catch.rowcount = sum(fishcount)), by = split.factors]
  
  names(catch.agg) <- gsub(x= names(catch.agg), pattern = "VPUE", 
                           replacement = value.var)
  #0 where the sampling was done
  #NA where the sampling was not done
  if(returnlist == T){
    return(list(VPUE.per.net = catch.sampl[, c("samplingid", 
                                               split.factors, "VPUE"), with = F], VPUE.per.factors  = catch.agg ))
  }else{
    return(catch.agg)
  }
  }



#function to cross join two data.tables
CJ.tables <- function(X,Y){ 
  unique_name <- last(make.unique(c(colnames(X),colnames(Y),"temporalColumnForCJtables"))) 
  return(X[,c(setNames(1,unique_name),.SD)][Y[,c(setNames(1,unique_name),.SD)],on=unique_name,allow.cartesian=TRUE]) 
}
#round posixct time to specific bin in seconds
round.POSIXct <- function(x, bin.secs){
  if(class(x)[1] != "POSIXct") stop("parameter x has to be type POSIXct")
  library(lubridate)
  as.POSIXct(round(as.numeric(x)/bin.secs)*bin.secs,origin= "1970-01-01", tz = tz(x))
}
floor.POSIXct <- function(x, bin.secs){
  if(class(x)[1] != "POSIXct") stop("parameter x has to be type POSIXct")
  library(lubridate)
  as.POSIXct(floor(as.numeric(x)/bin.secs)*bin.secs,origin= "1970-01-01", tz = tz(x))
}
ceiling.POSIXct <- function(x, bin.secs){
  if(class(x)[1] != "POSIXct") stop("parameter x has to be type POSIXct")
  library(lubridate)
  as.POSIXct(ceiling(as.numeric(x)/bin.secs)*bin.secs,origin= "1970-01-01", tz = tz(x))
}

#Functionto merge 12 and 4 size mesh panel into 16
#load 12
#load 4
#check if stratum, locality is the same for the pairs, 
#check if 
#divide 4 by 4
#union fish
#take date from 12
#leave GPS, notes, depthmin, depthmax etc. as it is. Its up to every one to 


#Shanon index
#Function computing shanon index 
#input: vector defining a species, vector defining abundances of each species(optional)
#if the abundance is not supplied, abundance is computed as number of occurences of each species in vector species
#NOTE: one species can occure in species vector more than once (can be used with pure catch table)
#output:Shanon index
#Example:
#catch[, .(Shannon.idx = computeShannon(species, abundance), by = .(reservoir, locality)]

computeShannon <- function(species,  abundance = NULL){
  if(is.null(abundance)){
    x <- data.table(species = species)
    x.sp <- x[, .(spcount = .N), by = species]
  }else{
    x <- data.table(species = species, abundance = abundance)
    x.sp <- x[, .(spcount = sum(abundance)), by = species]
  }
  x.sp[, totcount := sum(spcount)]
  x.sp[,  p.sp := spcount / totcount]
  return(x.sp[, -sum(p.sp*log(p.sp))])
}

#roll join function by constant width time window
#input: x - vector of values to be smoothed
#time.vec - vector of timestamps
#win.widt - width of the window in seconds
#FUN - function to apply on values in the window
rollTimeWindow <- function(x, time.vec, win.width, FUN){
  x.out <- vector(mode = class(x), length = length(x))
  x.out[1:length(x.out)] <- NA
  if(is.unsorted(time.vec)){
    stop("Time vector is not in increasing order!")
  }
  FUN <- match.fun(FUN)
  time.vec.minus <- time.vec-win.width
  time.vec.plus <- time.vec+win.width
  for(i in 1:length(x.out)){
    x.out[i] <- FUN(x[which(time.vec > time.vec.minus[i] & time.vec < time.vec.plus[i])])
  }
  return(x.out)
}

#Formula to change NaN to 0 in data tables#### 
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

#Negative Binomial GLM test and summary with 6 different pseudo-r2 computations CPUE
#Adding stars as symbols for the statistical significance

star <- function(pval) {
  if (pval <= 0.001) {
    return("***")
  }
  if (pval <= 0.01) {
    return("**")
  }		
  if (pval <= 0.05) {
    return("*")
  }
  else {return("non-significant")
  }
}
glm_nb_summary <- function(dataset, formula){
  glm_nb <- glm.nb(formula = formula, data = dataset)
  r2_kl <- round(rsq(fitObj = glm_nb, adj = F, type = c("kl")), 2)
  r2_zhang <- round(rsq(glm_nb),2)
  pseudo <- modEvA::RsqGLM(glm_nb)
  r2_CoxSnell <- round(pseudo$CoxSnell,2)
  r2_Nagelkerke <- round(pseudo$Nagelkerke,2)
  r2_McFadden <- round(pseudo$McFadden,2)
  r2_sqPearson <- round(pseudo$sqPearson,2)
  glm_aic <- round(AIC(glm_nb), 2)
  mnull <- update(glm_nb, . ~ 1)#Creating a null model
  fit.test <- round(pchisq(2 * (logLik(glm_nb) - logLik(mnull)), df = 1, lower.tail = FALSE), 2)#Comparing the null model with the full model. Attention to the df, it varies according to the number of predictors used
  return(list(Null_mdl_comparison = fit.test,
              r2_zhang = r2_zhang,
              r2_kl_divergence_based = r2_kl,
              r2_CoxSnell = r2_CoxSnell,
              r2_Nagelkerke = r2_Nagelkerke,
              r2_McFadden = r2_McFadden,
              r2_sqPearson = r2_sqPearson,
              aic = glm_aic,
              intercept = round(coef(glm_nb)[1],2),
              slope_year = round(coef(glm_nb)[2],2),
              p_year = round(coef(summary(glm_nb))[2,4],2),
              significance_year = star(round(coef(summary(glm_nb))[2,4],2)),
              slope_protection = round(coef(glm_nb)[3],2),
              p_protection = round(coef(summary(glm_nb))[3,4],2),
              significance_protection = star(round(coef(summary(glm_nb))[3,4],2))))
}

#General functions used for BOLEN project
Sys.setenv(TZ='UTC')

getDetectionDataMultipleReader <- function(detraw){
  #EXTRACT#
  #get rid of more than one blank space
  detraw <- gsub(x= detraw, pattern = " +", replacement = " ")
  detraw.df <- data.table(detraw)
  #mark tag detection lines
  detraw.df[, valid := grepl(pattern = "^D", x = detraw)]
  
  #get datetime
  detraw.df[valid == T, datetime := sub(x = detraw, replacement = "\\1", pattern = "^D ([0-9]+-[0-9]+-[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+) .*$")]
  #get number of detections
  detraw.df[valid == T, number_of_detections := sub(x = detraw, pattern = ".* A[0-9]+ ([0-9]+) ([0-9]+|\\.)$", replacement = "\\1")]
  #get detection time
  detraw.df[valid == T, detection.time := sub(x = detraw, pattern = "^D [0-9]+-[0-9]+-[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ ([0-9]+:[0-9]+:[0-9]+\\.[0-9]+) .*$", replacement = "\\1")]
  #get antenna
  detraw.df[valid == T, antenna := paste("M", sub(x = detraw, pattern = ".* A([0-9]+) [0-9]+ ([0-9]+|\\.)$", replacement = "\\1"), sep = "")]
  
  #REMOVE SECTION#
  #get tags - number between "bigChar number trailingnumber1 trailingnumber2 end"  
  detraw.df[valid == T, tag.number := sub(x = detraw, replacement = "\\1", pattern = ".* [A-Z]+ (.+) A[0-9]+ [0-9]+ ([0-9]+|\\.)$")]
  
  #remove detections with chinese (not ascii character)
  detraw.df[grepl(x = tag.number, pattern = "[^\x20-\x7E]"),  valid := F]
  #remove detections with any character (not a number) but not the R or A at the end of tag_number
  detraw.df[grepl(x = tag.number, pattern = "^.+[a-zA-Z]$"), valid := F]
  return(detraw.df)
}

#Function parsing antenna files
#'@param detraw Character vector containing lines of antenna logs
#'@return alksdjfl
getDetectionDataSingleReader <- function(detraw){
  require(data.table)
  #EXTRACT#
  #get rid of more than one blank space
  detraw <- gsub(x= detraw, pattern = " +", replacement = " ")
  detraw.df <- data.table(detraw)
  #the logs look different during antenna download - first is tag sn
  detraw.df[, down.detection := grepl(pattern = "^>*[A-Z] .+ [0-9]+/[0-9]+/[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ .*$", x = detraw)]
  #get datetime
  detraw.df[, datetime := sub(x = detraw, replacement = "\\1", pattern = "^>*([0-9]+/[0-9]+/[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+) .*$")]
  detraw.df[down.detection ==T, datetime := sub(x = detraw, pattern = "^>*[A-Z] .+ ([0-9]+/[0-9]+/[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+) [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ .*$", replacement = "\\1")]
  #get number of detections
  detraw.df[down.detection ==F, number_of_detections := sub(x = detraw, pattern = ".*[A-Z] .+ ([0-9]+) [0-9]+$", replacement = "\\1")]
  #get detection time
  detraw.df[down.detection ==F, detection.time := sub(x = detraw, pattern = "^>*[0-9]+/[0-9]+/[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ ([0-9]+:[0-9]+:[0-9]+\\.[0-9]+) .*$", replacement = "\\1")]
  detraw.df[down.detection ==T, detection.time := sub(x = detraw, pattern = "^>*[A-Z] .+ [0-9]+/[0-9]+/[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ ([0-9]+:[0-9]+:[0-9]+\\.[0-9]+) .*$", replacement = "\\1")]
  
  #REMOVE SECTION#
  detraw.df$valid <- F
  #remove lines which are not containing detections - length of detection is used for this
  detraw.df[down.detection == T | grepl(x = detraw, pattern = "^>*[0-9]+/[0-9]+/[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ .*$"), valid := T]
  #get tag serial number
  #get tags - number between "bigChar number trailingnumber1 trailingnumber2 end"  
  detraw.df$tag.number <- sub(x = detraw.df$detraw, replacement = "\\1", pattern = ".*[A-Z] (.+) [0-9]+ [0-9]+$")
  #some antenna logged 900_... logs very strangly : 04/12/2018 09:10:47.16  00:00:00.10 A 00000 0 900 230000023879    2     9
  detraw.df[grepl(pattern = "^[0-9]+/[0-9]+/[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ [A-Z] [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+$", x = detraw), 
            tag.number :=  sub(x = detraw, replacement = "\\1_\\2",  pattern = "^[0-9]+/[0-9]+/[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ [A-Z] [0-9]+ [0-9]+ ([0-9]+) ([0-9]+) [0-9]+ [0-9]+$")]
  #for detection during download
  detraw.df[down.detection == T, tag.number :=  sub(x = detraw, replacement = "\\1", pattern = "^>*[A-Z] (.+) [0-9]+/[0-9]+/[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ [0-9]+:[0-9]+:[0-9]+\\.[0-9]+ .*$")]
  
  #remove detections with chinese (not ascii character)
  detraw.df[grepl(x = tag.number, pattern = "[^\x20-\x7E]"),valid := F]
  #remove detections with any character (not a number) at the end of tag_number
  detraw.df[grepl(x = tag.number, pattern = "^.+[a-zA-Z]"), valid := F]
  return(detraw.df)
}



#New function####
#calculations of size at age (Function!)
#input: age - num[1] age of the fish
#age.length - num[1] length of fish at age "age"
#increments - increments in percentages increments[i] is increment from age i to age i+1
getAgeLength <- function(age, age.length, increments){
  #create new vector of lengths (fish.lengths[i] is length of fish at age i)
  fish.lengths <- vector(length = length(increments))
  #assign knows size
  fish.lengths[age] <- age.length
  #backcalculate length
  for(i in age:1){
    if(i >1){
      fish.lengths[i-1] <- fish.lengths[i]/(increments[i-1]/100 + 1)
    }
  }
  #calculate future length
  for(i in age:(length(increments)-1)){
    if(i < (length(increments))){
      fish.lengths[i+1] <- fish.lengths[i]+fish.lengths[i] * (increments[i]/100)
    }
  }
  #export
  return(fish.lengths)
}



fillGender <- function(gender, presence){
  #create new vector of gender
  gender.new <- vector(length = length(gender))
  i <- 1
  
  while(is.na(gender[i]) & i <= length(gender)){ #skip year before first capture
    gender.new[i] <- NA
    i <- i + 1
  }
  #if this is the end of vector, return it
  if(i >= length(gender)){return(gender)}
  
  #otherwise iterate to the end
  while(i < length(gender)){
    while(!is.na(gender[i])){ #skip years with some gender
      gender.new[i] <- gender[i]
      i <- i+1
    }
    #if this is the end of vector, return it
    if(i > length(gender)){return(gender)}
    
    #if it is not end of vector, store i-1 as left
    left <- i-1
    #get right
    right <- i
    while(is.na(gender[right]) & right <= length(gender)){
      right <- right + 1
    }
    #if there are NAs till the end of vecotor, set left
    if(right > length(gender)){
      gender.new[left:(right-1)] <- gender[left]
      return(gender.new)
    }
    #if the right is not more then length (There is gap in the vector)
    if(right <= length(gender)){
      #get presence between left and right
      first.pres.idx <- which(presence[(left+1):(right-1)] == T)
      if(length(first.pres.idx)==0){
        gender.new[left:(right-1)] <- gender[left]
        gender.new[right] <- gender[right]
      }else{
        first.pres <- left+min(first.pres.idx)
        gender.new[first.pres:(right)] <- gender[right]
        gender.new[left:(first.pres-1)] <- gender[left]
      }
    }
    i <- right
  }
  return(gender.new)
}

fixGenderBetween <- function(gender, priority.genders = c("M", "F")){
  i <- 1
  first.gender <- which(gender %in% priority.genders)
  if(length(first.gender) > 0){
    first.gender <- min(first.gender)
    gender[first.gender:length(gender)] <- gender[first.gender]
  }
  return(gender)
}

#' @inheritParams getPresenceMinTable
#input: vector of timestamps,
#character string indicating what should be the start of the day. Possible values "sunset", "sunrise" or string of time of format "%H:%M:%S"
#roll.forward
getShiftedDate <- function(time.vec, day.start, lat, lon, roll.forward = F, return.timestamp= F){
  if(!day.start %in% c("sunset", "sunrise") & !grepl(day.start, pattern = "[0-9]{2}:[0-9]{2}:[0-9]{2}")){stop("Format of day.start is not valid. Either 'sunset', 'sunrise' or string of format '%H:%M:%S'")}
  if(day.start %in% c("sunset", "sunrise") & (is.na(lat) | is.na(lon))){stop("No lat, lon provided when sunset/sunrise is required.")}
  from <- min(time.vec)-86400
  to <- max(time.vec)+86400
  if(day.start %in% c("sunset", "sunrise")){
    sunrise.sunset.t <- sunrise.set(lat, lon, from , num.days = difftime(time1 = to, time2 = from, units = "days")+2)
    if(roll.forward){
      time.vec.proj <- sunrise.sunset.t[,day.start][findInterval(time.vec, sunrise.sunset.t[,day.start])+1]
    }else{
      time.vec.proj <- sunrise.sunset.t[,day.start][findInterval(time.vec, sunrise.sunset.t[,day.start])]
    }
    
  }else{
    #extract start and end of observation period (add 1 day at each end) - set start of day to 8:00 AM UTC
    from.ss <- as.POSIXlt(strptime(paste(strftime(from - 86400, "%Y-%m-%d"), day.start, sep = " "), '%Y-%m-%d %H:%M:%S'))
    to.ss <- as.POSIXlt(strptime(paste(strftime(to + 86400, "%Y-%m-%d"), day.start, sep = " "), '%Y-%m-%d %H:%M:%S'))
    #get shifted day for data - assign days starting at time specified
    #vector of new midnights
    startdaytime.intervals <- seq(from, to, by = "24 hour")
    #for each line, find newly mapped date
    if(roll.forward){
      time.vec.proj <- startdaytime.intervals[findInterval(time.vec, startdaytime.intervals)+1]
    }else{
      time.vec.proj <- startdaytime.intervals[findInterval(time.vec, startdaytime.intervals)]
    }
    
  }
  if(return.timestamp){return(time.vec.proj)}else{return(as.Date(time.vec.proj))}
}

getNightDay <- function(x, lat = 49.5765639, lon = 14.6637706, xtz = "UTC", returnSunrSuns = F, day = "day", night = "night"){
  if(length(x) == 0) stop("length of input is 0")
  if(all(is.na(x))) return(as.character(x))
  
  library(StreamMetabolism)
  library(data.table)
  from <- min(x, na.rm = T)
  to <- max(x, na.rm = T)
  sunrise.sunset <- as.data.table(sunrise.set(lat, lon, from , num.days = difftime(time1 = to, time2 = from, units = "days", tz = xtz)+2, timezone = xtz))
  sunrise.sunset$Date <- as.Date(sunrise.sunset$sunrise, tz = xtz)
  res <- merge(data.table(x, xorder = 1:length(x), Date = as.Date(x)), sunrise.sunset, by = "Date", all.x = T)
  res[, Diel.period := as.character(ifelse(x > sunrise & x < sunset, "day", "night"))]
  setkey(res, xorder)
  if(returnSunrSuns){
    return(res[, .(Diel.period, sunrise, sunset)])
  }else{
    return(res$Diel.period)  
  }
}

round.POSIXct <- function(x, bin.secs){
  if(class(x)[1] != "POSIXct") stop("parameter x has to be type POSIXct")
  library(lubridate)
  as.POSIXct(round(as.numeric(x)/bin.secs)*bin.secs,origin= "1970-01-01", tz = tz(x))
}



findArrivalsDepartures <- function(x, threshold.s){
  if(any(x != x[order(x)])){print("df is not ordered by time!"); stop()}
  l <- length(x)
  #x <- x[order(x)]
  i <- 1 #index in time vector
  arrival <- x[1]#first is arrival
  arr.dep.l <- list()
  j <- 1
  #iterate over whole vector, if there is a gap bigger than treshold, set the position departure and next position arrival
  while(i < l){
    if(as.numeric(difftime(x[i+1],x[i], units = "secs")) > threshold.s){
      arr.dep.l[[j]] <- data.table(arrival = arrival, departure = x[i])
      arrival <- x[i+1]
      j <- j+1
    }
    i <- i + 1
  }
  arr.dep.l[[j]] <- data.table(arrival = arrival, departure = x[l])
  export.dt <- do.call(rbind, arr.dep.l)
  return(export.dt)
}


#' Get fish presence - each row one minute
#' @param k numeric vector of lenght one (seconds) defining minimal time between consecutive detections to classify it as non-presence
#' @param by_antenna_deployments logical vector of length one defining if the presence should be computed in each antenna separately (TRUE) or on all antennas (FALSE)
#' @param Years - specify years  (numeric vector) from which you want to load fish minutes
#' @param remove_visits_shorter_than numeric vector of length one (seconds) defining threshold on visits to be removed 
#' @param extend_visits_shorter_than numeric vector of length one (seconds) defining which fish visits should be extended 
#' @param extend_visits_by numeric vector of length one (seconds) defining how much he short visits should be extended (visits are extended at the end)
#' @param daystart How to cut the time in order to compute distance of presence from sunset/sunrize. Possible values "sunset", "sunrise" or string of time of format "%H:%M:%S"
#' @param lon longitude of the project place (needed for sunset/sunrise computation)
#' @param lat latitute of the project place (needed for sunset/sunrise computation)
#' @return table of fish minutes with many desriptive information (fishid, time, species, )
getPresenceMinTable <- function(k, Years,  by_antenna_deployments  = F, remove_visits_shorter_than = 0, extend_visits_shorter_than = 15*60,  extend_visits_by = 15*60, daystart = "sunrise", lon =  15.25583, lat = 49.57833)
{
  #get presence table
  if(by_antenna_deployments){
    presence <- data.table(dbGetQuery(con, paste("Select * from asp.presence_per_antenna where k =", k)))
  }else{
    presence <- data.table(dbGetQuery(con, paste("Select * from asp.presenceallk where k =", k)))
  }
  if(nrow(presence) == 0) stop("The result from database is empty. Did you specify proper k?")
  
  presence[, Year := year(timestamp_arr_utc)]
  #subset table by desired years
  presence <- presence[Year %in% Years]
  #captures table
  captures <- data.table(dbGetQuery(con, "Select * from asp.capture"))
  captures[, ":="(ca_timestamp_catch_utc = force_tz(ca_date_catch_utc, "UTC"), ca_timestamp_release_utc = force_tz(ca_date_release_utc, "UTC"))]
  
  #gender year table
  gender_year <- data.table(dbGetQuery(con, "Select * from asp.gender_year"))
  fish <-  data.table(dbGetQuery(con, "Select * from asp.fish"))
  gender_year[, Year := year]
  gender_year[, fi_fishid := fish_id]
  
  gender_year <- merge( fish[, .(fi_fishid, fi_species)], gender_year[, Year := year], by = "fi_fishid", all.x= T)#TODO
  #load capture-tag table
  capture_tag <- data.table(dbGetQuery(con, "Select * from asp.capture_tag natural join asp.tag_table"))
  #detections TODO: these are needed?
  # detections <- data.table(dbGetQuery(con, "Select * from asp.detection_fish_tagging"))
  # detections[, dt_timestamp_utc := force_tz(dt_timestamp_utc, "UTC")]
  
  #get length of stay for each fish
  presence[, visitlength := as.numeric(difftime(timestamp_dep_utc, timestamp_arr_utc, "secs"))]
  
  #get overview of visit lengths
  #timespent.y.id <- presence[, .(timespent.y = sum(visitlength)), by = .(fi_fishid = fi_fishid, Year)]
  
  #update visits if the length of stay is too short
  presence <- presence[visitlength > remove_visits_shorter_than]
  presence <- presence[visitlength < extend_visits_shorter_than, timestamp_dep_utc := timestamp_dep_utc + extend_visits_by]
  
  # save.image("~/asp.Rdata")
  #load("~/asp.Rdata")
  #for each visit, get time from last tagging  add T/F if the fish was tagged in the same season
  #get captures which were with taggings (not recaptures)
  cap.taggings <- merge(captures[month(ca_timestamp_catch_utc) %in% c(3,4)], capture_tag[ct_tagging == T], by.x= "ca_captureid", by.y = "ca_captureid")
  first_antenna_tagging_year <- cap.taggings[tt_type == "antenna", .(first_antenna_tagging_year = year(min(ca_date_release_utc))), by = .(fi_fishid)]
  
  #roll join to last tagging of fish - use start of visit and release time for join
  cap.taggings[, tm := ca_timestamp_release_utc]
  presence[, tm := timestamp_arr_utc]
  setkey(presence, fi_fishid, tm)
  setkey(cap.taggings, fi_fishid, tm)
  presence <- cap.taggings[,.(tm, fi_fishid, ca_timestamp_release_utc)][presence, , roll = Inf]
  tagging.y <- unique(cap.taggings[,.(fi_fishid, tagging.year = year(ca_timestamp_catch_utc))])
  tagging.y$tagging.yl <- T
  #tagging.y[, ntaggigs := .N, fi_fishid]
  presence[, Year := year(timestamp_arr_utc)]
  presence <- merge(presence, tagging.y, by.x =c("fi_fishid","Year"), by.y = c("fi_fishid", "tagging.year"), all.x = T)
  presence[is.na(tagging.yl), tagging.yl := F]
  
  #for each visit, add proper gender (not a gender from capture but from gender_year table (special recomputed table with new gender))
  presence <- merge(presence, gender_year[, .(fi_fishid, Year, gender, fi_species)], by = c("fi_fishid", "Year"), all.x = T)
  presence <- merge(presence, first_antenna_tagging_year, by = c("fi_fishid"), all.x = T)
  # Subset data ####
  presence.sub <- copy(presence)
  #remove all visits shorter than 1 s (the fish was detected only once in 8 hour window!)
  #presence.sub <- presence.sub[visitlength > 1]
  # presence.sub[visitlength < 5 * 60, timestamp_dep_utc := timestamp_arr_utc +20*60*60]
  #remove asp visits of asp which were tagged that spawning season (leave other fi_species)
  #presence.sub  <- presence.sub[!(tagging.yl == T & fi_species == "asp")]
  #remove bleaks visits which started less then 1 day after tagging
  #presence.sub <- presence.sub[!(fi_species == "bleak" & timestamp_arr_utc-86400 < ca_timestamp_release_utc)]
  #create one line for each fish minute
  #vector defining names of columns creating groups for which the minites are expanded
  minute_columns <- c("fi_fishid", "Year", "visitid", "visitlength", "tagging.yl", "gender"," fi_species", "first_antenna_tagging_year")
  if(by_antenna_deployments) minute_columns <- c(minute_columns, "ad_antdep_id")
  presence.min <- presence.sub[, .(presencem = as.POSIXct(seq.POSIXt(from = timestamp_arr_utc, to = timestamp_dep_utc, by = "min"))) , by = .(minute_columns)]
  
  #FISH MINUTES####
  #FISH MINUTES-adding shifted day + sunsetsunrise####
  #add date which the minute false into (shifted days starting at 10:00 am UTC)
  #get time limits (first and last detection)
  
  #add sunset sunrise - each shifted day will have sunset and subsequent sunrise assigned TODO: what happens if the new midnight will be in between sunset and subsequent sunrise?
  #get sunrise/sunset times for whole observation period (+1 day)
  
  from <- min(presence.min$presencem)-84600
  to <- max(presence.min$presencem)+86400
  #get times of sunset and sunrise and get shifted date for both 
  sunrise.sunset <- data.table(sunrise.set(lat = lat, long = lon, from , num.days = difftime(time1 = to, time2 = from, units = "days")+2, timezone = "UTC"))
  sunrise.sunset[, sunset.shiftedDate := getShiftedDate(time.vec = sunset, lat = lat, lon = lon, day.start = daystart, return.timestamp= F)]
  sunrise.sunset[, sunrise.shiftedDate := getShiftedDate(time.vec = sunrise, lat = lat, lon = lon, day.start = daystart, return.timestamp= F)]
  #get shifted date for presence
  presence.min[, shifteddate := getShiftedDate(time.vec = presencem, lat = lat, lon = lon, day.start = daystart, return.timestamp= F)]
  #add proper sunset to data
  presence.min <- merge(presence.min, sunrise.sunset[,.(sunrise, shifteddate = sunrise.shiftedDate)], by = "shifteddate", all.x = T)
  presence.min <- merge(presence.min, sunrise.sunset[,.(sunset, shifteddate = sunset.shiftedDate)], by = "shifteddate", all.x = T)
  return(presence.min)
}



readGPX_tracks <- function(file){
  #get track points and standardize them
  tracks.gpx.i <-
    readGPX(
      gpx.file = file,
      metadata = F,
      tracks = T
    )
  tracks2 <-
    lapply(
      tracks.gpx.i$tracks,
      FUN = function(x)
        do.call("rbind", x)
    )
  tracks3 <-
    lapply(tracks2, function(x)
      if ("time" %in% names(x))
        x
      else
        NULL)
  #rbind all points from gpx files
  tracks <- data.table(do.call("rbind", tracks3))
  #fix time in order to be POSIXct format in "UTC" timezone
  tracks[, time := as.POSIXct(x = time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")]
  return(tracks)
}

#function to cross join two data.tables
CJ.tables <- function(X,Y){ 
  unique_name <- last(make.unique(c(colnames(X),colnames(Y),"temporalColumnForCJtables"))) 
  res <- X[,c(setNames(1,unique_name),.SD)][Y[,c(setNames(1,unique_name),.SD)],on=unique_name,allow.cartesian=TRUE]
  res$temporalColumnForCJtables <-  NULL
  return(res) 
}
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

glm_nb_island <- function(dataset, formula){
  glm_nb <- glm.nb(formula = formula, data = dataset)
  r2_kl <- round(rsq(fitObj = glm_nb, adj = F, type = c("kl")), 2)
  r2_zhang <- round(rsq(glm_nb),2)
  pseudo <- modEvA::RsqGLM(glm_nb)
  r2_CoxSnell <- round(pseudo$CoxSnell,2)
  r2_Nagelkerke <- round(pseudo$Nagelkerke,2)
  r2_McFadden <- round(pseudo$McFadden,2)
  r2_sqPearson <- round(pseudo$sqPearson,2)
  glm_aic <- round(AIC(glm_nb), 2)
  mnull <- update(glm_nb, . ~ 1)#Creating a null model
  fit.test <- round(pchisq(2 * (logLik(glm_nb) - logLik(mnull)), df = 1, lower.tail = FALSE), 2)#Comparing the null model with the full model. Attention to the df, it varies according to the number of predictors used
  return(list(Null_mdl_comparison = fit.test,
              r2_zhang = r2_zhang,
              r2_kl_divergence_based = r2_kl,
              r2_CoxSnell = r2_CoxSnell,
              r2_Nagelkerke = r2_Nagelkerke,
              r2_McFadden = r2_McFadden,
              r2_sqPearson = r2_sqPearson,
              aic = glm_aic,
              intercept = round(coef(glm_nb)[1],2),
              slope_Treatment = round(coef(glm_nb)[2],2),
              p_Treatment = round(coef(summary(glm_nb))[2,4],2),
              significance_Treatment = star(round(coef(summary(glm_nb))[2,4],2))))
}

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

#combine graphs with shared legend
grid_arrange_shared_legend <- 
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }

wordCounter <- function(invec, word, ...) {
  vapply(regmatches(invec, gregexpr(word, invec, ...)), length, 1L)
}

#Function for legal size####
status <- list() 
legal.size <- function(data, latin.name, total.length.mm) {
  for(i in 1:nrow(data)){
    if (latin.name[i] == 'Perca fluviatilis' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#Perch
    else if (latin.name[i] == 'Perca fluviatilis' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#Perch
    else if (latin.name[i] == 'Abramis brama' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#Bream
    else if (latin.name[i] == 'Abramis brama' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#Bream
    else if (latin.name[i] == 'Alburnus alburnus' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#Bleak
    else if (latin.name[i] == 'Alburnus alburnus' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#Bleak
    else if (latin.name[i] == 'Anguilla anguilla' & total.length.mm[i] >= 500) {status[i] <- 'legal'}#Eel
    else if (latin.name[i] == 'Anguilla anguilla' & total.length.mm[i] < 500) {status[i] <- 'ilegal'}#Eel
    else if (latin.name[i] == 'Blicca bjoerkna' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#White bream
    else if (latin.name[i] == 'Blicca bjoerkna' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#White bream
    else if (latin.name[i] == 'Carassius gibelio' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#Prussian carp
    else if (latin.name[i] == 'Carassius gibelio' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#Prussian carp
    else if (latin.name[i] == 'Cyprinus carpio' & total.length.mm[i] >= 400) {status[i] <- 'legal'}#Carp
    else if (latin.name[i] == 'Cyprinus carpio' & total.length.mm[i] < 400) {status[i] <- 'ilegal'}#Carp
    else if (latin.name[i] == 'Esox lucius' & total.length.mm[i] >= 500) {status[i] <- 'legal'}#Pike
    else if (latin.name[i] == 'Esox lucius' & total.length.mm[i] < 500) {status[i] <- 'ilegal'}#Pike
    else if (latin.name[i] == 'Gobio gobio' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#Gudgeon
    else if (latin.name[i] == 'Gobio gobio' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#Gudgeon
    else if (latin.name[i] == 'Gymnocephalus cernua' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#Ruffe
    else if (latin.name[i] == 'Gymnocephalus cernua' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#Ruffe
    else if (latin.name[i] == 'Leuciscus aspius' & total.length.mm[i] >= 400) {status[i] <- 'legal'}#Asp
    else if (latin.name[i] == 'Leuciscus aspius' & total.length.mm[i] < 400) {status[i] <- 'ilegal'}#Asp
    else if (latin.name[i] == 'Oncorhynchus mykiss' & total.length.mm[i] >= 250) {status[i] <- 'legal'}#Trout
    else if (latin.name[i] == 'Oncorhynchus mykiss' & total.length.mm[i] < 250) {status[i] <- 'ilegal'}#Trout
    else if (latin.name[i] == 'Rutilus rutilus' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#Roach
    else if (latin.name[i] == 'Rutilus rutilus' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#Roach
    else if (latin.name[i] == 'Scardinius erythrophthalmus' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#Rudd
    else if (latin.name[i] == 'Scardinius erythrophthalmus' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#Rudd
    else if (latin.name[i] == 'Silurus glanis' & total.length.mm[i] >= 700) {status[i] <- 'legal'}#Wels
    else if (latin.name[i] == 'Silurus glanis' & total.length.mm[i] < 700) {status[i] <- 'ilegal'}#Wels
    else if (latin.name[i] == 'Tinca tinca' & total.length.mm[i] >= 200) {status[i] <- 'legal'}#Tench
    else if (latin.name[i] == 'Tinca tinca' & total.length.mm[i] < 200) {status[i] <- 'ilegal'}#Tench
    else if (latin.name[i] == 'Sander lucioperca' & total.length.mm[i] >= 450) {status[i] <- 'legal'}#Pikeperch
    else if (latin.name[i] == 'Sander lucioperca' & total.length.mm[i] < 450) {status[i] <- 'ilegal'}#Pikeperch 
    else if (latin.name[i] == 'Leuciscus idus' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#Idus
    else if (latin.name[i] == 'Leuciscus idus' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#Idus
    else if (latin.name[i] == 'Hypophthalmichthys molitrix' & total.length.mm[i] >= 150) {status[i] <- 'legal'}#Whitefish
    else if (latin.name[i] == 'Hypophthalmichthys molitrix' & total.length.mm[i] < 150) {status[i] <- 'ilegal'}#Whitefish
    else {status[i] <- 'unknown'}
    status.l <- data.table(legal_size = do.call(rbind, status))
    setnames(x = status.l, old = 'legal_size.V1', new = 'Status')
  }
  return(status.l)
}

