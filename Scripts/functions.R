process_data <- function(trait_filename = "katie_data.csv")
{
  dat <- read.csv(paste("data/",trait_filename,sep=""),stringsAsFactors = FALSE)
  colnames(dat)[1] <- "SAM"
  colnames(dat) <- gsub(pattern = "_",replacement = "",x = colnames(dat))
  traits <- colnames(dat)[-1]
  envs <- c("Dry","logdiff","Wet")
  
  if(is.numeric(dat$SAM))
  {
    dat$SAM <- paste("SAM",sapply(3-nchar(dat$SAM),function(X) paste("",rep("0",X),sep="",collapse="")),dat$SAM,sep="")
  }
  if(!any(grepl("SAM001",dat$SAM))) stop("dat must be formatted as integers (SAM lines) or of the format e.g. SAM001, SAM002, etc.")
  load("data/lines.RDat")
  
  dat <- dat[match(x = lines,table = dat$SAM),]
  dat$SAM <- lines
  
  writeLines(text = colnames(dat)[-1],"traits_to_run.txt")
  writeLines(text = envs,"environments_to_run.txt")
  
  new_names <- c("SAM",kronecker(X = colnames(dat)[-1],Y = envs,FUN = function(X,Y) paste(X,Y,sep="_")))
  dat <- cbind(dat,dat[,-1] + rnorm(n = length(unlist(dat[,-1])),sd = 100),dat[,-1] + rnorm(n = length(unlist(dat[,-1])),sd = 100))
  colnames(dat) <- new_names
  
  write.csv(x = dat,file = paste("data/",gsub(pattern = ".csv",replacement = "",x = trait_filename),"_for_pipeline.csv",sep=""),row.names = FALSE)
  
  prefs <- readLines("Scripts/original### Preferences ###")
  phen <-grep("replace1",prefs)
  prefs[phen] <- gsub(pattern = "replace1",replacement = paste(gsub(pattern = ".csv",replacement = "",trait_filename),"_for_pipeline.csv",sep=""),x = prefs[phen])
  writeLines(text = prefs,con = "Scripts/### Preferences ###")
}

set_threshold <- function(method)
{
  if(method==1)
  {
    ntests <- 20562
  } else if(method==2)
  {
    ntests <- 50
  } else
  {
    if(is.na(as.numeric(method)))
    {
      stop("******************************************************************************************************
       Please set the argument method to one of the following values:\n\t1) for standard Bonferroni; assumes ntests = 20562  (p-value threshold = .05/ntests)
        2) Sets threshold to 0.001 by artificially assuming 50 tests.
        
        Alternatively, you can enter a custom number of tests in place of method.
       ***********************************************************************************************************")
    }
  }
  prefs <- readLines("Scripts/original### Preferences ###")
  thresh <-grep("replace2",prefs)
  prefs[thresh] <- gsub(pattern = "replace2",replacement = ntests,x = prefs[thresh])
  writeLines(text = prefs,con = "Scripts/### Preferences ###")
}

