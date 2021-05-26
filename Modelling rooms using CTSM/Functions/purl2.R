purl2 <- function(nm){
  ## First remove the old .R file
  file.remove(paste0(nm,".R"))
  ## Then make the new file
  purl(paste0(nm,".Rnw"), documentation=1)
  ## Read
  x <- scan(paste0(nm,".R"), what="character", sep="\n",
            blank.lines.skip=FALSE)
  ## Find chunk header lines
  ilines <- grep("## ----", x)
  ilines[length(ilines)+1] <- length(x)+1
  ## Get the content
  for(i in 1:(length(ilines)-1)){
    ## Chunck with eval=F?
    if(length(grep("eval=F", gsub(" ","",x[ilines[i]]))) > 0){
      itmp <- (ilines[i]+1):(ilines[i+1]-1)
      x[itmp] <- substr(x[itmp], 4, nchar(x[itmp]))
    }
  }
  ## Remove the header lines
  x <- x[-ilines[-length(ilines)]]
  ####
  ## Remove empty lines in the top and bottom
  tmp <- cumsum(gsub(" ","",x) != "")
  if(which(tmp == 1) > 1){
    ## Remove empty lines in the top
    x <- x[-(1:(which(tmp == 1)-1))]
  }
  ## 
  tmp <- cumsum(gsub(" ","",rev(x)) != "")
  if(length(which(tmp == 1)) > 0){
    ## Remove empty lines in the bottom
    x <- x[-(length(x)-((max(which(tmp == 1))-2):0))]
  }
  ####
  ## Write it to the .R file
  write(x, file=paste0(nm,".R"))
  ## Return if needed
  invisible(x)
}
