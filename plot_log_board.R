while(TRUE) {
  #time_flag <- as.numeric(strsplit(as.character(Sys.time()), "\\:")[[1]][3]) %% 10 == 0
  time_flag <- TRUE
  axis_flag <- FALSE
  ylim <- 3
  
  if(time_flag) {
    df <- read.table("log_file.txt", sep="|", row.names=NULL, header=TRUE)
    y <- df$result_val

    if(axis_flag) {
      plot(1:length(y), y, pch=20, cex=0.25)      
    } else {
      plot(1:length(y), y, pch=20, cex=0.25, ylim=c(0, ylim))
    }    
  }

  Sys.sleep(5)
}

size <- 3
val <- 1
1e2*((1/size)-val)

size <- 15
val <- 1
1e2*((1/size)-val)

size <- 80
val <- 1
1e2*((1/size)-val)


