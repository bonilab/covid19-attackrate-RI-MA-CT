
traj.from.params <- function(beta,
                             params = NULL,
                             const.params = NULL,
                             non.odesim.params=NULL,
                             introday = NULL,
                             tf,
                             odepath,
                             loc="RI",
                             symp = NULL,
                             mean.report.rate = 1.0
){
  # A function that takes beta values as input and calls the odesim program.
  #   Input:
  #       1) beta: vector of beta values used in the sim
  #               should be a vector of length tf-60, with one value per day
  #       2) tf: integer day to stop sim (1 = 20200101)
  #       3) t0=jan 1, 2020
  #   Output:
  #       Data frame with times column and DeltaJ (simulated cases) column
  
  
  ## remove non.odesim.params from params and const.params
  if(length(non.odesim.params)>0){
    idx.remove <- integer()
    
    for(k in 1:length(non.odesim.params)){
      idx.remove <- c(idx.remove,which(names(params)==non.odesim.params[k]))
    }
    
    if (length(idx.remove) > 0){
      params <- params[-idx.remove]
    } else {
      params <- NULL
    }
    
    # if(length(const.params)>0){
    #   
    #   idx.remove <- integer()
    #   for(k in 1:length(non.odesim.params)){
    #     idx.remove <- c(idx.remove, which(names(const.params) == non.odesim.params[k]))
    #   }
    #   
    #   if (length(idx.remove) > 0){
    #     const.params <- const.params[-idx.remove]
    #   } else {
    #     const.params <- NULL
    #   }
    # }
  }
  
  
  ### make introday backwards compatible ###
  
  if(length(which(c(names(params), names(const.params)) == "introday")) + length(introday) == 0){
    # if no introday, set to 55 by default
    introday = 55
  } else if (length(which(c(names(params), names(const.params)) == "introday")) + length(introday) > 1){
    # print warning (introday specified twice)
    cat("WARNING: introday specified in two or more places.","\n")
  }
  
  ### create string to run odesim with parameters ###
  cmd1 <- paste(c(odepath,"odesim none -binary-output -tf ",as.character(tf)),collapse="")
  if(length(introday)==1){
    cmd1.5 <- paste(c(" -introday"), as.character(introday), collapse = " ")
    cmd2 <- paste(c(" -beta", as.character(round(beta,9))), collapse = " ")
    cmd <- paste(cmd1, cmd1.5, cmd2, sep="")
  } else {
    cmd2 <- paste(c(" -beta", as.character(round(beta,9))), collapse = " ")
    cmd <- paste(cmd1, cmd2, sep="")
  }
  
  ### add constants ###
  
  # initialize to 0-length vector
  tv_idx <- integer() 
  tv_params <- c()
  
  # if constants exist, add to 'cmd'
  if(length(const.params) > 0){
    
    tv_idx <- grep("^tv-.*_[0-9]+$", names(const.params)) ## get indices of all tv- params
    if(length(tv_idx) > 0){
      tv_params <- const.params[tv_idx]
      # sort by name (ascending)
      tv_params <- tv_params[sort(names(tv_params))] 
      # params that are not tv-
      const.params <- const.params[ -tv_idx ]
    }
    
    if(length(const.params) > 0){
      cmd <- paste(cmd, paste(" -", names(const.params), " ", 
                              as.character(const.params), 
                              sep = "", collapse = ""), 
                   sep = "")
    }
  }
  
  ### add inferred parameters ###
  if (length(params) > 0){
    
    # parameter names
    names <- names(params)
    
    # separating time-varying params (anything starting with "tv-") in params vector
    tv_idx <- grep("^tv-", names) ## get indices of all tv- params
    # tv_params <- c()
    if(length(tv_idx) > 0){
      tv_params <- c(tv_params, params[tv_idx])
      
      # if (length(tv_idx_const) > 0){
      #   tv_params <- c(tv_params, const.params[tv_idx_const])
      #   tv.names <- names(tv_params)
      #   tv_params <- as.numeric(tv_params)
      #   names(tv_params) <- tv.names
      # }
      
      # sort by name (ascending)
      tv_params <- tv_params[sort(names(tv_params))] 
      # params that are not tv-
      params <- params[ -tv_idx ]
    }
    
    # append non-time-varying params to the command
    if (length(params) > 0){
      cmd <- paste(cmd, paste(" -", names(params), " ", 
                              as.character(round(params, 9)), 
                              sep = "", collapse = ""), 
                   sep = "")
    }
  }
  
  # process time-varying params 
  if(length(tv_params) > 0){
    # print(tv_params)
    
    tv_names <- names(tv_params)
    tv_nm_spl <- strsplit(tv_names, "_", fixed = TRUE) ## split names at "_"
    
    k <- 1
    while(k <= length(tv_params)){
      if( is.na(tv_nm_spl[[k]][2]) ){
        nm_prefix <- paste0("^",tv_nm_spl[[k]][1]) # tv_nm_spl[[k]][1] # paste("^",tv_nm_spl[[k]][1], sep="")
      } else{
        nm_prefix <- paste0("^",tv_nm_spl[[k]][1],"_") # paste("^",tv_nm_spl[[k]][1],"_", sep="")
      }
      # print(nm_prefix)
      sub_idx <- grep(nm_prefix, tv_names) ## find elements in tv_names matching the current tv param name
      
      cmd <- paste(c(cmd," -",tv_nm_spl[[k]][1]), collapse = "") ## append CLO name
      cmd <- paste(c(cmd, as.character( round(as.numeric(tv_params[sub_idx]),9) )), collapse = " ") ## append CLO values
      k <- k + length(sub_idx) 
    }
  }
  
  if (length(symp) > 0){
    # specify symp-frac (used when sf.choice == TRUE)
    cmd <- paste(cmd, symp, collapse = "")
  }
  
  # add mean reporting rate (for vaccination)
  cmd <- paste(cmd, "-rr", as.character(round(mean.report.rate,5)), collapse = "")
  # add location parameter
  cmd <- paste(cmd, "-loc", as.character(loc), "2>&1", collapse="")
  
  # cat(cmd, "\n")
  
  
  ncols  = 307 #Number of expected columns in output
  nlines = as.integer(700) #Upper limit on number of output lines
  #nlines = as.integer(tf - introday + 1) # Output starts on introday
  #nlines = as.integer(tf) # Output starts on 0
  
  
  # call the odesim program, with beta as input
  p = pipe(cmd, 'rb')
  bin_result = readBin(p, 'double', nlines*ncols)
  close.connection(p)
  
  ##  browser()
  # cat(bin_result, "\n")
  
  #Reshaping 1D result to (*-by-ncols)
  dim(bin_result) = c(ncols, length(bin_result) / ncols) # Because dim fills row-first :(
  bin_result = aperm(bin_result, c(2,1))
  
  bin_result
}
