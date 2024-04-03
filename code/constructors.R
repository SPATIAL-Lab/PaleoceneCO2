
tsvbl = function(iv = c(0, 1), exp_range = iv, er_dist = "g", 
                      bounds = NULL){
  tsvbl = list()
  
  if(!inherits(iv, "numeric") | !(length(iv) %in% c(1, 2))){
    stop("iv must be numeric of length 1 or 2")
  }
  tsvbl$iv = iv
  
  if(!inherits(exp_range, "numeric") | length(exp_range) != 2){
    stop("exp_range must be numeric of length 2")
  }
  tsvbl$exp_range = exp_range
  
  if(!inherits(er_dist, "character") | !(er_dist %in% c("g", "c"))){
    stop("er_dist must be one of g, c")
  }
  tsvbl$er_dist = er_dist
  
  if(!is.null(bounds)){
    if(!inherits(bounds, "numeric") | length(exp_range) != 2){
      stop("bounds must be numeric length 2")
    }
    tsvbl$bounds = bounds
  }
  
  class(tsvbl) = "tsvbl"
  return(tsvbl)
}

vbl = function(parms = c(0, 1), dist = "g", bounds = NULL, 
                    fixed = FALSE){
  vbl = list()
  
  if(!inherits(parms, "numeric") | length(parms) != 2){
    stop("parms must be numeric of length 2")
  }
  vbl$parms = parms
  
  if(!inherits(dist, "character") | !(dist %in% c("g", "u"))){
    stop("dist must be one of g, u")
  }
  vbl$dist = dist
  
  if(!is.null(bounds)){
    if(!inherits(bounds, "numeric") | length(exp_range) != 2){
      stop("bounds must be numeric length 2")
    }
  }
  vbl$bounds = bounds
  
  if(!(inherits(fixed, "logical"))){
    stop("fixed must be TRUE/FALSE")
  }
  vbl$fixed = fixed
  
  class(vbl) = "vbl"
  return(vbl)
}

link = function(depvbl, expr){
 link = list()
 
 if(!inherits(depvbl, "character") | !inherits(expr, "character")){
   stop("depvbl and expr must both be character strings")
 }
 link$depvbl = depvbl
 link$expr = expr
 
 class(link) = "link"
 return(link)
}

ts = function(...){
  n = ...length()
  
  #check input
  for(i in seq_len(n)){
    if((!is.numeric(...elt(i)))){
      stop("all arguments must be numeric")
    }
  }
  
  #merge elements
  ts = numeric()
  for(i in seq_len(n)){
    ts = append(ts, ...elt(i))
  }
  
  #sorted and unique
  ts = sort(unique(ts))
  
  #create index vectors
  ts_ind = list()
  for(i in seq_len(n)){
    ts_ind[[i]] = match(...elt(i), ts)
  }
  names(ts_ind) = ...names()
  
  ts = list("ts" = ts, "ts_ind" = ts_ind)
  class(ts) = "ts"
  
  return(ts)
}

dat = function(vbl_name, values, sd, ts_ind){
  if(!inherits(vbl_name, "character")){
    stop("vbl_name must be character")
  }
  
  dat = list()
  
  dat$data = values
  dat$pre = 1 / sd^2
  dat$ts_ind = ts_ind
  
  class(dat) = dat
  return(dat)
}

build = function(vbls, links, tser, dat){
  
  #model output file
  fc = tempfile(fileext = ".txt")
  f = file(fc, "w")
  
  #header
  cat("model{", file = f, sep = "\n")
  
  #parse variables list
  v = vbls[unlist(lapply(vbls, "inherits", "vbl"))]
  tsv = vbls[unlist(lapply(vbls, "inherits", "tsvbl"))]
  #this will hold variables which are not fixed in time
  vv = list()
  
  #write model code for time-independent variables
  for(i in seq_along(v)){
    if(v[[i]]$fixed){
      mline = names(v)[i]
      mline = paste(mline, "~")
      if(v[[i]]$dist == "g"){
        mline = paste0(mline, " dnorm(", v[[i]]$parms[1], ", ", 
                       1 / v[[i]]$parms[2]^2, ")")
        if(!is.null(v[[i]]$bounds)){
          mline = paste0(mline, "T(", v[[i]]$bounds[1], ", ", 
                         v[[i]]$bounds[2], ")")
        }
      } else if(v[[i]]$dist == "u"){
        mline = paste0(mline, " dunif(", v[[i]]$parms[1], ", ",
                       v[[i]]$parms[2], ")")
      }
      cat(mline, file = f, sep = "\n")
    } else{
      #stash the time-varying variables in vv
      vv = append(vv, v[i])
    }
  }
  
  #now deal with time-dependent
  if(length(vv) > 0){
    #the loop goes through all timesteps
    cat("for(i in 1:length(tser)){", file = f, sep = "\n")
    #no temporal correlation
    for(i in seq_along(vv)){
      mline = names(vv)[i]
      mline = paste0(mline, "[i] ~ ")
      #gaussian
      if(vv[[i]]$dist == "g"){
        mline = paste0(mline, "dnorm(", vv[[i]]$parms[1], ", ", 
            1 / vv[[i]]$parms[2]^2, ")")
        if(!is.null(vv[[i]]$bounds)){
          mline = paste0(mline, "T(", vv[[i]]$bounds[1], ", ", 
                         vv[[i]]$bounds[2], ")")
        }
      }
      #uniform
      if(vv[[i]]$dist == "u"){
        mline = paste0(mline, "dunif(", vv[[i]]$parms[1], ", ", 
                       vv[[i]]$parms[2], ")")
      }
      cat(mline, file = f, sep = "\n")
    }
    cat("}", file = f, sep = "\n")
  }
  
  if(length(tsv) > 0){
    #initial conditions
    for(i in seq_along(tsv)){
      mline = names(tsv)[i]
      if(length(tsv[[i]]$iv == 1)){
        mline = paste0(mline, "[1] = ", tsv[[i]]$iv)
      } else{
        mline = paste0(mline, "[1] ~ dunif(", tsv[[i]]$iv[1], ", ", 
                       tsv[[i]]$iv[2], ")")
      }
      cat(mline, file = f, sep = "\n")
    }
    #timeseries loop
    cat("for i in 2:length(tser){", file = f, sep = "\n")
    for(i in seq_along(tsv)){
      mline = paste0(names(tsv)[i], "[i] ~ ")
      #gaussian
      if(tsv[[i]]$er_dist == "g"){
        mline = paste0(mline, "dnorm(", names(tsv)[i], "[i-1], ", 
                       1 / (diff(tsv[[i]]$exp_range) / 3)^2, ")")
        if(!is.null(tsv[[i]]$bounds)){
          mline = paste0(mline, "T(", tsv[[i]]$bounds[1], ", ", 
                         tsv[[i]]$bounds[2], ")")
        }
      }
      #cauchey
      if(tsv[[i]]$er_dist == "c"){
        mline = paste0(mline, "dt(", names(tsv)[i], "[i-1], ", 
                       1 / (diff(tsv[[i]]$exp_range) / 3)^2, ", 1)")
        if(!is.null(tsv[[i]]$bounds)){
          mline = paste0(mline, "T(", tsv[[i]]$bounds[1], ", ", 
                         tsv[[i]]$bounds[2], ")")
        }
      }
      cat(mline, file = f, sep = "\n")
    }
    cat("}", file = f, sep = "\n")
  }
  
  #links
  #need to split those w/ time depednece and those without
  
  
  #data
}
