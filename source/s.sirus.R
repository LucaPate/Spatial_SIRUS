# The following code is an adaptation of the SIRUS functions available at https://github.com/cran/sirus/tree/master/R .
# This code can be used for the implementation of Spatial SIRUS (S-SIRUS).
# At the moment the S-SIRUS functions can be used ONLY for regression applications.
# The use of the main S-SIRUS functions is described in the companion paper 
# and a guided example is presented in the s.sirus_guided_example.R script file.


##### Fit the Spatial SIRUS rules
# Differently from the original function in the SIRUS package, s.sirus.fit requires:
# data = matrix of the features, the first two must be the coordinates
# y = vector of the response values
# incl.coords = logic to consider coordinates provided in data as predictors (TRUE) or not (FALSE)
# type = type of application, 'reg' for regression
# max.depth = must be equal to 2

s.sirus.fit <- function(data, y, incl.coords = NULL, type = 'reg', num.rule = 10, p0 = NULL, num.rule.max = 25, q = 10, 
                        discrete.limit = 10, num.trees.step = 1000, alpha = 0.05, mtry = NULL,
                        max.depth = 2, num.trees = NULL, num.threads = NULL, replace = TRUE,
                        sample.fraction = ifelse(replace, 1, 0.632), verbose = TRUE, seed = NULL) {
  
  # Check arguments
  # check data type
  data.y.check(data, y) 
  # check S-SIRUS parameters
  s.sirus.param.check(data, num.rule.max, q, num.trees.step, alpha, mtry)
  
  # check num.rule
  num.rule.valid <- is.numeric(num.rule)
  if (num.rule.valid){num.rule.valid <- (round(num.rule) == num.rule) & num.rule > 0}
  if (!is.null(num.rule)){
    if (!num.rule.valid){
      stop("Invalid num.rule. Number of rules must be a positive integer or NULL.")
    }else{
      if (num.rule > 100){
        warning('Warning num.rule: S-SIRUS is designed to output short list of rules (typically < 100 rules).')
      }
    }
  }
  
  # check p0
  if (!is.null(p0)){
    p0.valid <- is.numeric(p0) & (p0 >= 0) & (p0 <= 1)
    if (!p0.valid){
      stop("Invalid p0. p0 must be a numeric value between 0 and 1 or NULL.")
    }
  }
  # check either p0 or num.rule not null
  if (is.null(p0) & is.null(num.rule)){
    stop('Invalid p0 and num.rule: Either p0 or num.rule has to be provided.')
  }
  
  # set default mtry
  if (is.null(mtry) & incl.coords == TRUE){
    mtry.ratio <- 1/3
    mtry <- max(1, floor((ncol(data))*mtry.ratio))
  }
  if (is.null(mtry) & incl.coords == FALSE){
    mtry.ratio <- 1/3
    mtry <- max(1, floor((ncol(data)-2)*mtry.ratio))
  }
  
  # set type
  type.valid <- type %in% c('reg')
  if (type.valid){
    if (type == 'reg'){
      if (length(unique(y)) == 2){
        warning('Warning type: y takes only two distinct values, classification cannot be implemented.')
      }
    }
  }else{
    stop('Invalid type. Type should reg (for regression).')
  }
  
  # Data info
  if(incl.coords == TRUE){
    data.names <- colnames(data)}    
  
  if(incl.coords == FALSE){
    data.names <- colnames(data[,-c(1:2)])}
  
  mean.out <- mean(y)
  # Coords
  coords <- data[,1:2]
  
  # Data binning
  bins.list <- lapply(data, get.bins, y = y, q = q, discrete.limit = discrete.limit)
  num.cat.valid <- sapply(bins.list, function(bins){
    if (bins$type =='categorical'){length(bins$levels) < 0.05*nrow(data)}else{TRUE}
  })
  if (!all(num.cat.valid)){
    names.warning <- paste0(data.names[!num.cat.valid], collapse = ', ')
    warning(paste0('Some categorical variables have a high number of categories. 
    Spatial SIRUS is likely to identify irrelevant rules, overfit, and to have long running times.
    It is recommended to discard or transform these variables: ', names.warning, '.'))
  }
  data.bin <- as.data.frame(sapply(1:ncol(data), function(j){
    binarize.X(X = data[,j], bins = bins.list[[j]], q = q)}))
  data.bin.y <- cbind(data.bin, y)
  
  if (nrow(data) > 5){
    # Grow forest
    forest <- rfgls.stab(data.bin.y, num.trees.step, alpha, mtry, max.depth, num.trees,
                           num.threads, replace, sample.fraction, verbose, seed, coords = coords, 
                           incl.coords)
    # forest <- result
    paths <- forest$paths[-1]
    proba <- forest$proba[-1]
    num.trees <- forest$num.trees
  }else{
    warning('Minimum sample size to run Spatial SIRUS is 5.')
    paths <- list()
    proba <- 1
    num.trees <- 0
  }
  
  if (length(paths) > 0){
    
    # path selection with p0
    if (!is.null(p0)) {
      selector <- proba > p0
      paths <- paths[selector]
      proba <- proba[selector]
      num.rule <- num.rule.max
    }
    
    # path post-treatment
    paths <- lapply(paths, function(path){
      lapply(path, function(split){
        if (split[2] == round(split[2])){
          split[2] <- split[2] - 0.5
        }
        split
      })
    })
    if (max.depth <= 2){
      paths.ftr <- paths.filter.2(paths, proba, num.rule)
    }else{
      paths.ftr <- paths.filter.d(paths, proba, num.rule, data.bin, incl.coords)
    }
    paths <- paths.ftr$paths
    proba <- paths.ftr$proba
    
    if (incl.coords == FALSE){
      bins.list = bins.list[-c(1:2)]
      data = data[,-c(1:2)]
    }
    
    # format paths
    paths <- lapply(paths, format.path, bins.list = bins.list)
    
    # build rules from paths
    rules <- lapply(paths, get.rule, bins.list = bins.list, data.names = data.names)
    data.rule.supp <- get.rule.support(data, rules)
    rules.out <- get.rule.outputs(data.rule.supp, y)
    
    # symbolic paths
    paths <- lapply(paths, function(path){lapply(path, function(split){
      if (bins.list[[as.numeric(split[1])]]$type %in% c('continuous', 'discrete')){
        split <- split[1:3]}; return(split)})
    })
    
    # for regression: fit rule coefficients
    if (type == 'reg' & length(paths) > 1){
      data.rule <- sapply(1:length(paths), function(j){
        X <- data.rule.supp[, j]
        X[X == 1] <- rules.out[[j]]$outputs[1]
        X[X == 0] <- rules.out[[j]]$outputs[2]
        X
      })
      lambda <- cv.glmnet(data.rule, y, lower.limits = rep(0, ncol(data.rule)),
                          alpha = 0, standardize = F)$lambda.min
      rule.glm <- glmnet(data.rule, y, lower.limits = rep(0, ncol(data.rule)),
                         lambda = lambda, alpha = 0, standardize = F)
      rule.weights <- as.vector(rule.glm$beta)
    }else{
      rule.glm <- NULL
      rule.weights <- rep(1/length(rules), length(rules))
    }
    
  }else{
    
    rules <- list()
    rules.out <- list()
    rule.glm <- NULL
    rule.weights <- 1
    
  }
  
  return(list(rules = rules, rules.out = rules.out, proba = proba, paths = paths,
              rule.weights = rule.weights, rule.glm = rule.glm, type = type, 
              num.trees = num.trees, data.names = data.names, mean = mean.out, bins = bins.list))
  
}


##### Print the list of rules output by Spatial SIRUS
s.sirus.print <- function(s.sirus.m, digits = 3){
  
  # check s.sirus.m is a valid Spatial SIRUS model
  s.sirus.m.valid <- s.sirus.model.check(s.sirus.m)
  if (!s.sirus.m.valid){
    stop('Invalid Spatial SIRUS model.')
  }
  digits.valid <- is.numeric(digits) & digits > 0
  if (!digits.valid){
    stop('Invalid number of digits.')
  }
  
  rules <- s.sirus.m$rules
  rules.out <- s.sirus.m$rules.out
  rule.weights <- sapply(s.sirus.m$rule.weights, signif, digits = digits)
  
  # format Spatial SIRUS output in a readable format
  if (length(rules) > 0){
    rules.print <- paste0(lapply(1:length(rules), function(ind) {
      rule <- rules[[ind]]
      rule.paste <- paste0(lapply(rule, function(split){
        if (split[2] %in% c('<', '>=')){
          split[3] <- signif(as.numeric(split[3]), digits = digits)
          split.paste <- paste0(split, collapse = ' ')
        }
        if (split[2] == '='){
          split.paste <- paste0(split[1], ' in {', paste0(split[3:length(split)], collapse = ', '), '}')
        }
        split.paste
      }), collapse = ' & ')
      rule.out <- rules.out[[ind]]
      out.true <- signif(rule.out$outputs[1], digits = digits)
      out.false <- signif(rule.out$outputs[2], digits = digits)
      size.true <- rule.out$supp.size[1]
      size.false <- rule.out$supp.size[2]
      weight <- signif(rule.weights[ind], digits = digits)
      paste0(c('if ', rule.paste, ' then ', out.true, ' (n=', size.true, ') else ', out.false, ' (n=', size.false, ')'), collapse = '')
    }))
    
      rules.print <- rules.print[rule.weights > 0]
      rule.weights <- rule.weights[rule.weights > 0]
      mean.print <- paste0('Mean of output y = ', signif(s.sirus.m$mean, digits = digits),
                           ' - Sample size n = ', sum(rules.out[[1]]$supp.size))
      rules.print <- c(mean.print, rules.print)
      if (length(rules) > 1){intercept <- s.sirus.m$rule.glm$a0}else{intercept <- 0}
      intercept <- paste('Intercept =', signif(intercept, digits = digits))
      rule.weights <- c(intercept, rule.weights)
      rules.print <- cbind(rule.weights, rules.print)
      colnames(rules.print) <- c('Weights', 'Rules')
    
  }else{
    out <- signif(s.sirus.m$mean, digits = digits)
    rules.print <- paste0('Empty rule set. Constant output = ', out)
  }
  
  return(rules.print)
}


##### Compute Spatial SIRUS predictions for the large scale (aka regression function)
s.sirus.predict <- function(s.sirus.m, data.test){
  
  # check s.sirus.m is a valid Spatial SIRUS model
  s.sirus.m.valid <- s.sirus.model.check(s.sirus.m)
  if (!s.sirus.m.valid){
    stop('Invalid Spatial SIRUS model.')
  }
  # check data.test
  data.check(data.test)
  # check names
  if (!is.null(colnames(data.test))){
    names.valid <- all(colnames(data.test) == s.sirus.m$data.names)
  }else{names.valid <- FALSE}
  if (!names.valid){
    stop('Invalid variable names for data.test.')
  }
  
  rules <- s.sirus.m$rules
  rules.out <- s.sirus.m$rules.out
  
  if (length(rules) > 0){
    data.rule <- get.data.rule(data.test, rules, rules.out)
    if (s.sirus.m$type == 'classif' | length(rules) == 1){
      pred <- apply(data.rule, 1, mean)
    }else{
      pred <- as.vector(predict(s.sirus.m$rule.glm, data.rule))
    }
  }else{
    pred <- rep(s.sirus.m$mean, nrow(data.test))
  }
  
  return(pred)
  
}


##### Tune the optimal hyperparameter p0 used to select the most frequent rules in s.sirus.fit
# Differently from the original function in the SIRUS package, s.sirus.cv requires:
# data = matrix of the features, the first must be the coordinates
# y = vector of the response values
# incl.coords = logic to consider coordinates provided in data as predictors (TRUE) or not (FALSE)
# type = type of application, 'reg' for regression
# max.depth = must be equal to 2

s.sirus.cv <- function (data, y, incl.coords = NULL, type = "reg", nfold = 10, ncv = 10, num.rule.max = 25, 
                        q = 10, discrete.limit = 10, num.trees.step = 1000, alpha = 0.05, 
                        mtry = NULL, max.depth = 2, num.trees = NULL, num.threads = NULL, 
                        replace = TRUE, sample.fraction = NULL, verbose = TRUE, seed = NULL) {
  
  data.y.check(data, y)
  nfold.valid <- is.numeric(nfold)
  if (nfold.valid) {
    nfold.valid <- (round(nfold) == nfold) & nfold >= 2
  }
  if (!nfold.valid) {
    stop("Invalid nfold. Number of cross-validation folds has to be an integer greater than 2.")
  }else{
    if (nfold > nrow(data)) {
      nfold <- nrow(data)
      warning(paste0("Warning nfold: nfold is greater than the sample size (=", 
                     nrow(data), "),\n                      and is then automatically set to ", 
                     nfold, "."))
    }
  }
  type.valid <- type %in% c("reg")
  if (type.valid) {
    if (type == "reg") {
      if (length(unique(y)) == 2) {
        warning("Warning type: y takes only two distinct values, classification is more appropriate.")
      }
    }
  }else{
    stop("Invalid type. Type should be reg (for regression).")
  }
  ncv.valid <- is.numeric(ncv)
  if (ncv.valid) {
    ncv.valid <- (round(ncv) == ncv) & ncv >= 1
  }
  if (!ncv.valid) {
    stop("Invalid ncv. Number of cross-validations has to be an integer greater than 1.")
  }else{
    if (ncv == 1) {
      warning("Warning ncv: It is recommended to run multiple cross-validations for a robust estimation of p0.")
    }
  }
  s.sirus.param.check(data, num.rule.max, q, num.trees.step, 
                    alpha, mtry)
  
  if (!is.null(mtry)){
    mtry.ratio <- mtry/ncol(data)
  }
  if (is.null(mtry) & incl.coords == TRUE){
    mtry.ratio <- 1/3
    mtry <- max(1, floor((ncol(data))*mtry.ratio))
  }
  if (is.null(mtry) & incl.coords == FALSE){
    mtry.ratio <- 1/3
    mtry <- max(1, floor((ncol(data)-2)*mtry.ratio))
  }
  
  if (is.null(sample.fraction)) {
    sample.fraction <- ifelse(replace, 1, 0.632)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # p0.grid
  p0.grid <- exp(seq(log(mtry.ratio), log(1/1000), length.out = 500))
  
  # run cross-validation
  error.grids <- lapply(1:ncv, function(iter, p0.grid) { 
    
    if (verbose == TRUE) {
      print(paste0("Running cross-validation iteration ", iter, "/", 
                   ncv, " ..."))
    }
    
    # create K-folds
    ndata <- nrow(data)
    ind <- cut(seq(1,ndata), breaks = nfold, labels = F)
    ind <- sample(ind, size = ndata, replace = F)
    folds.ind <- lapply(1:nfold, function(fold, ind){
      test <- which(ind == fold)
      train <- setdiff(1:ndata, test)
      return(list(train = train, test = test))
    }, ind  = ind)
    
      num.rule.cv <- 2*num.rule.max
    
    # fit Spatial SIRUS for each fold and compute prediction for all rule selections
    pred.cv <- lapply(1:nfold, function(fold){
      
      if (verbose == TRUE){
        print(paste0('~~~~~~~~~~ Running SPATIAL SIRUS for fold ', fold, '/', nfold, ' ...'))
      }
      
      data.train <- data[folds.ind[[fold]]$train, , drop = F]
      y.train <- y[folds.ind[[fold]]$train]
      
      check_out_folds.ind = list(folds.ind=folds.ind)
      
      s.sirus.cv <- s.sirus.fit(data.train, y.train, incl.coords, type = type, num.rule = num.rule.cv, p0 = NULL, q = q,
                              num.trees.step = num.trees.step, alpha = alpha, mtry = mtry, max.depth = max.depth,
                              num.trees = num.trees, num.threads = num.threads, replace = replace,
                              sample.fraction = sample.fraction, verbose = FALSE, seed = seed)
      
      check_out_s.sirus.cv = list(s.sirus.cv = s.sirus.cv)
      
      data.test <- data[folds.ind[[fold]]$test, , drop = F]
      
      if (incl.coords == FALSE){
        data.test = data.test[,-c(1:2)]
        data.train = data.train[,-c(1:2)]
      }
      
      if (length(s.sirus.cv$rules) > 0) {
        data.rule.test <- get.data.rule(data.test, s.sirus.cv$rules, 
                                        s.sirus.cv$rules.out)
        if (type == "reg") {
          data.rule.train <- get.data.rule(data.train, 
                                           s.sirus.cv$rules, s.sirus.cv$rules.out)
        }
      }
      else {
        data.rule.test <- as.data.frame(rep(s.sirus.cv$mean, 
                                            nrow(data.test)))
        if (type == "reg") {
          data.rule.train <- as.data.frame(rep(s.sirus.cv$mean, 
                                               nrow(data.train)))
        }
      }
      
      pred.df <- lapply(1:ncol(data.rule.test), function(ind) {
          if (ind > 1) {
            
            lambda <- cv.glmnet(data.rule.train[, 1:ind, 
                                                drop = F], y.train, lower.limits = rep(0, 
                                                                                       ind), alpha = 0, standardize = F)$lambda.min
            s.sirus.glm <- glmnet(data.rule.train[, 1:ind, 
                                                drop = F], y.train, lower.limits = rep(0, 
                                                                                       ind), lambda = lambda, alpha = 0, standardize = F)
            beta.tmp <- rep(0, num.rule.cv)
            beta.tmp[1:ind][as.vector(s.sirus.glm$beta) > 
                              0] <- 1
            pred <- predict(s.sirus.glm, data.rule.test[, 
                                                      1:ind, drop = F])
            beta <- beta.tmp
          }
          else {
            beta.tmp <- c(1, rep(0, num.rule.cv - 1))
            pred <- data.rule.test[, 1, drop = F]
            beta <- beta.tmp
          }
        list(pred = pred, beta = beta)
      })
      beta.df <- lapply(pred.df, function(pred.ind) {
        pred.ind$beta
      })
      pred.df <- matrix(sapply(pred.df, function(pred.ind) {
        pred.ind$pred
      }), nrow = nrow(data.rule.test))
      pred.df <- cbind(rep(s.sirus.cv$mean, nrow(pred.df)), 
                       pred.df)
      proba.cv <- s.sirus.cv$proba
      paths.cv <- s.sirus.cv$paths
      return(list(pred.df = pred.df, beta.df = beta.df, 
                  proba = proba.cv, paths = paths.cv))
    })
    proba.cv <- lapply(pred.cv, function(pred) {
      pred$proba
    })
    pred.df <- lapply(pred.cv, function(pred) {
      pred$pred.df
    })
    y.test <- unlist(lapply(1:nfold, function(fold) {
      y[folds.ind[[fold]]$test]
    }))
    folds.ncol <- sapply(1:nfold, function(fold) {
      proba <- proba.cv[[fold]]
      proba <- c(proba, 0)
      fold.ncol <- rep(0, length(p0.grid))
      for (ind in 2:(min(length(proba), num.rule.cv + 1))) {
        fold.ncol[proba[ind] <= p0.grid & proba[ind - 
                                                  1] > p0.grid] <- ind - 1
      }
      fold.ncol
    })
    error.grid <- unlist(lapply(1:length(p0.grid), function(ind) {
      pred <- unlist(lapply(1:nfold, function(fold) {
        pred.df[[fold]][, folds.ncol[ind, fold] + 1]
      }))
        sum((pred - y.test)^2)/sum((y.test - mean(y.test))^2)
    }))
    stab.df <- lapply(1:(nfold - 1), function(fold1) {
      pred.cv1 <- pred.cv[[fold1]]
      lapply((fold1 + 1):nfold, function(fold2) {
        pred.cv2 <- pred.cv[[fold2]]
        sapply(p0.grid, function(p0) {
          if (sum(pred.cv1$proba > p0)) {
            paths1 <- pred.cv1$paths[pred.cv1$proba > 
                                       p0][pred.cv1$beta.df[[folds.ncol[which(p0.grid == 
                                                                                p0), fold1]]] > 0]
          }
          else {
            paths1 <- list()
          }
          if (sum(pred.cv2$proba > p0)) {
            paths2 <- pred.cv2$paths[pred.cv2$proba > 
                                       p0][pred.cv2$beta.df[[folds.ncol[which(p0.grid == 
                                                                                p0), fold2]]] > 0]
          }
          else {
            paths2 <- list()
          }
          len <- (length(paths1) + length(paths2))/2
          if (len > 0) {
            c(length(intersect(paths1, paths2))/len, 
              len)
          }
          else {
            c(1, 0)
          }
        })
      })
    })
    num.rule.df <- lapply(1:(nfold - 1), function(fold1) {
      sapply(1:(nfold - fold1), function(fold2) {
        stab.df[[fold1]][[fold2]][2, ]
      })
    })
    num.rule.df <- do.call("cbind", num.rule.df)
    num.rules <- apply(num.rule.df, 1, mean)
    stab.df <- lapply(1:(nfold - 1), function(fold1) {
      sapply(1:(nfold - fold1), function(fold2) {
        stab.df[[fold1]][[fold2]][1, ]
      })
    })
    stab.df <- do.call("cbind", stab.df)
    stab.grid <- apply(stab.df, 1, mean)
    return(list(error.grid = error.grid, stab.grid = stab.grid, 
                num.rules = num.rules))
  }, p0.grid = p0.grid)
  error.df <- sapply(error.grids, function(x) {
    x$error.grid
  })
  error.mean <- apply(error.df, 1, mean)
  if (ncv > 1) {
    error.sd <- apply(error.df, 1, sd)/sqrt(ncv)
  }
  else {
    error.sd <- rep(0, length(p0.grid))
  }
  stab.df <- sapply(error.grids, function(x) {
    x$stab.grid
  })
  stab.mean <- apply(stab.df, 1, mean)
  if (ncv > 1) {
    stab.sd <- apply(stab.df, 1, sd)/sqrt(ncv)
  }
  else {
    stab.sd <- rep(0, length(p0.grid))
  }
  num.rules.df <- sapply(error.grids, function(x) {
    x$num.rules
  })
  num.rules.mean <- apply(num.rules.df, 1, mean)
  if (ncv > 1) {
    num.rules.sd <- apply(num.rules.df, 1, sd)/sqrt(ncv)
  }
  else {
    num.rules.sd <- rep(0, length(p0.grid))
  }
  error.grid.p0 <- as.data.frame(cbind(p0.grid, num.rules.mean, 
                                       stab.mean, error.mean, num.rules.sd, stab.sd, error.sd))
  ind.1 <- which.min(abs(1 - num.rules.mean))[1]
  ind.max <- which.min(abs(min(num.rule.max, max(num.rules.mean)) - 
                             num.rules.mean))
  error.grid.p0 <- error.grid.p0[ind.1:ind.max, ]
  ind.min <- ind.1 - 1 + which.min(error.mean[ind.1:ind.max])
  ind.pred <- max(ind.1, min(which(error.mean <= error.mean[ind.min] + 
                                     2 * error.sd[ind.min])))
  p0.stab.cv <- sapply(error.grids, function(grid.cv) {
    criterion <- (grid.cv$stab.grid[ind.1:ind.max] - 0.9)^2 + 
      (0 - grid.cv$error.grid[ind.1:ind.max])^2
    p0.grid[ind.1:ind.max][which.min(criterion)]
  })
  p0.stab <- median(p0.stab.cv)
  return(list(p0.pred = p0.grid[ind.pred], p0.stab = p0.stab, 
              error.grid.p0 = error.grid.p0, type = type))
}


##### Plot the Spatial SIRUS cross-validation path: error and stability versus the number of rules when p0 varies.
s.sirus.plot.cv <- function(s.sirus.cv.grid, p0.criterion = NULL, num.rule.max = 25){
  
  if (is.null(p0.criterion)){
    if (s.sirus.cv.grid$type == 'reg'){p0.criterion <- 'stab'}else{p0.criterion <- 'pred'}
  }else{
    if (!p0.criterion %in% c('pred', 'stab')){
      stop('Invalid p0 criterion. Its value has to be pred or stab.')
    }
  }
  
  # max num of rules
  num.rule.valid <- is.numeric(num.rule.max) & num.rule.max > 0
  if (!num.rule.valid){
    stop('Invalid maximum number of rules. Its value should be a positive integer.')
  }
  
  # filter cv grid & performance metrics
  error.grid.p0 <- s.sirus.cv.grid$error.grid.p0
  num.rule.max <- min(num.rule.max, max(error.grid.p0$num.rules.mean))
  grid.index <- sapply(1:num.rule.max, function(ind){which.min(abs(ind - s.sirus.cv.grid$error.grid.p0$num.rules.mean))[1]})
  if (p0.criterion == 'stab'){
    p0 <- s.sirus.cv.grid$p0.stab
  }else{
    p0 <- s.sirus.cv.grid$p0.pred
  }
  ind.p0 <- which.min(abs(s.sirus.cv.grid$error.grid.p0$p0.grid - p0))
  grid.cv <- s.sirus.cv.grid$error.grid.p0[ind.p0,]
  num.rules <- grid.cv$num.rules.mean
  error <- grid.cv$error.mean
  stab <- grid.cv$stab.mean
  error.grid.p0 <- s.sirus.cv.grid$error.grid.p0[sort(unique(c(grid.index,ind.p0))),]
  
  # declare variables
  num.rules.mean <- NULL
  error.mean <- NULL
  error.sd <- NULL
  stab.mean <- NULL
  stab.sd <- NULL
  
  # plot error vs number of rules
  label.error <- if (s.sirus.cv.grid$type == 'reg'){'Unexplained variance'}else{'1-AUC'}
  tag.y <- (min(error.grid.p0$error.mean) + max(error.grid.p0$error.mean))/2
  plot.error <- ggplot(error.grid.p0, aes(x = num.rules.mean, y = error.mean)) +
    geom_line(linewidth = 0.8) + geom_point(size=1) +
    geom_errorbar(aes(ymin=error.mean - 2*error.sd,
                      ymax=error.mean + 2*error.sd),  width=0.7, linewidth=0.5)+
    geom_hline(yintercept = error,
               linetype = 'dashed', color = 'blue', linewidth = 0.7) +
    geom_vline(xintercept = num.rules,
               linetype = 'dashed', color = 'blue', linewidth = 0.7) +
    annotate("text", x=num.rules, y=tag.y, label = 'Optimal p0', angle = '90',
              vjust=1.5, color='blue', size = 7) +
    annotate("text", x=num.rule.max - 5, y=error, label = paste0('S-SIRUS Error: ', round(error,2)),
              vjust=-1, color='blue', size = 7) +
    xlab('Number of rules') +
    ylab(label.error) +
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'black'),
          axis.line.y = element_line(colour = 'black'),
          text = element_text(size=20),
          plot.title = element_text(hjust = 0.5, size=18, face="italic"))
  
  # plot stability vs number of rules
  stab.extremum <- c(min(error.grid.p0$stab.mean), max(error.grid.p0$stab.mean))
  tag.p0 <- (stab.extremum[which.max(abs(stab.extremum - stab))] + stab)/2
  plot.stab <- ggplot(error.grid.p0, aes(x = num.rules.mean, y = stab.mean)) +
    geom_line(linewidth = 0.8) + geom_point() +
    geom_errorbar(aes(ymin=stab.mean - stab.sd,
                      ymax=stab.mean + stab.sd),  width=0.7, linewidth=0.5)+
    geom_hline(yintercept = stab,
               linetype = 'dashed', color = 'blue', linewidth = 0.7) +
    geom_vline(xintercept = num.rules,
               linetype = 'dashed', color = 'blue', linewidth = 0.7) +
    annotate("text", x = num.rules, y = tag.p0, label = "Optimal p0",
             angle = '90', vjust=1.5, color='blue', size = 7) +
    annotate("text",x = num.rules - 3, y = stab, label = paste0('S-SIRUS Stability: ', round(stab,2)),
              vjust=-1, color='blue', size = 7) +
    xlab('Number of rules') +
    ylab('Stability') +
    theme_classic() + ylim(c(0.9*min(error.grid.p0$stab.mean), 1.1*max(error.grid.p0$stab.mean))) +
    theme(axis.line.x = element_line(colour = 'black'),
          axis.line.y = element_line(colour = 'black'),
          text = element_text(size=20),
          plot.title = element_text(hjust = 0.5, size=18, face="italic"))
  
  return(list(error = plot.error, stability = plot.stab))
  
}
