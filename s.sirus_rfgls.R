# The following code is an adaptation of the SIRUS functions available at https://github.com/cran/sirus/tree/master/R .
# This code can be used for the implementation of the functions for Spatial SIRUS (S-SIRUS).

##### Run RFGLS and prepare the paths from the grown trees

rfgls <- function(formula = NULL, data = NULL, num.trees = 500, mtry = NULL,
                     importance = "none", write.forest = TRUE, probability = FALSE,
                     min.node.size = NULL, max.depth = NULL, replace = TRUE, 
                     sample.fraction = ifelse(replace, 1, 0.632), 
                     case.weights = NULL, class.weights = NULL, splitrule = NULL, 
                     num.random.splits = 1, alpha = 0.5, minprop = 0.1,
                     split.select.weights = NULL, always.split.variables = NULL,
                     respect.unordered.factors = NULL,
                     scale.permutation.importance = FALSE,
                     keep.inbag = FALSE, inbag = NULL, holdout = FALSE,
                     quantreg = FALSE, oob.error = TRUE,
                     num.threads = NULL, save.memory = FALSE,
                     verbose = TRUE, seed = NULL, 
                     dependent.variable.name = NULL, status.variable.name = NULL, 
                     classification = NULL,
                     coords = NULL,
                     incl.coords = NULL) {
  
  ## Formula interface. Use whole data frame if no formula provided and depvarname given
  if (is.null(formula)) {
    if (is.null(dependent.variable.name)) {
      stop("Error: Please give formula or dependent variable name.")
    }
    if (is.null(status.variable.name)) {
      status.variable.name <- ""
      response <- data[, dependent.variable.name, drop = TRUE]
    } else {
      response <- survival::Surv(data[, dependent.variable.name], data[, status.variable.name]) #data[, c(dependent.variable.name, status.variable.name)]
    }
    data.selected <- data
  } else {
    formula <- formula(formula)
    if (!inherits(formula, "formula")) {
      stop("Error: Invalid formula.")
    }
    data.selected <- parse.formula(formula, data, env = parent.frame())
    response <- data.selected[, 1]
  }
  
  ## Check missing values
  if (any(is.na(data.selected))) {
    offending_columns <- colnames(data.selected)[colSums(is.na(data.selected)) > 0]
    stop("Missing data in columns: ",
         paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }
  
  treetype <- 3
  
  if (!is.null(formula)) {
    dependent.variable.name <- names(data.selected)[1]
    status.variable.name <- ""
    independent.variable.names <- names(data.selected)[-1]
  } else {
    independent.variable.names <- colnames(data.selected)[colnames(data.selected) != dependent.variable.name &
                                                            colnames(data.selected) != status.variable.name]
  }
  
  ## Recode characters as factors and recode factors if 'order' mode
  if (!is.matrix(data.selected) && !inherits(data.selected, "Matrix")) {
    character.idx <- sapply(data.selected, is.character)
    
    ## Recode characters only
    data.selected[character.idx] <- lapply(data.selected[character.idx], factor)
  }
  
  data.final <- data.matrix(cbind(coords,data.selected))

  variable.names <- colnames(data.final)
  
  all.independent.variable.names <- independent.variable.names
  
  ## Error if no covariates
  if (length(all.independent.variable.names) < 1) {
    stop("Error: No covariates found.")
  }
  
  ## Number of trees
  if (!is.numeric(num.trees) || num.trees < 1) {
    stop("Error: Invalid value for num.trees.")
  }
  
  ## mtry
  if (is.null(mtry)) {
    mtry <- 0
  } else if (!is.numeric(mtry) || mtry < 0) {
    stop("Error: Invalid value for mtry")
  }
  
  ## Seed
  if (is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }
  
  ## Num threads
  ## Default 0 -> detect from system in C++.
  if (is.null(num.threads)) {
    num.threads = 0
  } else if (!is.numeric(num.threads) || num.threads < 0) {
    stop("Error: Invalid value for num.threads")
  }
  
  ## Minumum node size
  if (is.null(min.node.size)) {
    min.node.size <- 0
  } else if (!is.numeric(min.node.size) || min.node.size < 0) {
    stop("Error: Invalid value for min.node.size")
  }
  
  ## Tree depth
  if (is.null(max.depth)) {
    max.depth <- 0
  } else if (!is.numeric(max.depth) || max.depth < 0) {
    stop("Error: Invalid value for max.depth. Please give a positive integer.")
  }
  
  splitrule <- "variance"
  splitrule.num <- 1
  
  ## Maxstat splitting
  if (alpha < 0 || alpha > 1) {
    stop("Error: Invalid value for alpha, please give a value between 0 and 1.")
  }
  
  ## Clean up
  rm("data.selected")
  
  y = data.final[,3]
  
  if(incl.coords == TRUE){
    X = as.matrix(data.final[,-c(1:3)]) # Coords included in large scale!!!
  }
  if(incl.coords == FALSE){
    X = as.matrix(data.final[,-c(1:5)]) # NOT Coords included in large scale!!!
  }
  
  suppressWarnings({
    result = RFGLS_estimate_spatial(coords = as.matrix(coords),
                                    y = y,
                                    X = X, Xtest = X,
                                    nthsize = 25, mtry = mtry, ntree = num.trees,
                                    cov.model = "exponential",
                                    h = 1,
                                    param_estimate = TRUE, verbose = FALSE) })
  
  tree_List_RFGLS <- NULL
  tree_List_RFGLS$ntree <- num.trees
  tree_List_RFGLS$list <- vector("list", num.trees)
  
  
  result$RFGLS_object$ndbigtree = c()
  
  for (i in 1:tree_List_RFGLS$ntree){
    result$RFGLS_object$ndbigtree[i] = length(which(result$RFGLS_object$nodestatus[,i] !=0))
  }
  
  col_names = c("left daughter", "right daughter",  "split var", "split point", "status", "prediction")
  
  for (i in 1:tree_List_RFGLS$ntree) {
    tree_List_RFGLS$list[[i]] <- cbind(result$RFGLS_object$ldaughter[,i],
                                       result$RFGLS_object$rdaughter[,i],
                                       result$RFGLS_object$mbest[,i],
                                       result$RFGLS_object$upper[,i],
                                       result$RFGLS_object$nodestatus[,i],
                                       result$RFGLS_object$avnode[,i]
    )[1:result$RFGLS_object$ndbigtree[i], ]
    colnames(tree_List_RFGLS$list[[i]]) = col_names
  }
  
  
  paths_all = list()
  paths = list()
  
  for(i in 1:num.trees){
    
    app = as.data.frame(tree_List_RFGLS$list[[i]])
    app = app[1:(max.depth + 1),]
    
    paths[[1]] <- list(c(app[1,3],app[1,4],0))
    paths[[2]] <- list(c(app[1,3],app[1,4],1))
    
    paths[[3]] <- list(unlist(paths[[1]]),c(app[2,3],app[2,4],0))
    paths[[4]] <- list(unlist(paths[[1]]),c(app[2,3],app[2,4],1))
    
    paths[[5]] <- list(unlist(paths[[2]]),c(app[3,3],app[3,4],0))
    paths[[6]] <- list(unlist(paths[[2]]),c(app[3,3],app[3,4],1))
    
    paths_all = c(paths_all,paths)
    
  }
  
  
  paths_all = paths_all[c(which(!(fun_search(paths_all, list(c(0,0,0))) == TRUE)))]
  paths_all = paths_all[c(which(!(fun_search(paths_all, list(c(0,0,1))) == TRUE)))]
  
  
  paths_all = c(list(list(c(0,0,0))),paths_all)
  
  paths.proba = c()
  paths.proba[1] = num.trees
  
  paths_all_collapsed = lapply(paths_all, function(x) {paste0(x, collapse = "")})
  
  for (i in 2:length(paths_all_collapsed)){
    paths.proba = c(paths.proba,fun_search_count(paths_all_collapsed, paths_all_collapsed[[i]]))
  }
  
  path_all_count = cbind(unlist(paths_all_collapsed),paths.proba)
  if(is_empty(which(duplicated(path_all_count)))){
    paths.proba = as.integer(path_all_count[,2])
  }else{
    paths.proba = as.integer(path_all_count[-c(which(duplicated(path_all_count))),2])
  } 
  
  if(is_empty(which(duplicated(path_all_count)))){
    paths_all = paths_all
  }else{
    paths_all = paths_all[-c(which(duplicated(path_all_count)))]
  }
  
  y_hat = result$Predicted_Test
  prediction.error = mean((y_hat - y)^2)
  
  result = list(as.double(result$Predicted_Test),num.trees, num.independent.variables = ncol(X), # ncol(data.final)-1,
                mtry, min.node.size = 20,prediction.error, paths_all, paths.proba)
  
  names(result) = c("predictions","num.trees","num.independent.variables","mtry",
                    "min.node.size","prediction.error","paths","paths.proba")
  
  
  if (length(result) == 0) {
    stop("User interrupt or internal error.")
  }
  
  ## Splitrule
  result$splitrule <- splitrule
  
  result$treetype <- "Regression"

  result$call <- sys.call()
  result$importance.mode <- importance
  result$num.samples <- nrow(data.final)
  result$replace <- replace
  
  class(result) <- "ranger"
  
  
  return(result)
}

