#################
#Function for implementing peeling methods
#"peeling" denotes peeling methods "MP" or "WMP", "formula" takes standard formula input, "coeflearn" takes input "Breiman" (used in the thesis and peeling article) "Freund" and "Zhu"
#This function is inspired by the original boosting() function of the adabag package and repeats the original function based on the resulting noise criteria

PeelBoost<- function (formula, data, mfinal1 = 10, mfinal2 = 200, coeflearn = "Breiman", 
                         control, peeling = "MP", wl=0.5 , ...) 
{
  if (!(as.character(coeflearn) %in% c("Freund", "Breiman", 
                                       "Zhu"))) {
    stop("coeflearn must be 'Freund', 'Breiman' or 'Zhu' ")
  }
  formula <- as.formula(formula)
  vardep <- data[, as.character(formula[[2]])]
  n <- length(data[, 1]) #
  nclases <- nlevels(vardep)
  weights <- rep(1/n, n)
  weightbank <- array(0, c(n, mfinal1))
  w <- rep(1/n, n) # 
  data <- data.frame(weights, data)
  trees <- list() 
  pond <- rep(0, mfinal1)
  pred <- data.frame(rep(0, n)) #
  trees[[1]] <- rpart(formula, data = data[, -1], control = rpart.control(minsplit = 1, 
                                                                          cp = -1, maxdepth = 30))
  nvar <- dim(varImp(trees[[1]], surrogates = FALSE, competes = FALSE))[1] #
  imp <- array(0, c(mfinal1, nvar)) # 
  for (m in 1:mfinal1) { 
    
    w <<- weights 
    fit <- rpart(formula = formula, data = data[, -1], 
                 weights = w, control = control)
    flearn <- predict(fit, data = data[, -1], type = "class")
    ind <- as.numeric(vardep != flearn)
    err <- sum(weights * ind)
    
    c <- log((1 - err)/err)
    if (coeflearn == "Breiman") {
      c <- (1/2) * c
    }
    if (coeflearn == "Zhu") {
      c <- c + log(nclases - 1)
    }
    weightbank[, m] <- weights
    weights <- weights * exp(c * ind)
    weights <- weights/sum(weights)
    maxerror <- 0.5
    eac <- 0.001
    if (coeflearn == "Zhu") {
      maxerror <- 1 - 1/nclases
    }
    if (err >= maxerror) {
      weights <- rep(1/n, n)
      maxerror <- maxerror - eac
      c <- log((1 - maxerror)/maxerror)
      if (coeflearn == "Breiman") {
        c <- (1/2) * c
      }
      if (coeflearn == "Zhu") {
        c <- c + log(nclases - 1)
      }
    }
    if (err == 0) {
      weights <- rep(1/n, n)
      c <- log((1 - eac)/eac)
      if (coeflearn == "Breiman") {
        c <- (1/2) * c
      }
      if (coeflearn == "Zhu") {
        c <- c + log(nclases - 1)
      }
    }
    trees[[m]] <- fit
    pond[m] <- c
    if (m == 1) {
      pred <- flearn
    }
    else {
      pred <- data.frame(pred, flearn)
    }
    if (length(fit$frame$var) > 1) {
      k <- varImp(fit, surrogates = FALSE, competes = FALSE)
      imp[m, ] <- k[sort(row.names(k)), ]
    }
    else {
      imp[m, ] <- rep(0, nvar)
    }
  }
  classfinal <- array(0, c(n, nlevels(vardep)))
  for (i in 1:nlevels(vardep)) {
    classfinal[, i] <- matrix(as.numeric(pred == levels(vardep)[i]), 
                              nrow = n) %*% as.vector(pond)
  } 
  if(peeling == "MP")  
  { 
    margin<- c(rep(0, n))
    margin<- (apply(as.matrix.noquote(pred),2,as.numeric) %*% (pond))/sum(pond)*as.numeric(levels(vardep))[vardep]
    peel<-(margin)>0 
    peelsum <- n - sum(peel) 
    
  } 
  if(peeling == "WMP")  
  { 
    M<- matrix(0, nrow(pred), ncol(pred)) 
    R<- matrix(0, nrow(pred), ncol(pred)) 
    for(i in 1:ncol(pred)){ 
      R[,i] <- pred[,i] == vardep 
    }
    r <- colSums(R) / n # 
    for(i in 1:ncol(pred)){ 
      M[,i] <- pred[,i] != vardep 
    }
    w<-(M%*%r)/sum(r)
    peel <- w<wl #wl=0.5  
    peelsum <-n - sum(peel)
  }
  ###############################################  
  W <- rep(1/(n - peelsum), (n - peelsum))  
  if (peelsum > 0) { 
    data<-data[peel,] 
    n <- length(data[, 1]) 
    formula <- as.formula(formula)
    vardep <- data[, as.character(formula[[2]])] 
    nclases <- nlevels(vardep) # 
    weightbank <- array(0, c(n, mfinal2))
    data <- data.frame(data) 
    weights <- rep(1/n, n)
    data$weights <- weights # 
    trees <- list() 
    pond <- rep(0, mfinal2)
    pred <- data.frame(rep(0, n))
    trees[[1]] <- rpart(formula, data = data[, -1], control = rpart.control(minsplit = 1, 
                                                                            cp = -1, maxdepth = 30))
    nvar <- dim(varImp(trees[[1]], surrogates = FALSE, competes = FALSE))[1] #
    imp <- array(0, c(mfinal2, nvar))
    ###############################################  
    for (m in 1:mfinal2) {
      W <<- weights
      fit <- rpart(formula = formula, data = data[, -1], 
                   weights = W, control = control) #### SISTA GÅNGEN VI SER DEN
      flearn <- predict(fit, data = data[, -1], type = "class")
      ind <- as.numeric(vardep != flearn)
      err <- sum(weights * ind)
      
      c <- log((1 - err)/err)
      if (coeflearn == "Breiman") {
        c <- (1/2) * c
      }
      if (coeflearn == "Zhu") {
        c <- c + log(nclases - 1)
      }
      weightbank[, m] <- weights
      weights <- weights * exp(c * ind)
      weights <- weights/sum(weights)
      maxerror <- 0.5
      eac <- 0.001
      if (coeflearn == "Zhu") {
        maxerror <- 1 - 1/nclases
      }
      if (err >= maxerror) {
        weights <- rep(1/n, n)
        maxerror <- maxerror - eac
        c <- log((1 - maxerror)/maxerror)
        if (coeflearn == "Breiman") {
          c <- (1/2) * c
        }
        if (coeflearn == "Zhu") {
          c <- c + log(nclases - 1)
        }
      }
      if (err == 0) {
        weights <- rep(1/n, n)
        c <- log((1 - eac)/eac)
        if (coeflearn == "Breiman") {
          c <- (1/2) * c
        }
        if (coeflearn == "Zhu") {
          c <- c + log(nclases - 1)
        }
      }
      trees[[m]] <- fit
      pond[m] <- c
      if (m == 1) {
        pred <- flearn
      }
      else {
        pred <- data.frame(pred, flearn)
      }
      if (length(fit$frame$var) > 1) {
        k <- varImp(fit, surrogates = FALSE, competes = FALSE)
        imp[m, ] <- k[sort(row.names(k)), ]
      }
      else {
        imp[m, ] <- rep(0, nvar)
      }
    }
    classfinal <- array(0, c(n, nlevels(vardep)))
    for (i in 1:nlevels(vardep)) {
      classfinal[, i] <- matrix(as.numeric(pred == levels(vardep)[i]), 
                                nrow = n) %*% as.vector(pond)
    }
  }
  ################################################
  predclass <- rep("O", n)
  ########################
  SELECT <- function (fila, vardep.summary, ...) 
  {
    if (length(which(fila == max(fila))) > 1) {
      predclass <- names(vardep.summary[which(fila == max(fila))])[order(vardep.summary[which(fila == 
                                                                                                max(fila))], decreasing = TRUE)[1]]
    }
    else {
      predclass <- as.character(names(vardep.summary)[(order(fila, 
                                                             decreasing = TRUE)[1])])
    }
    predclass
  }
  ######################## 
  predclass[] <- apply(classfinal, 1, FUN = SELECT, vardep.summary = summary(vardep)) # 
  imppond <- as.vector(as.vector(pond) %*% imp)
  imppond <- imppond/sum(imppond) * 100
  names(imppond) <- sort(row.names(k))
  votosporc <- classfinal/apply(classfinal, 1, sum) #
  ans <- list(formula = formula, trees = trees, weights = pond, 
              votes = classfinal, prob = votosporc, class = predclass, 
              importance = imppond, peelsum = peelsum, peel = peel, margins = margin)
  attr(ans, "vardep.summary") <- summary(vardep, maxsum = 700) 
  mf <- model.frame(formula = formula, data = data)
  terms <- attr(mf, "terms") 
  ans$terms <- terms
  ans$call <- match.call()
  class(ans) <- "boosting"
  ans
}   