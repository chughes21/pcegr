
#' The Check For Clean Data Function
#'
#' The function [ceg::CheckForCleanData]
#'
#' @param data.frame A data set
#'
#' @return A logical value indicating whether the data set is clean, along with a message if it isn't
#'
#' @examples
CheckForCleanData1<-function(data.frame){
  out <- TRUE
  data.frame <- as.data.frame(apply(data.frame, 2, function(x) gsub("^\\s*$",
                                                                    NA, x)))
  if (anyNA(data.frame)) {
    message("Your data contains rows with <NA>, absent, or whitespace  values")
    out <- FALSE
  }
  return(out)
}

#' The Check And Clean Data Function
#'
#' The function [ceg::CheckAndCleanData]
#'
#' @param data.frame A data set
#'
#' @return A cleaned data set
#'
#' @examples
CheckAndCleanData1<-function(data.frame){
  if (!methods::is(data.frame, "data.frame")) {
    message("Input is not a data frame, please check it")
    return(NULL)
  }
  data.frame <- as.data.frame(apply(data.frame, 2, function(x) gsub("^\\s*$",
                                                                    NA, x)))
  if (anyNA(data.frame)) {
    message("The data frame has rows with <NA> or absent values (\"\")")
    message("All these rows were removed")
  }
  else message("Your data does not contain rows with <NA>, absent, or whitespace  values")
  out <- data.frame[!(rowSums(is.na(data.frame))), ]
  return(out)
}

#' The Label Stage Function
#'
#' The function [ceg::LabelStage]
#'
#'
#' @param k
#' @param num.variable
#' @param num.situation
#' @param label.category
#' @param num.category
#'
#' @examples
LabelStage1<-function(k, num.variable, num.situation, label.category, num.category){
  if (k > num.variable) {
    var <- k - num.variable
  }
  else {
    var <- k
  }
  label <- c(1, rep("NA", sum(num.situation[1:k]) - 1))
  label <- c(label, rep(label.category[[var]], num.situation[k]))
  if (k < (num.variable - 1)) {
    label <- c(label, TruncatedPath1(num.variable, k, var,
                                     num.category, num.situation, label.category))
  }
  return(label)
}

#' The Label Stage Function
#'
#' The function [ceg::TruncatedPath]
#'
#' @param ref
#' @param k
#' @param var
#' @param num.category
#' @param num.situation
#' @param label.category
#'
#' @return
#'
#' @examples
TruncatedPath1<-function (ref, k, var, num.category, num.situation, label.category)
{
  if (ref < k + 2)
    return(c())
  return(c(TruncatedPath1(ref - 1, k, var, num.category, num.situation,
                          label.category), rep(label.category[[var]], each = num.situation[ref]/num.situation[k +
                                                                                                                1], num.situation[k + 1]/num.category[var])))
}

#' The Stratified Event Tree Function
#'
#' This is a slightly modified version of the [ceg::stratified.event.tree()] constructor method.
#'
#' @param x A data frame
#' @param ...
#'
#' @return A Statified.event.tree object
#' @export
#'
#' @examples
#' set(knee_pain_obs[,1:3])
set<-function(x = "data.frame", ...)
{
  .local <- function (x = "data.frame")
  {
    data.frame <- x
    if (!CheckForCleanData1(data.frame)) {
      stop("Consider using CheckAndCleanData1() function")
    }
    num.variable <- length(data.frame[1, ])
    num.slice <- ncol(data.frame)/num.variable
    label.category <- lapply(1:num.variable, function(x) levels(data.frame[, x])) #changed for better labels
    num.category <- c()
    num.category <- sapply(label.category, length)
    num.situation <- c(1, cumprod(num.category[1:(num.variable -  1)]))
    begin.stage <- c(0, cumsum(num.situation[1:(num.variable -  1)]))
    mergedlist <- sapply(1:(num.variable - 1), function(x) LabelStage1(x, num.variable, num.situation, label.category, num.category))
    mergedlist <- lapply(1:(num.variable), function(x) {
      lapply(seq_len(num.situation[x]), function(y) mergedlist[y + begin.stage[x], ])
    })
    hyper.stage <- list()
    return(new("Stratified.event.tree", num.variable,
               num.category, label.category, num.situation, num.slice,
               mergedlist, hyper.stage))
  }
  .local(x, ...)
}

#' The StratifiedCEGPosition1 Function
#'
#' @param stage
#' @param num.category
#' @param num.situation
#'
#' @return
#'
#' @examples
StratifiedCegPosition1<-function (stage, num.category, num.situation)
{
  num.level <- length(num.category)
  result <- list()
  length(result) <- num.level
  result[[num.level]] <- PositionLevel1(stage[[num.level]]@cluster,
                                        0, num.situation[num.level])
  for (level in (num.level - 1):2) {
    result[[level]] <- PositionLevel1(stage[[level]]@cluster,
                                      num.category[level], num.situation[level + 1], result[[level + 1]])
  }
  result[[1]] <- list(list(1))
  return(result)
}

#' The Position Level Function
#'
#' @param stage.list
#' @param num.category
#' @param num.situation.next
#' @param pos.next.level
#'
#' @return
#'
#' @examples
PositionLevel1<-function(stage.list, num.category, num.situation.next, pos.next.level = list())
{
  aux <- which(!is.na(stage.list))
  N <- length(aux)
  if (num.category == 0) {
    stage.list <- lapply(1:N, function(x) list(stage.list[[aux[x]]]))
    return(stage.list)
  }
  stage.list <- lapply(1:N, function(x) stage.list[[aux[x]]])
  pos.next.level <- PositionVector1(num.situation.next, pos.next.level)
  result <- lapply(1:N, function(x) PositionStage1(stage.list[[x]],
                                                   num.category, pos.next.level))
  return(result)
}

#' The Position Vector Function
#'
#' @param num.situation
#' @param pos.list
#'
#' @return
#' @export
#'
#' @examples
PositionVector1<-function (num.situation, pos.list)
{
  num.situation <- length(pos.list)
  pos.vec <- rep(0, num.situation)
  count <- 1
  for (stage in 1:num.situation) {
    for (pos in 1:length(pos.list[[stage]])) {
      pos.vec[pos.list[[stage]][[pos]]] <- count
      count <- count + 1
    }
  }
  return(pos.vec)
}

#' The Position Stage Function
#'
#' @param stage.vector
#' @param num.category
#' @param pos.next.level
#'
#' @return
#'
#' @examples
PositionStage1<-function (stage.vector, num.category, pos.next.level)
{
  stage.vector <- sort(stage.vector)
  result <- list()
  count <- 1
  stop <- FALSE
  N <- length(stage.vector)
  if (N == 1)
    return(list(stage.vector))
  while (stop == FALSE) {
    aux.stage <- sapply(2:N, function(x) PairwisePosition1(c(stage.vector[1],
                                                             stage.vector[x]), num.category, pos.next.level))
    aux.stage <- c(TRUE, aux.stage)
    result[[count]] <- stage.vector[aux.stage]
    count <- count + 1
    stage.vector <- stage.vector[!aux.stage]
    N <- length(stage.vector)
    if (N == 1) {
      stop <- TRUE
      result[[count]] <- stage.vector
    }
    else if (N == 0)
      stop <- TRUE
  }
  return(result)
}

#' The Pairwise Position Function
#'
#' @param pair.situation
#' @param num.category
#' @param pos.next.level
#'
#' @return
#'
#' @examples
PairwisePosition1<-function (pair.situation, num.category, pos.next.level)
{
  situation.1 <- (pair.situation[1] - 1) * num.category + 1:num.category
  situation.2 <- (pair.situation[2] - 1) * num.category + 1:num.category
  situation.1 <- pos.next.level[situation.1]
  situation.2 <- pos.next.level[situation.2]
  aux <- sum(situation.1 == situation.2)
  if (aux != num.category)
    return(FALSE)
  else return(TRUE)
}

#' The Stratified CEG Function
#'
#' This is the [ceg::stratified.event.tree()] constructor method.
#'
#' @param stratified.staged.tree An object of class "Stratified.staged.tree"
#'
#' @return An object of class "Stratified.ceg.model"
#' @export
#'
#' @examples
#' mod1<-pceg(knee_pain_obs,2,TRUE,TRUE)
#'  staged.tree<-staged.tree.creator(mod1)
#'  mod2<-sceg(staged.tree)
#'  plot(mod2)
sceg<-function (stratified.staged.tree)
{
  position <- StratifiedCegPosition1(stratified.staged.tree@stage.structure,
                                     stratified.staged.tree@event.tree@num.category, stratified.staged.tree@event.tree@num.situation)
  return(new("Stratified.ceg.model", stratified.staged.tree,
             position))
}
