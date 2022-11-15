#' The Label Stage Function
#'
#' The function [ceg::LabelStage()]
#'
#' @param k An integer value
#' @param num.variable An integer value
#' @param num.situation An integer value
#' @param label.category A string detailing the label of the category
#' @param num.category An integer value
#'
#' @return A vector of labels
#'
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


#' The Truncated Path Function
#'
#' The function [ceg::TruncatedPath()]
#'
#' @param ref An integer value
#' @param k An integer value
#' @param var An integer value
#' @param num.category An integer value
#' @param num.situation An integer value
#' @param label.category A string detailing the label of the category
#'
#' @return A vector
#'
TruncatedPath1<-function (ref, k, var, num.category, num.situation, label.category)
{
  if (ref < k + 2)
    return(c())
  return(c(TruncatedPath1(ref - 1, k, var, num.category, num.situation,
                          label.category), rep(label.category[[var]], each = num.situation[ref]/num.situation[k +
                                                                                                                1], num.situation[k + 1]/num.category[var])))
}


#' The StratifiedCEGPosition1 Function
#'
#' The hidden function StratifiedCegPosition() from the ceg package.
#'
#' @param stage A list detailing the stage structure of a Stratified.staged.tree model at each level.
#' @param num.category An integer vector detailing the number of categories for each variable.
#' @param num.situation An integer vector detailing the number of situations at each level.
#'
#' @return A list detailing the position structure at each level.
#'
StratifiedCegPosition1<-function (stage, num.category, num.situation)
{
  num.level <- length(num.category)
  result <- list()
  length(result) <- num.level
  result[[num.level]] <- PositionLevel1(stage[[num.level]],
                                        0, num.situation[num.level])
  for (level in (num.level - 1):2) {
    result[[level]] <- PositionLevel1(stage[[level]],
                                      num.category[level], num.situation[level + 1], result[[level + 1]])
  }
  result[[1]] <- list(list(1))
  return(result)
}

#' The Position Level Function
#'
#' The function [ceg::PositionLevel()]
#'
#' @param stage.list A list detailing the stage structure for a single level.
#' @param num.category An integer value detailing the number of stages at a single level.
#' @param num.situation.next An integer value detailing the number of situations at the next level.
#' @param pos.next.level A list detailing the position structure for the next level.
#'
#' @return A list detailing the position structure for a single level.
#'
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
#' The function [ceg::PositionVector()]
#'
#' @param num.situation An integer value detailing the number of situations at a single level.
#' @param pos.list A list detailing the position structure for a single level.
#'
#' @return An integer vector detailing the position structure for a single level.
#'
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
#' The function [ceg::PositionStage()]
#'
#' @param stage.vector An integer vector detailing the stages merged at a single level.
#' @param num.category An integer value for the number of categories for the associated variable.
#' @param pos.next.level An integer vector detailing the position structure for a single level.
#'
#' @return A list
#'
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
#' The hidden function PairwisePosition() from the ceg package.
#'
#' @param pair.situation An integer vector
#' @param num.category An integer value
#' @param pos.next.level An integer vector
#'
#' @return A logical value indicating whether two situations are in the same position.
#'
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


