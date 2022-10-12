#' The time minimiser
#'
#' A function that enforces a minimum value of time for an observed time.
#'
#' @param data A data set where the last column is the observed time vector.
#' @param t_min A numeric minimum time value.
#'
#' @return A data set where the observed times now have a minimum value.
#' @export
#'
#' @examples
#' n<-10
#' t<-runif(n)
#' y<-rpois(n,2*t)
#' data<-data.frame(y=y,t=t)
#' data
#' data.min<-time_minimiser(data,t_min=0.5)
#' data.min
time_minimiser<-function(data,t_min){
  n<-dim(data)[2]
  t<-data[,n]
  ind<-which(t<t_min)
  data[ind,n]<-t_min
  return(data)
}

#' The Zero Checker
#'
#' A function that checks how many elements of a given vector are zero.
#'
#' @param v A numeric vector
#'
#' @return An integer how how many elements of the vector are zero.
#'
#' @examples zero_checker(c(0,2,0,1,1))
zero_checker<-function(v){
  return(sum(v==0))
}
