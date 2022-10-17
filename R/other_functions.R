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
#' @return An integer value specifying how how many elements of the vector are zero.
#'
#' @examples zero_checker(c(0,2,0,1,1))
zero_checker<-function(v){
  return(sum(v==0))
}

#' The Nonzero Checker
#'
#' A function that checks how many elements of a given vector are nonzero.
#'
#' @param v A numeric vector
#'
#' @return An integer value specifying how how many elements of the vector are zero.
#'
#' @examples
#' nonzero_checker(c(0,2,0,1,1))
nonzero_checker<-function(v){
  return(sum(abs(v)>0))
}

#' The Keeper Function
#'
#' A function that, when given a specific value, keeps all values in a vector that are equal as 1, and any value that is different as 0.
#'
#' @param v A numeric vector
#' @param n A numeric value
#'
#' @return A vector of the same length where values equal to n are 1, and values different to n are 0
#'
#' @examples
#' keeper(c(0,2,0,1,1),1)
keeper<-function(v,n){
  output<-v
  m<-length(output)
  output[which(output != n)]=0
  output[which(output == n)]=1
  return(output)
}

#' The vec_from function
#'
#' A function to create a descending vector of natural numbers
#'
#' @param x A positive integer
#'
#' @return A vector of integers from $x-1$ to 0.
#'
#' @examples
#' vec_from(7)
vec_from<-function(x){
  return(c((x-1):0))
}

#' The vec_to function
#'
#' A function to create an ascending vector of natural numbers
#'
#' @param x A positive integer
#'
#' @return A vector of integers from 0 to $x-1$.
#'
#' @examples
#' vec_to(7)
vec_to<-function(x){
  return(c(0:(x-1)))
}

