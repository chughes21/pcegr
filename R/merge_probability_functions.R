
#' The Merge Probability function
#'
#' Calculate the probability of two true stages being merged
#'
#' Given a sample size, at least one stage probability, path probabilities and a prior structure, returns the probability of the two stages being merged. If only one probability is specified, the second probability is assumed to be the same. If alpha is provided, the priors are based on a uniform path prior and the stage probabilities.
#'
#' @param n An integer value for the total sample size.
#' @param p1 A numeric value between 0 and 1 specifying the assumed true probability for the first stage.
#' @param p2 A numeric value between 0 and 1 specifying the assumed true probability for the second stage. If NULL, p1 is used by default.
#' @param l1 A numeric value between 0 and 1 specifying the path probability up to the first stage.
#' @param l2 A numeric value between 0 and 1 specifying the path probability up to the second stage.
#' @param alpha A numeric value specifying the equivalent sample size used for setting the prior. If NULL, the inputted prior values are used.
#' @param a11 A numeric value specifying the first prior component for the first stage. If NULL, alpha is used for the prior instead.
#' @param a12 A numeric value specifying the second prior component for the first stage. If NULL, alpha is used for the prior instead.
#' @param a21 A numeric value specifying the first prior component for the second stage. If NULL, alpha is used for the prior instead.
#' @param a22 A numeric value specifying the second prior component for the second stage. If NULL, alpha is used for the prior instead.
#'
#' @return A list containing i) A matrix of the joint probabilities of each possible count combination; ii) A matrix of whether each count combination merges or not; iii) the total probability of merging, the sum of the elementwise product of (i) and (ii).
#' @export
#'
#' @examples
#' merge_prob(100,p1=0.4,p2=0.6,l1=1/8,l2=1/8,alpha=1)
#' merge_prob(100,p1=0.4,p2=0.6,l1=1/8,l2=1/8,a11=1/16,a12=1/16,a21=1/16,a22=1/16)
#' merge_prob(100,p1=0.4,p2=0.6,l1=1/8,l2=1/8,a11=1/20,a12=3/40,a21=3/40,a22=1/20)
merge_prob<-function(n, p1,p2=NULL, l1, l2, alpha=NULL, a11 = NULL, a12 = NULL, a21 = NULL, a22 = NULL){
  n1<-round(n*l1,0)
  n2<-round(n*l2,0)

  v1<-c(0:n1)
  v2<-c(0:n2)

  if(length(alpha)==0 & length(a11)==0){
    stop("At least one of alpha or prior must be specified")
  }

  if(max(length(a11),length(a12),length(a21),length(a22))>0 & min(length(a11),length(a12),length(a21),length(a22))==0){
    stop("Either supply no stage priors or all")
  }

  if(length(p2)==0){
    p2<-p1
  }

  if(max(p1,p2)>1){
    stop("Probabilities should be less than 1")
  }

  if(min(p1,p2)<0){
    stop("Probabilities should be greater than 0")
  }


  if(length(a11)==0){
    a11=alpha*l1*p1
    a12=alpha*l1*(1-p1)
    a21=alpha*l2*p2
    a22=alpha*l2*(1-p2)
  }

  x1<-dbinom(v1,n1,p1)
  x2<-dbinom(v2,n2,p2)

  M_prob<-outer(x1,x2)

  M_ahc<-matrix(nrow=n1+1,ncol=n2+1)

  a1<-c(a11,a12)
  a2<-c(a21,a22)

  for(i in 0:n1){
    for(j in 0:n2){
      y1<-c(i,n1-i)
      y2<-c(j,n2-j)

      M_ahc[i+1,j+1]<-ahc_merge(y1,y2,a1,a2)

    }
  }

  M_out<-M_ahc>0

  return(list(prob_all=M_prob,ahc_result=M_out,prob=sum(M_prob*M_out)))

}
