#' The Binary Resizer Function
#'
#' Reparameterise any data set into one where each non-binary variable is replaced by dummy variables to account for the number of categories.
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#'
#' @return A data set where each variable is binary.
#' @export
#'
#' @examples
#' v1<-factor(sample(c(0,1,2),10000,TRUE))
#' data<-cbind(v1,knee_pain_obs)
#' summary(data)
#' summary(binary_resizer(data))
binary_resizer<-function(data,poisson_response = TRUE, variable_time = TRUE){

  if(!poisson_response & variable_time){
    stop("Variable time requires a Poisson response")
  }

  n<-dim(data)[2]-1*poisson_response-1*variable_time
  p<-dim(data)[1]

  numbcat<-sapply(data[,1:n],FUN=nlevels)

  if(max(numbcat)==2){
    stop("Dataset is already binary")
  }

  excess<-numbcat-rep(2,n)
  cum_excess<-cumsum(excess)
  bin_n<-n+sum(excess)

  ind<-which(excess>0)
  bin_ind<-which(excess==0)
  bin_ind_new<-bin_ind+cum_excess[bin_ind]

  col_names<-colnames(data[,1:n])

  new_col_names<-character(bin_n)
  new_col_names[bin_ind_new]<-col_names[bin_ind]

  for(i in 1:length(ind)){

    var<-ind[i]

    if(i==1){
    k<-var
    }else{
      k<-k+var-ind[i-1]
    }


    lev<-levels(data[,var])

    temp_ind<-which(data[,var]==lev[1])

    temp<-factor(rep("other",p),levels=c(lev[1],"other"))
    temp[temp_ind]<-lev[1]

    if(var == 1){
      newdata<-data.frame(temp)
    }else if(i == 1){
      newdata<-data.frame(data[,(1:(var-1))],temp)
    }else if(i == length(ind)){
      if(var-1>=ind[i-1]+1){
      newdata<-cbind(newdata,data[,c((ind[i-1]+1):(var-1))],temp)
      }else{
        newdata<-cbind(newdata,temp)
      }
    #}else if(var+1 == ind[i+1]){
     # newdata<-cbind(newdata,temp) #don't know if this right or necessary
    }else{
      newdata<-cbind(newdata,data[,((ind[i-1]+1):(var-1))],temp)
    }

    new_col_names[k]<-paste0(col_names[var],"-",lev[1])

    for(j in 1:excess[var]){

      k<-k+1

      if(j==1){
      temp_ind_NA<-temp_ind
      }else{
      temp_ind_NA<-unique(c(temp_ind_NA,temp_ind))
      }

      temp_ind<-which(data[,var]==lev[j+1])

      if(j==excess[var]){
      temp<-factor(rep(lev[j+2],p),levels=c(lev[(j+1):(j+2)]))
      }else{
      temp<-factor(rep("other",p),levels=c(lev[j+1],"other"))
      }

      temp[temp_ind]<-lev[j+1]
      temp[temp_ind_NA]<-NA

      newdata<-cbind(newdata,temp)

      new_col_names[k]<-paste0(col_names[var],"-",lev[j+1])

    }

  }

  if(max(ind)!=n){
    newdata<-cbind(newdata,data[,((max(ind)+1):n)])
  }

  colnames(newdata)<-new_col_names

  if(poisson_response){
    if(variable_time){
      newdata<-cbind(newdata,data[,(n+1):(n+2)])
      colnames(newdata)[(n+1):(n+2)]<-colnames(data)[(n+1):(n+2)]
    }else{
      newdata<-cbind(newdata,data[,n+1])
      colnames(newdata)[n+1]<-colnames(data)[n+1]
    }
  }

  return(newdata)

}
