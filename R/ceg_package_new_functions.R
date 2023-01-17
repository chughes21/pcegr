#' The Stage Structure Function
#'
#' This function takes a list of inputs and creates a stage structure that is compatible with the StagedTree S3 class.
#'
#' @param mod A list created in [pceg()].
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#' @param remove_risk_free A logical value indicating whether the risk free leaves and edges should be removed (TRUE) or not (FALSE).
#'
#' @return A list specifying the stage structure calculated in [pceg()] in a way that is compatible with the StagedTree S3 class.
#'
stage_structure<-function(mod,zip=FALSE,remove_risk_free = FALSE){
  output<-list()
  M<-mod$merged
  comparisonset<-mod$comparisonset
  numb<-mod$numb

  for(i in 1:length(numb)){

    N<-numb[[i]]

    if(length(M)==0){
      for(j in 1:N){
        cluster[[j]]=j
      }
    }else{

    m<-as.matrix(M[,which(M[3,]==(i-1))])

    if(i>1){
      comp<-comparisonset[[i-1]]
      start<-sum(numb[1:(i-1)])
      comp<-comp-start
    }else{
      start<-1
      comp<-1
    }

    m[1:2,]<-m[1:2,]-start

    if(length(m)>0){
      check<-pcegr:::merged_list_extractor(m)
    }else{
      check<-NA
    }

    final_ind <-(i==length(numb))
    zip_ind<-final_ind & zip

    if(final_ind & remove_risk_free){
      no_risk_ind<-NA
      m[1:2,]<-m[1:2,]/2
      comp<-comp[-1]/2
      N<-N/2
      for(j in 1:length(check)){
        check[[j]]<-check[[j]]/2
      }
    }else if(zip_ind & (!remove_risk_free)){
      no_risk_ind<-seq(from=1,to=numb[i]-1,by=2)
    }else{
      no_risk_ind<-NA
    }

    k<-1

    cluster<-vector(mode="list",length=N)

    for(j in 1:N){

      if((j==1)&zip_ind&(!remove_risk_free) ){
        cluster[[j]]<-no_risk_ind
      } else if((j %in% no_risk_ind[-1])&zip_ind){
        cluster[[j]]<-NA
      }else if(all(is.na(check))){
        cluster[[j]]<-j
      }else if(!(j %in% comp)){
        cluster[[j]]<-NA
      }else if(!(j%in%m[1:2,])){
        cluster[[j]]<-j
      }else{
        cluster[[j]]<-as.vector(unlist(check[k]))
        k<-k+1
      }
     }
    }
    output[[i]]=cluster

  }

  return(output)

}

#' The Output List Converter Function
#'
#' A function to transform the list of prior and data outputs into the S3 class versions of the prior and posterior.
#'
#' @param mod A list created in [pceg()].
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#'
#' @return A list for the prior distribution and a list for the posterior distribution.
output_list_converter<-function(mod,poisson_response=TRUE,levels = NA){
  prior<-mod$prior
  data<-mod$data
  numb<-mod$numb
  n<-length(numb)

  numb_temp<-c(1,numb)

  nodes_start<-cumsum(numb_temp)[1:n]
  nodes_end<-cumsum(numb)

  nodes<-nodes_end-nodes_start + 1

  prior_out<-vector(mode="list",length=n)
  data_out<-vector(mode="list",length=n)
  post_out<-vector(mode="list",length=n)

  for(i in 1:n){
    prior_mat_temp<-matrix(nrow=nodes[i],ncol=dim(prior[[nodes_start[i]]])[2])
    data_mat_temp<-matrix(nrow=nodes[i],ncol=dim(prior[[nodes_start[i]]])[2])

    nodes_temp<-seq(nodes_start[i],nodes_end[i],by=1)
    for(j in 1:nodes[i]){
      prior_mat_temp[j,]<-prior[[nodes_temp[j]]]
      data_mat_temp[j,]<-data[[nodes_temp[j]]]
    }

    final_poiss_ind<-(i==n)&poisson_response

    if(final_poiss_ind){
    colnames(prior_mat_temp)<-c("a","b")

    colnames(data_mat_temp)<-c(levels[[n]],"t")
    prior_out[[i]]<-prior_mat_temp
    data_out[[i]]<-data_mat_temp

    sum_out<-prior_mat_temp+data_mat_temp

    post_out[[i]]<-sum_out[,1]/sum_out[,2]

    }else{
    colnames(prior_mat_temp)<-levels[[i]]
    colnames(data_mat_temp)<-levels[[i]]

    prior_out[[i]]<-prior_mat_temp
    data_out[[i]]<-data_mat_temp

    sum_out<-prior_mat_temp+data_mat_temp

    post_out[[i]]<-sum_out/rowSums(sum_out)
    }

    }
  return(list(prior=prior_out,data=data_out,posterior=post_out))
}


