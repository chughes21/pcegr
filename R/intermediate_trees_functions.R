#' Category Extractor Function
#'
#' Looks at each possible category combination in a StagedTree model, and returns what stage they are in for the final level.
#'
#' Intended to be used for Intermediate CEGs, this function looks at each path in the StagedTree model up to the final level, and records what stage it is in for this level.This will define the categories for an Intermediate CEG.
#'
#' @param mod A StagedTree object.
#'
#' @return A data frame where the last column corresponds to the stage for each path.
#'
#' @examples
#' mod<-pceg(knee_pain_obs[,1:3],2,FALSE,FALSE)
#' cat_extractor(mod)
cat_extractor<-function(mod){
  stages<-mod$result
  n_stages<-length(stages)
  posterior<-mod$posterior.expectation
  n_level<-length(posterior)
  n_final_stage<-sum(!is.na(posterior[[n_level]][,1]))
  stage_seq<-n_stages-c((n_final_stage-1):0)

  stages<-stages[stage_seq]
  n_cov<-n_level-1

  situations<-sapply(stages,length)/n_cov
  n_situations<-sum(situations)

  data<-data.frame(matrix(nrow=n_situations,ncol=n_cov),category=rep(1,n_situations))

  counter<-1

  for(i in 1:n_final_stage){
    for(j in 1:situations[i]){

      if(situations[i]==1){
        stages_temp<-matrix(stages[[i]],ncol=1)
      }else{
        stages_temp<-stages[[i]]
      }

      data[counter,1:n_cov]<-t(stages_temp[,j])
      data[counter,n_cov+1]<-i

      counter<-counter+1

    }
  }

  if(dim(data)[2]==2){
    colnames(data)[1]<-"x1"
  }

  return(data)

}

#' Background Imputer Function
#'
#' Imputes a new categorisation of the background variables.
#'
#' Intended to be used for Intermediate CEGs, this function takes a data set and a definition of new categorisations of certain background variables in the data. This categorisation can be obtained from the [[cat_extractor()]] function.
#'
#' @param data A data set, where the observed count vector and time vector (if included) are the last two columns.
#' @param background.data A data set where each column besides the last corresponds to the same column as
#'
#' @return A data set where the variables in background.data are replaced by their corresponding category for a new variable "background", and the rest of the variables in data are unchanged.  of
#'
#' @examples
#  mod<-pceg(knee_pain_obs[,1:3],2,FALSE,FALSE)
#' backdata<-cat_extractor(mod)
#' background_imputer(knee_pain_obs, backdata)
background_imputer<-function(data,background.data){

  N<-dim(data)[1]
  n<-dim(data)[2]

  cutpoint<-dim(background.data)[2]

  if(cutpoint==n){
    warning("Data and background.data have same dimension - this corresponds to a recategorisation of the entire data set")
  }

  data.pre<-data[,1:(cutpoint-1)]
  data.post<-data[,cutpoint:n]

  background<-factor(rep(1,N),levels=c(1:max(background.data[,cutpoint])))

  which_not_1<-which(background.data[,cutpoint]>1)

  for(j in which_not_1){
    v<-background.data[j,-cutpoint]
    if(length(v)==1){
      ind<-which(data.pre==v)
    }else{
      ind<-which(row.match(data.pre,v)==1)
    }
    background[ind]<-background.data[j,cutpoint]

  }

  data.out<-data.frame(cbind(background=background,data.post))
  return(data.out)

}

#' Background Extractor
#'
#' Takes a background model which produces a new categorisation of certain variables, and replaces them in data.
#'
#' Intended for use with Intermediate CEGs, this function takes a data set and background model which has been fit on a subset of the data in the same order, and replaces the background variables with their categorisation from the background model.
#'
#' @param data A data set, where the observed count vector and time vector (if included) are the last two columns.
#' @param mod.background A StagedTree object fit on a subset of data, from the first variable to some cutpoint variable.
#'
#' @return A data set where the variables from the background model, besides the cutpoint variable, have been replaced based on the stage structure for the cutpoint variable.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs[,1:3],2,FALSE,FALSE)
#' background_extractor(knee_pain_obs,mod)
background_extractor<-function(data,mod.background){
  background.data<-cat_extractor(mod.background)
  return(background_imputer(data,background.data))
}

#' Category Replacer
#'
#' Replace the categories in X with the corresponding final column in Y
#'
#' @param X A data set.
#' @param Y A data set with one more column than X.
#'
#' @return A data set where the categories in X are placed by the final column in Y.
#'
cat_replacer<-function(X,Y){

  if(is.null(dim(X))){
    X<-matrix(X)
  }

  if(dim(X)[2]!=dim(Y)[2]-1){
    stop("Y must one more column than X")
  }

  nrow<-dim(X)[1]
  ncol<-dim(X)[2]

  output<-matrix(nrow=nrow,ncol=ncol)

  for(i in 1:nrow){
    v<-X[i,]
    if(ncol>1){
      ind<-which(row.match(Y[,-(ncol+1)],v)==1)
      output[i,]<-rep(Y[ind,ncol+1],ncol)
    }else{
      ind<-which(v==Y[,1])
      output[i,]<-rep(Y[ind,ncol+1],ncol)
    }
  }

  return(output)

}

#' Saturated Model Computer
#'
#' A function to quickly compute a saturated StagedTree model from a data set, with a cutpoint variable specifying where the saturation begins.
#'
#' @param data A data set, where the observed count vector and time vector (if included) are the last two columns.
#' @param cutpoint.variable An integer value specifying after which variable the tree should be saturated. If NULL, the entire tree will be saturated.
#' @param equivsize A numeric value specifying the equivalent sample size for the prior, a measure of confidence in the prior.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param gamma_alpha A numeric vector for the shape hyperparameters of the Gamma priors for the Poisson rates, if applicable. If a single value, each prior takes the same value.
#' @param gamma_beta A numeric vector for the rate hyperparameters of the Gamma priors for the Poisson rates, if applicable. If a single value, each prior takes the same value.
#' @param prior_input A list where each element is a matrix corresponding to a Dirichlet prior distribution for each level of the tree.
#'
#' @return A StagedTree model where the variables up to the cutpoint variable are merged normally, and any subsequent variables are saturated.
#' @export
#'
#' @examples
#' saturated_model_computer(knee_pain_obs,3,2,TRUE,TRUE)
saturated_model_computer<-function(data, cutpoint.variable=NULL, equivsize=2,poisson_response = FALSE, variable_time = FALSE, gamma_alpha = 1, gamma_beta = 2,prior_input=NULL){

  if(is.null(cutpoint.variable)){
    cutpoint.variable<-1
  }

  num.variable<-dim(data)[2]-1*variable_time

  mod.sat<-pceg(data,equivsize,poisson_response,variable_time,gamma_alpha=gamma_alpha,gamma_beta = gamma_beta,saturated=c((cutpoint.variable+1):num.variable),prior_input = prior_input)

  return(mod.sat)

}

#Doesn't update the result yet - should

#' Stage Updater
#'
#' Updates a model by merging two stages.
#'
#' @param mod A StagedTree object.
#' @param level An integer value specifying what level of the tree the merging takes place at.
#' @param ref1 An integer value specifying what the first stage at the chosen level to be merged is.
#' @param ref2 An integer value specifying what the second stage at the chosen level to be merged is.
#'
#' @return A StagedTree object which has merged the two stages.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2);plot(mod)
#' mod2<-stage_updater(mod,2,1,2);plot(mod2)
#' mod3<-stage_updater(mod,3,1,3);plot(mod3)
stage_updater<-function(mod,level,ref1,ref2){

  num.variable<-mod$event.tree$num.variable
  poisson.response<-mod$event.tree$poisson.response

  if(level>num.variable | level<1){
    stop("Please input a level between 1 and the number of variables")
  }

  if(level==num.variable & poisson.response){
    poisson<-TRUE
  }else{
    poisson<-FALSE
  }

  num.situation<-mod$event.tree$num.situation
  start.situations<-c(1,cumsum(num.situation[1:(num.variable-1)])+1)
  end.situations<-cumsum(num.situation)

  situations<-c(start.situations[level]:end.situations[level])

  sit1<-situations[ref1]
  sit2<-situations[ref2]

  d1<-mod$data.summary[[level]][ref1,]
  d2<-mod$data.summary[[level]][ref2,]

  a1<-mod$prior.distribution[[level]][ref1,]
  a2<-mod$prior.distribution[[level]][ref2,]

  mod$data.summary[[level]][ref1,]<-mod$data.summary[[level]][ref1,]+mod$data.summary[[level]][ref2,]
  mod$prior.distribution[[level]][ref1,]<-mod$prior.distribution[[level]][ref1,]+mod$prior.distribution[[level]][ref2,]

  mod$data.summary[[level]][ref2,]<-NA
  mod$prior.distribution[[level]][ref2,]<-NA

  if(poisson){
  data.temp<-list(d1,d2)
  prior.temp<-list(a1,a2)

  mod$model.score<-mod$model.score+as.numeric(bayes_factor(data.temp,prior.temp,1,2))

  mod$posterior.expectation[[level]][ref1]<-(mod$prior.distribution[[level]][ref1,1]+mod$data.summary[[level]][ref1,1])/(mod$prior.distribution[[level]][ref1,2]+mod$data.summary[[level]][ref1,2])
  mod$posterior.expectation[[level]][ref2]<-NA

  }else{
  mod$model.score<-mod$model.score+ahc_merge(d1,d2,a1,a2)

  mod$posterior.expectation[[level]][ref1,]<-(mod$prior.distribution[[level]][ref1,]+mod$data.summary[[level]][ref1,])/sum(mod$prior.distribution[[level]][ref1,]+mod$data.summary[[level]][ref1,])
  mod$posterior.expectation[[level]][ref2,]<-NA

  }

  mod$merged<-cbind(mod$merged,c(sit1,sit2,level))

  mod$stage.structure[[level]][[ref1]]<-c(mod$stage.structure[[level]][[ref1]],mod$stage.structure[[level]][[ref2]])
  mod$stage.structure[[level]][[ref2]]<-NA

  return(mod)
}


#' Model Combiner
#'
#' Takes a background model and response model and combines them in one full model.
#'
#' Intended for use with Intermediate CEGs, this creates a StagedTree object which is equivalent to combining the background process in mod.background, with the response process in mod. response. This combined model can be used to compare to other methods of model selection, such as the score.
#'
#' @param data A data set, where the observed count vector and time vector (if included) are the last two columns.
#' @param mod.background A StagedTree object modelling the background process.
#' @param mod.response A StagedTree object modelling the response process.
#' @param prior_input A list where each element is a matrix corresponding to a Dirichlet prior distribution for each level of the tree. If NULL, the default prior is used.
#'
#' @return A StagedTree object which is equivalent to combining mod.background and mod.response.
#' @export
#'
#' @examples
#' mod1<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' mod.back<-pceg(knee_pain_obs[,1:3],2,FALSE,FALSE)
#' newdata<-background_extractor(knee_pain_obs,mod.back)
#' mod.resp<-pceg(newdata,2,TRUE,TRUE)
#' mod2<-model_combiner(knee_pain_obs,mod.back,mod.resp)
model_combiner<-function(data,mod.background,mod.response,prior.input=NULL){

  resp.variable<-mod.response$event.tree$num.variable
  cutpoint.variable<-mod.background$event.tree$num.variable
  num.variable<-cutpoint.variable+resp.variable-2

  poisson.response<-mod.response$event.tree$poisson.response
  variable.time<-mod.response$event.tree$variable.time

  if(num.variable != dim(data)[2]-1*variable.time){
    stop("Number of variables in data set does not match models provided")
  }

  equivsize<-sum(mod.response$prior.distribution[[1]],na.rm = TRUE)

  if(poisson.response){
    gamma_alpha<-sum(mod.response$prior.distribution[[resp.variable]][,1],na.rm = TRUE)/length(mod.response$prior.distribution[[resp.variable]][,1])
    gamma_beta<-sum(mod.response$prior.distribution[[resp.variable]][,2],na.rm = TRUE)/length(mod.response$prior.distribution[[resp.variable]][,2])
  }else{
    gamma_alpha<-1
    gamma_beta<-2
  }

  mod.sat<-saturated_model_computer(data,cutpoint.variable,equivsize, poisson.response,variable.time,gamma_alpha,gamma_beta,prior_input = prior.input)

  resp.merge<-mod.response$merged

  situations<-mod.sat$event.tree$num.situation
  start.situations<-c(1,cumsum(situations[1:(num.variable-1)])+1)
  end.situations<-cumsum(situations)

  resp.situations<-mod.response$event.tree$num.situation
  resp.start.situations<-c(1,cumsum(resp.situations[1:(mod.response$event.tree$num.variable-1)])+1)
  resp.end.situations<-cumsum(resp.situations)

  data.resp<-background_extractor(data,mod.background)

  cats.sat<-lapply(data,levels)
  cats.resp<-lapply(data.resp,levels)

  #need to fix, doesn't work when more than one covariate in background

  situations.list<-list()
  for(i in 1:num.variable){

    #resp.situations.temp<-c(resp.start.situations[i]:resp.end.situations[i])

    situations.list[[i]]<-cbind(c(start.situations[i]:end.situations[i]),NA)
    if(i <= cutpoint.variable){
      situations.list[[i]][,2]<-situations.list[[i]][,1]
    }

    if(i>cutpoint.variable){

      i.resp<-i-cutpoint.variable+2
      resp.situations.temp<-c(resp.start.situations[i.resp]:resp.end.situations[i.resp])

      tree.sat<-rev(expand.grid(rev(cats.sat[-(i:(num.variable+1*variable.time))])))
      tree.resp<-rev(expand.grid(rev(cats.resp[-(i.resp:(resp.variable+1*variable.time))])))

      tree.sat[,1:(cutpoint.variable-1)]<-cat_replacer(tree.sat[,1:(cutpoint.variable-1)],cat_extractor(mod.background))

      if(cutpoint.variable>2){
        #tree.resp<-cbind(rep(tree.resp$background,cutpoint.variable-2),tree.resp) old version
        tree.resp<-cbind(matrix(rep(tree.resp$background,cutpoint.variable-2),ncol=cutpoint.variable-2,byrow=FALSE),tree.resp)
      }

      for(j in 1:situations[i]){
        v<-tree.sat[j,]
        ind<-which(row.match(tree.resp,v)==1)
        situations.list[[i]][j,2]<-resp.situations.temp[ind]
      }

      #first merge all equal situations together

      for(j in 1:resp.situations[i.resp]){

        num<-resp.situations.temp[j]

        ind.num<-which(situations.list[[i]][,2]==num)
        p<-length(ind.num)
        ref1<-ind.num[1]
        if(p>1){
          for(k in 2:p){
            ref2<-ind.num[k]

            mod.sat<-stage_updater(mod.sat,i,ref1,ref2)

            situations.list[[i]][ref2,]<-rep(NA,2)

          }
        }
      }
    }
  }

  for(i in 2:(resp.variable-1)){

    #now perform the actual merging

    resp.merge.temp<-resp.merge[1:2,which(resp.merge[3,]==i)]

    if(length(resp.merge.temp)>0){

    variable<-cutpoint.variable+i-1

    situations.temp<-situations.list[[variable]]

    for(j in 1:dim(resp.merge.temp)[2]){

      merge.top<-resp.merge.temp[1,j]
      merge.bottom<-resp.merge.temp[2,j]

      ref1<-which(situations.temp[,2]==merge.top)
      ref2<-which(situations.temp[,2]==merge.bottom)

      mod.sat<-stage_updater(mod.sat,variable,ref1,ref2)

      situations.list[[variable]][ref2,]<-rep(NA,2)

      #fix result also
    }
   }
  }

  return(mod.sat)

}
