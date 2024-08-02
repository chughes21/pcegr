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
  mod$model.score<-mod$model.score+ahc_merge(d1,d2,a1,a2,FALSE)

  mod$posterior.expectation[[level]][ref1,]<-(mod$prior.distribution[[level]][ref1,]+mod$data.summary[[level]][ref1,])/sum(mod$prior.distribution[[level]][ref1,]+mod$data.summary[[level]][ref1,])
  mod$posterior.expectation[[level]][ref2,]<-NA

  }

  mod$merged<-cbind(mod$merged,c(sit1,sit2,level))

  mod$stage.structure[[level]][[ref1]]<-c(mod$stage.structure[[level]][[ref1]],mod$stage.structure[[level]][[ref2]])
  mod$stage.structure[[level]][[ref2]]<-NA

  #for result, we need the stage number (ordered)

  stages<-mod$stages
  #stages.level<-stages[which((stages>=start.situations[level])&(stages<=end.situations[level]))]
  #stages.pre.level<-stages[which(stages<start.situations[level])]
 # num.stages.pre.level<-length(stages.pre.level)

  true.sit1<-start.situations[level]+ref1-1
  true.sit2<-start.situations[level]+ref2-1

 # ind.sit1<-which(stages.level==true.sit1)
 # ind.sit2<-which(stages.level==true.sit2)

  ind.stage1<-which(stages==true.sit1)
  ind.stage2<-which(stages==true.sit2)

  mod$result[[ind.stage1]]<-cbind(mod$result[[ind.stage1]],mod$result[[ind.stage2]])
  mod$result[[ind.stage2]]<-NULL

  mod$stages<-stages[-ind.stage2]

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
#' @param background.order An integer vector detailing how the default backgrounds from the [[background_extractor()]] function have been reordered in the augmented dataset. If integer $m$ is in position $n$ of the vector, than the background $m$ from [[background_extractor()]] is now background $n$ in the augmented data set.
#' @param back.prior.input A list where each element is a matrix corresponding to a Dirichlet prior distribution for each level of the background tree before model selection. If NULL, the prior is extracted from the background model and assumed equally distributed.
#' @param resp.prior.input A list where each element is a matrix corresponding to a Dirichlet prior distribution for each level of the response tree before model selection. If NULL, the prior is extracted from the response model and assumed equally distributed.
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
model_combiner<-function(data,mod.background,mod.response,background.order=NULL, back.prior.input=NULL,resp.prior.input=NULL){

  resp.variable<-mod.response$event.tree$num.variable #how many variables in the response tree
  cutpoint.variable<-mod.background$event.tree$num.variable #how many variables in background tree
  num.variable<-cutpoint.variable+resp.variable-2 #total variables in the tree - subtract two because of background variavble and double counting intermediate var

  poisson.response<-mod.response$event.tree$poisson.response
  variable.time<-mod.response$event.tree$variable.time

  if(num.variable != dim(data)[2]-1*variable.time){
    stop("Number of variables in data set does not match models provided")
  }

  equivsize<-sum(mod.response$prior.distribution[[1]],na.rm = TRUE)

  #the priors should be taken from the given models-TO DO

  if(poisson.response){
    gamma_alpha<-sum(mod.response$prior.distribution[[resp.variable]][,1],na.rm = TRUE)/length(mod.response$prior.distribution[[resp.variable]][,1])
    gamma_beta<-sum(mod.response$prior.distribution[[resp.variable]][,2],na.rm = TRUE)/length(mod.response$prior.distribution[[resp.variable]][,2])
  }else{
    gamma_alpha<-1
    gamma_beta<-2
  }

  #compute a saturated model without account for cutpoint.variable
  mod.sat<-saturated_model_computer(data,cutpoint.variable=NULL,equivsize, poisson.response,variable.time,gamma_alpha,gamma_beta,prior_input = NULL)

  back.merge<-mod.background$merged
  resp.merge<-mod.response$merged #the situations merged in the response ceg

  data.resp<-background_extractor(data,mod.background)

  if(length(background.order)>0){
    n_group<-nlevels(data.resp[,1])
    if(length(unique(background.order))!=n_group){
      stop("Background order input should have same number of unique elements as background groups")
    }else if(max(background.order)>n_group){
      stop("Maximum value in background order vector should be less than or equal to number of groups")
    }else if(min(background.order)<1){
      stop("Minimum value in background order vector should be 1")
    }else{
      data.temp<-data.resp
      for(i in 1:n_group){
        ind.temp<-which(data.resp[,1]==i)
        data.temp[ind.temp,1]<-rep(which(background.order==i),length(ind.temp))
      }
    }

    data.resp<-data.temp

  }

  if(length(back.prior.input)>0){
    back.prior<-back.prior.input
    back.prior.check<-TRUE
    if(length(back.prior.input)!=length(mod.background$prior.distribution)){
      stop("Length of prior input doesn't match background model")
    }
  }else{
    back.prior<-mod.background$prior.distribution
    back.prior.check<-FALSE
  }


  if(length(resp.prior.input)>0){
    resp.prior<-resp.prior.input
    resp.prior.check<-TRUE
    if(length(resp.prior.input)!=length(mod.response$prior.distribution)){
      stop("Length of prior input doesn't match response model")
    }
  }else{
    resp.prior<-mod.response$prior.distribution
    resp.prior.check<-FALSE
    }


  #create a saturated tree
  #first need levels
  cats.sat<-lapply(data,levels)
  cats.resp<-lapply(data.resp,levels)

  if(variable.time){
    ind.excl.sat<-c(num.variable:(num.variable+1))
    ind.excl.resp<-c(resp.variable:(resp.variable+1))
  }else{
    ind.excl.sat<-num.variable
    ind.excl.resp<-resp.variable
  }

  tree.sat<-rev(expand.grid(rev(cats.sat[-ind.excl.sat]))) #don't need the last column as its leaves and also wanna exclude variable time
  tree.resp<-rev(expand.grid(rev(cats.resp[-ind.excl.resp])))

  #now create a tree that is the saturated tree except with background imputed
  tree.background<-tree.sat
  tree.background<-pcegr:::cat_replacer(tree.sat[,1:(cutpoint.variable-1)],pcegr:::cat_extractor(mod.background))

  tree.ref<-cbind(background=tree.background[,1],tree.sat[-(1:(cutpoint.variable-1))])

  if(length(background.order)>0){
    tree.temp<-tree.ref
    for(i in 1:n_group){
      ind.temp<-which(tree.ref[,1]==i)
      tree.temp[ind.temp,1]<-rep(which(background.order==i),length(ind.temp))
    }
    tree.ref<-tree.temp
  }

  #don't know if i need the following two lines
  situations.sat<-mod.sat$event.tree$num.situation
  start.situations.sat<-c(1,cumsum(situations.sat[1:(num.variable-1)])+1)

  back.stage<-mod.background$stages #the stages in the background

  #now see what levels the stages are at for the response model

  resp.stages<-mod.response$stages
  situations.resp<-mod.response$event.tree$num.situation
  start.situations.resp<-c(1,cumsum(situations.resp[1:(resp.variable-1)])+1)

  resp.stages<-rbind(resp.stages,sapply(resp.stages,findInterval,vec=start.situations.resp))

  #index isn't working - should start at num.variable

  for(i in num.variable:2){

    lev<-length(cats.sat[[i]])

    if(i>cutpoint.variable){
      i.resp<-i-cutpoint.variable+2
      situations<-c(1:length(tree.resp[,i.resp-1]))

      corr.sit<-numeric(length(situations)) #a vector showing the corresponding first situation in the big tree for the response tree
      count.sit<-0 #a counter checking the number of situations

      prior.sum<-colSums(resp.prior[[i.resp]],na.rm=TRUE) #the total column sums of the prior to be distributed equally

      for(j in situations){
        v<-tree.resp[j,] #the covariates
        ind.sit<-which(prodlim::row.match(tree.ref,v)==1) #which situations match those covariates
        corr.sit[j]<-ind.sit[1] #keep track of which elements situations correspond to
        p<-length(ind.sit)
        count.sit<-count.sit+p #keep track of total number of situations
        ref1<-ind.sit[1]

        if(resp.prior.check){
          mod.sat$prior.distribution[[i]][ref1,]<-resp.prior[[i.resp]][j,]/p #divide by p to distributed across all paths
        }else{
          prior.const<-prior.sum/length(situations)
          prior.frac<-prior.const/p #divide by the number of paths - assumes each path had the same prior before
          mod.sat$prior.distribution[[i]][ref1,]<-prior.frac
        }

        if(p>1){
          for(k in 2:p){
            ref2<-ind.sit[k]
            if(resp.prior.check){
              mod.sat$prior.distribution[[i]][ref2,]<-resp.prior[[i.resp]][j,]/p
            }else{
              mod.sat$prior.distribution[[i]][ref2,]<-prior.frac
            }
            #Updating the priors here should be enough

            # print(c(j,k))
            mod.sat<-stage_updater(mod.sat,i,ref1,ref2)
            #print(paste(j,ref1,ref2,sum(mod.sat$data.summary[[7]],na.rm=TRUE)))
            # print(mod.sat$model.score)
          }
        }
      }

      if(count.sit!=situations.sat[i]){
        stop(paste0("Not all situations accounted for at level ",i))
      }

      stage.ref<-resp.stages[1,which(resp.stages[2,]==i.resp)] #which stages are at this level
      stage.ref<-stage.ref- start.situations.resp[i.resp]+1 #adjust for indexing

      for(j in stage.ref){
        ind.stage<-mod.response$stage.structure[[i.resp]][[j]]
        n.stage<-length(ind.stage)
        ref1<-corr.sit[ind.stage[1]]

        if(n.stage>1){

          for(k in 2:n.stage){
            ref2<-corr.sit[ind.stage[k]]
            #print(paste(j,mod.sat$data.summary[[i]][ref1,],mod.sat$data.summary[[i]][ref2,]))
            mod.sat<-stage_updater(mod.sat,i,ref1,ref2)
          }
        }
      }

    }else{

      situations<-c(1:situations.sat[i])
      ind.stage.level<-which(back.stage>=start.situations.sat[i]) #index of which stages are at this level
      stage.ref<-back.stage[ind.stage.level]
      stage.ref<-stage.ref-start.situations.sat[i]+1

      if(back.prior.check){
        mod.sat$prior.distribution[[i]]<-back.prior[[i]]
      }else{
        prior.sum<-colSums(back.prior[[i]],na.rm=TRUE)
        prior.frac<-prior.sum/situations.sat[i]
        mod.sat$prior.distribution[[i]]<-matrix(rep(prior.frac,situations.sat[i]),ncol=lev,nrow=situations.sat[i],byrow=TRUE)
      }

      for(j in stage.ref){
        ind.stage<-mod.background$stage.structure[[i]][[j]]
        n.stage<-length(ind.stage)

        if(n.stage>1){
          ref1<-ind.stage[1]
          for(k in 2:n.stage){
            ref2<-ind.stage[k]
            mod.sat<-stage_updater(mod.sat,i,ref1,ref2)
          }
        }
      }
      back.stage<-back.stage[-ind.stage.level]
    }

    cut.ind<-seq(from=1,to=count.sit-1,by=lev) #a sequence to take only one path from every combo
    cut.ind.resp<-seq(from=1,to=length(situations)-1,by=lev) #a sequence to take only one path from every combo for response tree
    tree.ref<-tree.ref[cut.ind,]
    tree.resp<-tree.resp[cut.ind.resp,]

    #   situations.lev<-((start.situations.sat[i]+situations.sat[i]-1):start.situations.sat[i])
    #    for(k in situations.lev){
    #    if(is.na(mod.sat$result[[]]))
    #  }

  }

  #To do in stage_updater
  #1. Fix score output - doesn't work in current way


  return(mod.sat)

}
