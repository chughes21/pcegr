
#' The Marginal Effect Function
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param mod A StagedTree model fit to the data set, as produced by pceg() or zipceg()
#' @param input_variable An integer vector detailing the covariates whose marginal effect is to be analysed. The default is to analyse all covariates.
#' @param rel_output A non-positive integer value detailing the variable which is being affected by the input variable(s), relative to the response. When 0, this will analyse the response variable, but any variable that appears after the input variables in the tree can be analysed.
#' @param max_per_plot An integer value specifying the maximum number of leaves that can be shown in a single plot.
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param zip A logical value indicating whether the model specified is zero-inflated (TRUE) or not (FALSE).
#'
#' @return A marginal effect plot for each variable.
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' marginal_effect(knee_pain_obs,mod)
marginal_effect<-function(data,mod,input_variable = c(),rel_output=0,max_per_plot = 4,variable_time = TRUE,zip=FALSE){

  names<-colnames(data)
  n<-dim(data)[2]

  if(zip){
    names<-c(names[c(1:(n-1-1*variable_time))],"risk free",names[-c(1:(n-1-1*variable_time))])
  }

  poisson_response<-mod$event.tree$poisson.response
  remove_risk_free<-mod$remove.risk.free.edges

  numbvariables<-n + 1*zip-1*variable_time #total number of variables in the dataset including response and possibly ZIP, without time
  numbcovar<-numbvariables-1-1*zip #total number of covariates in the data

  if(rel_output>0){
    stop("Output variable relative to the response should be nonpositive") #default output variable is the response
  }

  if(!poisson_response & (variable_time)){
    stop("Poisson response needed for variable time")
  }

  output_variable <- numbvariables + rel_output

  if(output_variable <= numbcovar){ #if the output variable is one of the covariates
    poisson_response<-FALSE
    variable_time<-FALSE
    remove_risk_free<-FALSE
    zip<-FALSE
  }

  resp_out<-(rel_output==0) #is the response the output variable

  if(length(input_variable)==0){
    input_variable <- c(1:(output_variable-1-1*zip*resp_out)) #default input variables are all variables not including the risk state if ZIP
  }

  if(length(output_variable)>1){
    stop("Only one output variable possible")
  }

  if(max(input_variable)>numbvariables){
    stop("Input variable number exceeds variables in data")
  }

  if(max(input_variable)>= output_variable){
    stop("Output variable must come after input variable(s)")
  }

  cv<-output_variable-1-1*zip*resp_out #the number of covariates that will be used for this
  #don't include the risk state as a covariate cause it doesn't affect the parameter estimates differently (always 0 or not)
  cv_ind<-1:cv

  numbcat <-sapply(data[,1:cv],FUN=nlevels) #number of categories at each level
# numb<-c(1,cumprod(numbcat[1:(numbvariables-1)])) #number of nodes at each level
  labels<-lapply(data[,1:cv],levels) #the labels for the covariates, for plotting

  #below is copied from expected_counts and into quantile_band - if this changes, so should that
  #make into own function

  ov<-output_variable-1*zip*resp_out #an indicator of how far in the data set you want to look, should be all covariates up to output, and then the output.
  #if output is the Poisson response with variable time, then we include it
  #the risk state is not in the data so we need to subtract it from output_variable in the case that the response is the output

  path_details<-pcegr:::refactored_tree_matrix(data[,1:ov],variable_time = FALSE)#can assume time is false cause we don't need to include
  data_use<-path_details$data_use
  p<-path_details$p
  tree<-path_details$tree_matrix

  #a lot of the below is copied from expected_counts and into quantile_band - if this changes, so should that

  posterior<-mod$posterior.expectation
  stage.struct<-mod$stage.structure

 #  n1<-ov -1*variable_time +1*zip #a number to be used below

  output<-pcegr:::parameter_extractor(stage.struct,posterior,output_variable,poisson_response,remove_risk_free) #used to have n1 instead of output_variable, but now I think they're the same

  if(!(poisson_response&resp_out)){
    m<-dim(output)[2]-1
    output<-as.matrix(output[,1:m])
  }else{
    m<-1
    output<-as.matrix(output)
    }

  num_levels<-lapply(numbcat,pcegr:::vec_from)

  #below copied from Christmas tree plots

  plots<-list()

  if(length(input_variable) > max_per_plot){
    n_plot<-length(input_variable)%/%max_per_plot
    p_plot<-rep(max_per_plot, n_plot)
    rem<-length(input_variable)%%max_per_plot
    if(rem>0){
      p_plot<-c(p_plot,rem)
      n_plot<-n_plot+1
    }
    ind_plot_start<-c(1,cumsum(p_plot)+1)
    ind_plot_end<-ind_plot_start[-1]-1
    ind_plot_start<-ind_plot_start[-(n_plot+1)]
    print(paste0("Number of plots greater than maximum number of plots per page - ",n_plot," pages needed, use back arrow to see previous plots"))
  }else{
    n_plot<-1
    p_plot<-length(input_variable)
    ind_plot_start<-1
    ind_plot_end<-length(input_variable)}

  q<-1

  for(i in input_variable){

    non_input<-cv_ind[-i] #the covariates that are not being examined

    combinations<-expand.grid(num_levels[-i])

    h<-dim(combinations)[1]

    for(j in 1:h){
    v<-combinations[j,]
    if(length(non_input)>1){
    ind<-which(prodlim::row.match(tree[,non_input],v)==1 )
    }else{
      ind<-which(tree[,non_input]==v )
    }
    input_temp<-as.matrix(tree[ind,i])
    output_temp<-as.matrix(output[ind,])

    labels_ind<-as.matrix(numbcat[non_input]-v)

    string<-paste0(names[non_input[1]],"=",labels[[non_input[1]]][labels_ind[1]])
    if(cv>2){
      for(k in 2:(cv-1)){
        string<-paste0(string,",",names[non_input[k]],"=",labels[[non_input[k]]][labels_ind[k]])
      }
    }

    if(j==1){
    data_temp<-data.frame(input=input_temp,output=output_temp,label=string)
    }else{
      data_temp<-rbind(data_temp,data.frame(input=input_temp,output=output_temp,label=string))
    }
    }
    if(poisson_response & resp_out){
      ylabel<-"Rate"
    }else{
      ylabel<-"Prob"
    }

    data_temp$input<-as.factor(data_temp$input)

    max_lim<-ceiling(max(data_temp$output))

    plots[[q]]<-ggplot(data=data_temp,mapping=aes(x=input,y=output,col=label,group=label))+geom_point()+geom_line(position = position_dodge(width = 0.2))+
      xlab(names[i])+ylab(ylabel)+ylim(0,max_lim)+ggtitle(paste0("Effect of ",names[i]," on ",names[output_variable]," given covariates"))
    q<-q+1
  }

  #below copied from Christmas_tree_plots

  for(i in 1:n_plot){
    p.temp<-p_plot[i]
    start.temp<-ind_plot_start[i]
    end.temp<-ind_plot_end[i]
    ind.temp<-c(start.temp:end.temp)
    plots.temp<-list()
    for(j in 1:p.temp){
      plots.temp[[j]]<-plots[[ind.temp[j]]]
    }
    if(p.temp>1){ #added
    rowlim<-ceiling(p.temp/2)
    collim<-2
    }else{
      rowlim<-1
      collim<-1
    }
    print(do.call(ggpubr::ggarrange,list(plotlist=plots.temp,nrow=rowlim,ncol=collim)))
  }
}
