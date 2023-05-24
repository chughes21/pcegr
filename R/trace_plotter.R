#' The Traceplotter Function
#'
#' Produces a traceplot for chains produced in order to estimate the underlying parameters of the ZIP distribution.
#'
#' @param y_chain A list where each element is a vector corresponding to the chain used to estimate a parameter.
#' @param y_est A vector with the final estimate for each parameter. This will appear as a red horizontal line.
#' @param ind A vector of which chains to plot. If empty, this will produce a trace plot of each chain.
#' @param max_per_plot An integer value specifying the maximum number of chains that can be shown in a single plot.
#'
#' @return A traceplot for each chain provided, along with the final estimate as a red horizontal line.
#' @export
#'
#' @examples
#' gibbs_result<-gibbs_zip(knee_pain_obs)
#' l<-gibbs_result$lambda
#' l_est<-gibbs_result$summary$l_hat
#' trace_plotter(l,l_est)
trace_plotter<-function(y_chain,y_est,ind=c(),max_per_plot = 8){

  #some check for y being a list

  if(length(y_chain)!=length(y_est)){
    stop("Chain and final estimates must correspond")
  }

  p<-length(y_chain)

  leaves<-list()

  if(length(ind)==0){
    ind=c(1:p)
  }

  k<-length(ind)

  #copied from christmas tree functions

  if(k > max_per_plot){
    n_plot<-k%/%max_per_plot
    p_plot<-rep(max_per_plot, n_plot)
    rem<-k%%max_per_plot
    if(rem>0){
      p_plot<-c(p_plot,rem)
      n_plot<-n_plot+1
    }
    ind_plot_start<-c(1,cumsum(p_plot)+1)
    ind_plot_end<-ind_plot_start[-1]-1
    ind_plot_start<-ind_plot_start[-(n_plot+1)]
    print(paste0("Number of leaves greater than maximum number of plots per page - ",n_plot," pages needed, use back arrow to see previous plots"))
  }else{
    n_plot<-1
    p_plot<-k
    ind_plot_start<-1
    ind_plot_end<-k}


  for(i in 1:length(ind)){
    j<-ind[i]
    y<-y_chain[[j]]
    data<-data.frame(chain=y,est=y_est[[j]],index=c(1:length(y)))
    leaves[[i]]<-ggplot(data)+geom_line(mapping=aes(x=index,y=est),col="red")+geom_line(mapping=aes(x=index,y=chain))+
      ggtitle(paste0("Trace Plot - Chain ",j))+xlab("Iteration")+ylab("Estimate")
  }

  #then plot g

  for(i in 1:n_plot){
    p.temp<-p_plot[i]
    start.temp<-ind_plot_start[i]
    end.temp<-ind_plot_end[i]
    ind.temp<-c(start.temp:end.temp)
    leaves.temp<-list()
    for(j in 1:p.temp){
      leaves.temp[[j]]<-leaves[[ind.temp[j]]]
    }
    print(do.call(ggarrange,list(plotlist=leaves.temp,nrow=ceiling(p.temp/2),ncol=2)))
  }

}
