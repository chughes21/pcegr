#' The EventTree constructor function
#'
#' Create an object in the StagedTree S3 class from a data set.
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param zip A logical value indicating whether the process is zero-inflated (TRUE) or not (FALSE). If the process is zero-inflated but this has already been accounted for in the data input, this should still be FALSE.
#' @param remove_risk_free A logical value indicating whether the risk free leaves and edges should be removed (TRUE) or not (FALSE).
#'
#' @return An object of the S3 class EventTree
#' @export
#'
#' @examples
#' tree<-EventTree(knee_pain_obs,TRUE,TRUE)
#' plot(tree)
EventTree<-function(data,poisson_response = TRUE, variable_time=TRUE,zip = FALSE, remove_risk_free = FALSE){
            if(!poisson_response & variable_time){
              stop("Poisson response needed for variable time")
            }

             if(!poisson_response & zip){
               stop("Zero Inflated Poisson requires Poisson response")
            }

            n<-dim(data)[2]-1*poisson_response - 1*variable_time

            if(poisson_response){
              empty.resp<-factor(rep("y",length(data[,1])),levels=c("y"))
              data.temp<-data.frame(data[,1:n],resp=empty.resp)
            }else{
              data.temp<-data
            }

              if(remove_risk_free){
              state<-factor(rep("Risk",length(data[,1])),levels<-c("Risk"))
              n1<-n-1+1*zip
              data.final<-data.frame(data.temp[,1:n1], State = state, data.temp[,-(1:n)]) #n in the second one is supposed to be like that
              }else if(zip){
              state<-factor(rep("Risk",length(data[,1])),levels<-c("No Risk", "Risk"))
              data.final<-data.frame(data.temp[,1:n], State = state, data.temp[,-(1:n)])
              }else{
              data.final<-data.temp
            }

            num.variable <- length(data.final[1, ])
            num.slice <- ncol(data.final)/num.variable
            label.category <- lapply(1:num.variable, function(x) levels(data.final[, x])) #changed for better labels
            num.category <- c()
            num.category <- sapply(label.category, length)
            num.situation <- c(1, cumprod(num.category[1:(num.variable -  1)]))
            begin.stage <- c(0, cumsum(num.situation[1:(num.variable -  1)]))
            mergedlist <- sapply(1:(num.variable - 1), function(x) LabelStage1(x, num.variable, num.situation, label.category, num.category))
            mergedlist <- lapply(1:(num.variable), function(x) {
              lapply(seq_len(num.situation[x]), function(y) mergedlist[y + begin.stage[x], ])
            })
            return(structure( list(num.variable = num.variable,
                       num.category = num.category, label.category = label.category, num.situation = num.situation,
                       num.slice = num.slice, situation.list = mergedlist, poisson.response = poisson_response, variable.time=variable_time),class="EventTree"))
          }


#' The StagedTree constructor function
#'
#' Create an object in the StagedTree S3 class from certain inputs for a model. Some fo these inputs can be calculated directly in the [pceg()] function.
#'
#' @param data A data set where the observed response vector and time vector (if applicable and variable) are the last two columns.
#' @param prior A list detailing the prior distribution for the model.
#' @param counts A list specifying the edge counts or event counts and times for the model.
#' @param posterior A list specifying the posterior expectations for the model.
#' @param stage.struc A list specifying the stage structure for the model.
#' @param stages A numeric vector indicating which of the original situations are now stages.
#' @param merged A matrix specifying which situations have been merged into stages.
#' @param result A list detailing the covariate values in each stage.
#' @param poisson_response A logical value indicating whether the response variable is Poisson (TRUE) or categorical (FALSE).
#' @param variable_time A logical value indicating whether the observed time is uniform (FALSE) or variable (TRUE), if applicable.
#' @param zip A logical value indicating whether the process is zero-inflated (TRUE) or not (FALSE). If the process is zero-inflated but this has already been accounted for in the data input, this should still be FALSE.
#' @param remove_risk_free A logical value indicating whether the risk free leaves and edges should be removed (TRUE) or not (FALSE).
#' @param lik A numeric value specifying the log marginal likelihood for the model.
#'
#' @return An object of the S3 class StagedTree
#'
StagedTree<-function(data,prior,counts,posterior,stage.struc,stages,merged,result,poisson_response = TRUE, variable_time = TRUE,zip=FALSE,remove_risk_free = FALSE,lik = 0){
  if(!poisson_response & variable_time){
    stop("Variable time requires a Poisson response")
  }

  if(!poisson_response & zip){
    stop("Zero Inflated Poisson requires Poisson response")
  }

  event.tree<-EventTree(data,poisson_response,variable_time,zip,remove_risk_free)

  return(structure(list(event.tree = event.tree,stages = stages, remove.risk.free.edges = remove_risk_free, merged = merged, stage.structure = stage.struc,
                   result = result, prior.distribution=prior, data.summary = counts, posterior.expectation=posterior,
                   model.score=lik),class="StagedTree"))
}

#' The ChainEventGraph constructor function.
#'
#' Create an object in the S3 class ChainEventGraph using an object in the S3 class StagedTree.
#'
#' @param staged.tree An object in the S3 class StagedTree
#'
#' @return An object in the S3 class ChainEventGraphs
#' @export
#'
#' @examples
#' mod<-pceg(knee_pain_obs,2,TRUE,TRUE)
#' ceg<-ChainEventGraph(mod)
#' plot(ceg)
ChainEventGraph<-function(staged.tree){

  if(staged.tree$event.tree$num.variable<=2){
    stop("ChainEventGraph constructor only defined for StagedTrees with more than two variables.")
  }

  position <- StratifiedCegPosition1(staged.tree$stage.structure,
                                     staged.tree$event.tree$num.category, staged.tree$event.tree$num.situation)
  return(structure(list(staged.tree = staged.tree,
             position = position),class="ChainEventGraph"))
}

#' @export
#' @method plot EventTree
plot.EventTree<-function(x, ... ){
            event.tree.graph <- EventTreeGraph1(x)
            g <- new("graphNEL", nodes = event.tree.graph$node$nodes,
                     edgeL = event.tree.graph$edge$edges, edgemode = "directed")
            attrsAtt <- list()
            graphAtt <- list(rankdir = "LR", size = "18.0,24.0",
                             bgcolor = "white")
            edgeAtt <- list(color = "deepskyblue")
            nodeAtt <- list(fillcolor = "lightgreen", shape = "ellipse",
                            fixedsize = FALSE)
            attrsAtt <- list(node = nodeAtt, edge = edgeAtt, graph = graphAtt)
            nodes.label.list <- event.tree.graph$node$nodes
            names(nodes.label.list) <- graph::nodes(g)
            nAttrs <- list()
            nAttrs$label <- nodes.label.list
            edges.label.list <- event.tree.graph$edge$label
            names(edges.label.list) <- graph::edgeNames(g)
            eAttrs <- list()
            eAttrs$label <- edges.label.list
            nAttrs$fillcolor <- event.tree.graph$node$color
            names(nAttrs$fillcolor) <- graph::nodes(g)
            grDevices::graphics.off()
            graphics::par("mar")
            graphics::par(mar = c(1, 1, 1, 1))
            Rgraphviz::plot(g, main = "Event Tree Graph",
                            nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
         #  grDevices::pdf("./event-tree-plot.pdf",
         #                  width = 8, height = 6, title = "")
         #   Rgraphviz::plot(g, main = "Event Tree Graph",
         #                   nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
         #   grDevices::dev.off()
}


#' @export
#' @method plot StagedTree
plot.StagedTree<-function(x, ...){

  staged.tree.graph <- TreeGraph1(x)
  g <- new("graphNEL", nodes = staged.tree.graph$node$nodes,
           edgeL = staged.tree.graph$edge$edges, edgemode = "directed")
  attrsAtt <- list()
  graphAtt <- list(rankdir = "LR", size = "18.0,24.0",
                   bgcolor = "white")
  edgeAtt <- list(color = "deepskyblue")
  nodeAtt <- list(fillcolor = "lightgreen", shape = "ellipse",
                  fixedsize = FALSE)
  attrsAtt <- list(node = nodeAtt, edge = edgeAtt, graph = graphAtt)
  nodesLabelList <- staged.tree.graph$node$nodes
  names(nodesLabelList) <- graph::nodes(g)
  nAttrs <- list()
  nAttrs$label <- nodesLabelList
  edgesLabelList <- staged.tree.graph$edge$label
  names(edgesLabelList) <- graph::edgeNames(g)
  eAttrs <- list()
  eAttrs$label <- edgesLabelList
  nAttrs$fillcolor <- staged.tree.graph$node$color
  names(nAttrs$fillcolor) <- graph::nodes(g)
  Rgraphviz::plot(g, main = "Staged Tree Graph",
                  nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
 # grDevices::pdf("./staged.tree.graph.pdf", width = 8,
 #                 height = 6, title = "")
 # Rgraphviz::plot(g, main = "Staged Tree Graph",
 #                 nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
 # grDevices::dev.off()

}

#' @export
#' @method plot ChainEventGraph
plot.ChainEventGraph<-function(x, ...){
    ceg.graph.simple <- CegGraphSimple1(x$staged.tree,
                                       x$position)
    g <- new("graphNEL", nodes = ceg.graph.simple$node$nodes,
             edgeL = ceg.graph.simple$edge$edges, edgemode = "directed")
    attrsAtt <- list()
    graphAtt <- list(rankdir = "LR", size = "18.0,24.0",
                     bgcolor = "white")
    edgeAtt <- list(color = "deepskyblue")
    nodeAtt <- list(fillcolor = "lightgreen", shape = "ellipse",
                    fixedsize = FALSE)
    attrsAtt <- list(node = nodeAtt, edge = edgeAtt, graph = graphAtt)
    nodes.label.list <- ceg.graph.simple$node$nodes
    names(nodes.label.list) <- graph::nodes(g)
    nAttrs <- list()
    nAttrs$label <- nodes.label.list
    edges.label.list <- ceg.graph.simple$edge$label
    names(edges.label.list) <- graph::edgeNames(g)
    eAttrs <- list()
    eAttrs$label <- edges.label.list
   # Rgraphviz::plot(g, main = "Chain Event Graph (propagation analysis) ",
    #                nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
   # grDevices::pdf("./figure01ceg.pdf", width = 8,
    #               height = 6, title = "Chain Event Graph (propagation analysis)")
    #Rgraphviz::plot(g, main = "Chain Event Graph (propagation analysis)",
     #               nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
    #grDevices::dev.off()
    nAttrs$fillcolor <- ceg.graph.simple$node$color
    names(nAttrs$fillcolor) <- graph::nodes(g)
    Rgraphviz::plot(g, main = "Chain Event Graph",
                    nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
  #  grDevices::pdf("./figure02ceg.pdf", width = 8,
  #                  height = 6, title = "Chain Event Graph")
  #  Rgraphviz::plot(g, main = "Chain Event Graph",
  #                    nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
  #  grDevices::dev.off()

}

#' @export
#' @method summary StagedTree
summary.StagedTree<-function(object){

    post<-object$posterior.expectation

    m<-length(post)

    for(i in 1:m){
      temp<-as.matrix(post[[i]])
      ind<-which(is.na(temp[,1]))
      if(length(ind)>0){
      post[[i]]<-temp[-ind,]
      }
    }

    print("Merged Stages")
    print(object$result)
    print("Posterior Expectations")
    print(post)
    print("Log Marginal Likelihood")
    print(object$model.score)

}


#' @export
#' @method summary ChainEventGraph
summary.ChainEventGraph<-function(object){
  print("Position Structure")
  print(object$position)
  summary(object$staged.tree)
}
