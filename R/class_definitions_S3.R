EventTree<-function(data,poisson_response = TRUE, variable_time=TRUE){
            if(!poisson_response & variable_time){
              stop("Poisson response needed for variable time")
            }
            n<-dim(data)[2]-1*poisson_response - 1*variable_time

            if(poisson_response){
              empty.resp<-factor(rep("y",length(data[,1])),levels=c("y"))
              data.frame<-data.frame(data[,1:n],resp=empty.resp)
            }else{
              data.frame<-data
            }
            num.variable <- length(data.frame[1, ])
            num.slice <- ncol(data.frame)/num.variable
            label.category <- lapply(1:num.variable, function(x) levels(data.frame[, x])) #changed for better labels
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
                       num.slice = num.slice, mergedlist = mergedlist, poisson.response = poisson_response),class="EventTree"))
          }


StagedTree<-function(data,prior,posterior,stage.struc,merged,result,poisson_response = TRUE, variable_time = TRUE,zip=FALSE,lik = 0){
  if(!poisson_response & variable_time){
    stop("Variable time requires a Poisson response")
  }

  if(!poisson_response & zip){
    stop("Zero Inflated Poisson requires Poisson response")
  }

  if(zip){
    n<-dim(data)[2]-1*poisson_response -1*variable_time
    state<-factor(rep("No Risk",length(data[,1])),levels<-c("No Risk", "Risk"))
    data.final<-data.frame(data[,1:n], State = state, data[,-(1:n)])
  }else{
    data.final<-data
  }

  event.tree<-EventTree(data.final,poisson_response,variable_time)

  return(structure(list(event.tree = event.tree,
                   situation = list(), contingency.table = list(), merged = merged, stage.structure = stage.struc,
                   result = result,stage.probability = list(), prior.distribution=prior, posterior.distribution=posterior,
                   model.score=lik),class="StagedTree"))
}

ChainEventGraph<-function(staged.tree){
  position <- StratifiedCegPosition1(staged.tree$stage.structure,
                                     staged.tree$event.tree$num.category, staged.tree$event.tree$num.situation)
  return(structure(list(staged.tree = staged.tree,
             position = position),class="ChainEventGraph"))
}

plot.EventTree<-function(x){
            event.tree.graph <- EventTreeGraph1(x)
            g <- new("graphNEL", nodes = event.tree.graph$node$nodes,
                     edgeL = event.tree.graph$edge$edges, edgemode = "directed")
            attrsAtt <- list()
            graphAtt <- list(rankdir = "LR", size = "18.0,24.0",
                             bgcolor = "white")
            edgeAtt <- list(color = "cyan")
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
           grDevices::pdf("./event-tree-plot.pdf",
                           width = 8, height = 6, title = "")
            Rgraphviz::plot(g, main = "Event Tree Graph",
                            nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
            grDevices::dev.off()
}

plot.StagedTree<-function(x){
  staged.tree.graph <- TreeGraph1(x$event.tree, x$stage.structure)
  g <- new("graphNEL", nodes = staged.tree.graph$node$nodes,
           edgeL = staged.tree.graph$edge$edges, edgemode = "directed")
  attrsAtt <- list()
  graphAtt <- list(rankdir = "LR", size = "18.0,24.0",
                   bgcolor = "white")
  edgeAtt <- list(color = "cyan")
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
  grDevices::pdf("./staged.tree.graph.pdf", width = 8,
                 height = 6, title = "")
  Rgraphviz::plot(g, main = "Staged Tree Graph",
                  nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
  grDevices::dev.off()
}

plot.ChainEventGraph<-function(x){
    ceg.graph.simple <- CegGraphSimple1(x$staged.tree$event.tree,
                                       x$position)
    g <- new("graphNEL", nodes = ceg.graph.simple$node$nodes,
             edgeL = ceg.graph.simple$edge$edges, edgemode = "directed")
    attrsAtt <- list()
    graphAtt <- list(rankdir = "LR", size = "18.0,24.0",
                     bgcolor = "white")
    edgeAtt <- list(color = "cyan")
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
    grDevices::pdf("./figure02ceg.pdf", width = 8,
                  height = 6, title = "Chain Event Graph")
    Rgraphviz::plot(g, main = "Chain Event Graph",
                    nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrsAtt)
    grDevices::dev.off()

}

