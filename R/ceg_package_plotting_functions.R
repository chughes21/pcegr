
TreeGraph1<-function(stagedtree, name = c(), range.color = 1)
{
  tree<-stagedtree$event.tree
  solution<-stagedtree$stage.structure
  post<-stagedtree$posterior.expectation

  n<-tree$num.variable
  m<-tree$num.situation[n]

  if(tree$poisson.response){
  rates<-parameter_extractor(solution,post,n,tree$poisson.response,stagedtree$remove.risk.free.edges)
  }else{rates<-NA}

  if(tree$poisson.response & stagedtree$remove.risk.free.edges){
    probs<-parameter_extractor(solution,post,n-1,tree$poisson.response,stagedtree$remove.risk.free.edges)
    probs<-probs[,2]
  }else{probs<-NA}

  nodes <- NodeSet1(tree)
  edgeList <- EdgeList1(tree, nodes)
  node.label <- NodeLabel1(tree$num.variable, tree$num.situation,
                          tree$num.category, name)
  edge.label <- EdgeLabel1(tree$num.variable, tree$num.situation,
                          tree$label.category, poisson_response = tree$poisson.response,
                          remove_risk_free = stagedtree$remove.risk.free.edges, rates=rates,probs = probs)
  node.color <- NodeColor1(tree$num.variable, tree$num.situation,
                          tree$num.category, solution, range.color)
  graph <- list()
  graph$node <- list()
  graph$node$nodes <- nodes
  graph$node$label <- node.label
  graph$node$color <- node.color
  graph$edge <- list()
  graph$edge$edges <- edgeList
  graph$edge$label <- edge.label
  return(graph)
}

EventTreeGraph1<-function (event.tree)
{
  nodes <- NodeSet1(event.tree)
  edgeList <- EdgeList1(event.tree, nodes)
  node.label <- NodeLabel1(event.tree$num.variable, event.tree$num.situation,
                          event.tree$num.category, names(data))
  edge.label <- EdgeLabel1(event.tree$num.variable, event.tree$num.situation,
                          event.tree$label.category, FALSE, FALSE) #even when these are true, we don't need them for the event tree
  total.node <- cumsum(event.tree$num.situation)
  node.color <- rep("white", total.node[event.tree$num.variable] +
                      event.tree$num.category[event.tree$num.variable] * event.tree$num.situation[event.tree$num.variable])
  graph <- list()
  graph$node <- list()
  graph$node$nodes <- nodes
  graph$node$label <- node.label
  graph$node$color <- node.color
  graph$edge <- list()
  graph$edge$edges <- edgeList
  graph$edge$label <- edge.label
  graph$node$color <- rep("white", length(graph$node$color))
  return(graph)
}

NodeSet1<-function (tree)
{
  tree
  num.node <- sum(tree$num.situation) + tree$num.situation[tree$num.variable] *
    tree$num.category[tree$num.variable]
  node <- paste("s", 1:num.node - 1, sep = "")
  return(node)
}

EdgeList1<-function (event.tree, node)
{
  start.situation <- cumsum(event.tree$num.situation)
  max <- sum(event.tree$num.situation) + event.tree$num.situation[event.tree$num.variable] *
    event.tree$num.category[event.tree$num.variable]
  start.situation <- start.situation + 1
  start.situation <- c(1, start.situation)
  edge.list <- lapply(1:max, function(x) EdgeSituation1(x, start.situation,
                                                       event.tree$num.category))
  names(edge.list) <- node
  return(edge.list)
}

NodeLabel1<-function (num.variable, num.situation, num.category, label)
{
  result <- sapply(1:num.variable, function(x) rep(label[x],
                                                   num.situation[x]))
  result <- ListToVector1(result, num.variable)
  num.leaf <- num.category[num.variable] * num.situation[num.variable]
  aux <- sapply(1:num.leaf, function(x) paste("leaf",
                                              x))
  result <- c(result, aux)
  return(result)
}

EdgeLabel1<-function (num.variable, num.situation, label, poisson_response, remove_risk_free, rates = NA, probs = NA)
{
  result <- lapply(1:num.variable, function(x) rep(label[[x]],
                                                   num.situation[x]))

  rates<-round(rates,2)
  probs<-round(probs,2)
  if(poisson_response){
    result[[num.variable]]<-paste0(result[[num.variable]]," (",rates,")")
  }

  if(remove_risk_free & poisson_response){
    result[[num.variable-1]]<-paste0(result[[num.variable-1]]," (",probs,")")
  }

  result <- ListToVector1(result, num.variable)
  return(result)
}

NodeColor1<-function (num.variable, num.situation, num.category, stage.structure,
                      range.color)
{
  total.node <- cumsum(num.situation)
  result <- rep("white", total.node[num.variable] + num.category[num.variable] *
                  num.situation[num.variable])
  count <- 2
  if (range.color == 1) {
    color <- palette()
    color[1] <- "white"
    color<-c(color,"#B06500","#0B6623","#E39ff6")
    start_count<-2
  }
  else {
    if (range.color == 2) {
      color <- colors(1)
      color <- color[-21]
      start_count<-1
    }
  }

  max_count<-length(color)

  for (i in 2:num.variable) {
    for (j in 1:num.situation[i]) {
      if (!is.na(stage.structure[[i]][[j]][1])) {
        if (length(stage.structure[[i]][[j]]) == 1)
          result[j + total.node[i - 1]] <- "white"
        else {
          result[stage.structure[[i]][[j]] +
                   total.node[i - 1]] <- color[count]
          count <- count + 1
          if(count > max_count){
            count<-start_count
          }
        }
      }
    }
  }
  return(result)
}

EdgeSituation1<-function (situation, start.situation, num.category)
{
  max <- start.situation[length(start.situation)]
  if (situation >= max)
    return(list())
  aux <- findInterval(situation, start.situation)
  result <- start.situation[aux + 1] - 2 + (situation - start.situation[aux]) *
    num.category[aux] + 1:num.category[aux]
  result <- paste("s", result, sep = "")
  result <- list(edges = result)
  return(result)
}

ListToVector1<-function (x, n)
{
  if (n < 1)
    return(c())
  return(c(ListToVector1(x, n - 1), x[[n]]))
}

MergeLabels1<-function (edge.list, edge, level)
{
  aux.merge <- which(edge.list == edge)
  aux <- length(aux.merge)
  aux.label <- level[aux.merge[1]]
  if (aux > 1) {
    for (i in 2:aux) {
      aux.label <- paste0(aux.label, "-", level[aux.merge[i]])
    }
  }
  return(aux.label)
}

CegGraphSimple1<-function (staged.tree, position, range.color = 1)
{
  event.tree<-staged.tree$event.tree

  node.vector <- c()
  node.variable <- c()
  node.color <- c()
  edge.list <- list()
  edge.label <- c()
  edge.weight <- c()
  count.color <- 2
  count.pos <- -1
  if (range.color == 1) {
    color <- grDevices::palette()
    color[1] <- "white"
  }
  else {
    color <- colors(1)
    color <- color[-21]
  }
  for (var in 1:(event.tree$num.variable - 1)) {
    start.pos <- count.pos
    edge.var.list <- c()
    num.situation <- length(position[[var]])
    for (stage in 1:num.situation) {
      num.pos <- length(position[[var]][[stage]])
      pos.next.var <- PositionVector1(event.tree$num.situation[var + 1], position[[var + 1]])
      for (pos in 1:num.pos) {
        count.pos <- count.pos + 1
        node.vector <- c(node.vector, paste0("w",
                                             count.pos))
        node.variable <- c(node.variable, var)
        if (num.pos == 1 & (length(position[[var]][[stage]][[1]])==1))
          node.color <- c(node.color, color[1])
        else node.color <- c(node.color, color[count.color])
        aux <- (position[[var]][[stage]][[pos]][1] -
                  1) * event.tree$num.category[var]
        aux <- (aux + 1):(aux + event.tree$num.category[var])
        edge.var.list <- c(edge.var.list, pos.next.var[aux])
        edge.weight <- c(edge.weight, rep(round(1/event.tree$num.category[var],
                                                2), event.tree$num.category[var]))
      }
      if (num.pos != 1 | (length(position[[var]][[stage]][[1]])>1))
        count.color <- count.color + 1
      if(count.color>length(color)){ #added
        count.color<-2 #added
      }#added
    }
    edge.var.list <- edge.var.list + count.pos
    dim(edge.var.list) <- c(event.tree$num.category[var],
                            length(edge.var.list)/event.tree$num.category[var])
    edge.var.list <- as.matrix(edge.var.list)
    for (pos in start.pos:(count.pos - 1)) {
      edge.list[[pos + 2]] <- list()
      aux.edge.var.list <- unique(edge.var.list[, pos -
                                                  start.pos + 1])
      edge.label <- c(edge.label, sapply(1:length(aux.edge.var.list),
                                         function(x) MergeLabels1(edge.var.list[, pos -
                                                                                 start.pos + 1], aux.edge.var.list[x], event.tree$label.category[[var]])))
      edge.list[[pos + 2]]$edges <- paste0("w", aux.edge.var.list)
    }
  }
  var <- event.tree$num.variable
  num.pos <- length(position[[var]])
  node.vector <- c(node.vector, paste0("w", count.pos +
                                         1:num.pos))
  node.variable <- c(node.variable, rep(var, num.pos))
  aux.color<-c()
  for(i in 1:num.pos){
    if(length(position[[var]][[i]][[1]])>1){
      aux.color<-c(aux.color,color[count.color])
      count.color<-count.color+1
      if(count.color>length(color)){
        count.color<-2
      }
    }else{aux.color<-c(aux.color,color[1])}
  }
  node.color <- c(node.color,aux.color)
  aux.label <- event.tree$label.category[[var]][1]
  if(event.tree$poisson.response){
    rates<-round(staged.tree$posterior.expectation[[var]],2)
    rates<-rates[-which(is.na(rates))]
    if(staged.tree$remove.risk.free.edges){
      rates<-rates[-1]
    }
    aux.label<-rep(aux.label,num.pos)
    aux.label<-paste0(aux.label," (",rates,")")
    edge.label<-c(edge.label,aux.label)
  }else{
  for (i in 2:event.tree$num.category[var]) {
    aux.label <- paste0(aux.label, "-", event.tree$label.category[[var]][i])
    }
  edge.label <- c(edge.label, rep(aux.label, num.pos))
  }
  edge.weight <- c(edge.weight, rep(rep(round(1/event.tree$num.category[var],
                                              2), event.tree$num.category[var]), num.pos))
  ref <- count.pos + num.pos + 1
  for (pos in 1:num.pos) {
    edge.list[[pos + count.pos + 1]] <- list()
    edge.list[[pos + count.pos + 1]]$edges <- paste0("w",
                                                     ref)
  }
  ref <- ref + 1
  node.vector <- c(node.vector, paste0("w", ref - 1))
  node.color <- c(node.color, color[1])
  edge.list[[ref]] <- list()
  names(edge.list) <- node.vector
  graph <- list()
  graph$node <- list()
  graph$node$nodes <- node.vector
  graph$node$variable <- node.variable
  graph$node$color <- node.color
  graph$edge <- list()
  graph$edge$edges <- edge.list
  graph$edge$label <- edge.label
  graph$edge$weight <- edge.weight
  return(graph)
}
