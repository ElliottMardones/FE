#' @import igraph
#' @import parallel
#' @import boot
#' @import poweRlaw
#'

sample_mean <- function(data, indices){
  sample <- data[indices]
  bar <- mean(sample)
  return( (bar))
}

sample_mean_ln<-function(data,indices){
  remuestra<-data[indices]
  remuestra<- remuestra + 1e-6
  m_ln <- conlnorm$new(remuestra)
  m_ln$setXmin(min(remuestra))
  m_ln$pars<-estimate_pars(m_ln)
  return(exp(m_ln$pars[1]))
}

sample_mean_pl<-function(data,indices){
  remuestra<-data[indices]
  remuestra<- remuestra + 1e-6
  m_pl <- conpl$new(remuestra)
  est <- estimate_xmin(m_pl)
  m_pl$setXmin(est)
  m_pl$pars<-estimate_pars(m_pl)
  mediana <- 2^(1/m_pl$pars[1]-1)*m_pl$xmin
  return( mediana )
}

listDataframeToArray3D <- function(data){
  dimData <- dim(data[[1]])
  lengthData <- length(data)
  temp <- array(unlist(data), dim = c(dimData, lengthData))
  rownames(temp) <- rownames(data[[1]])
  colnames(temp) <- colnames(data[[1]])
  message("2D A 3D | cambiar texto y pasar a ingles")
  return(temp)
}

transformMatrix <- function( AA, AB, BB){
  nFilas <- nrow(AA) + nrow(BB)
  nColumnas <- ncol(AA) + ncol(BB)
  nexp <- dim(AA)[3]
  data_output <- array(rep(0,nFilas*nColumnas*nexp), dim = c(nFilas,nColumnas,nexp))
  for( i in seq_len(dim(AA)[3])){
    AAs <- as.data.frame(AA[,,i])
    ABs <- as.data.frame(AB[,,i])
    BBs <- as.data.frame(BB[,,i])
    first <- cbind( AAs, ABs )
    m_output <- as.data.frame(matrix(0, ncol = ncol(AA), nrow = nrow(BB) ))
    colnames(m_output) <- colnames(AA)
    rownames(m_output) <- rownames(BB)
    second <- cbind( m_output, BBs)
    final_matrix <- rbind( first, second)
    data_output[,,i] <- as.matrix(final_matrix)
  }
  allRownames <- c(c(rownames(AA)), c(rownames(BB)))
  allColnames <- c(c(colnames(AA)), c(colnames(BB)))
  rownames(data_output) <- (allRownames)
  colnames(data_output) <- (allColnames)
  return(data_output)
}

validations <- function(CC, CE, EE){
  A <- if( missing(CC))( message("Parameter CC is missing, its required."))else(flag <- "Ok")
  B <- if( missing(CE))( message("Parameter CE is missing, its required."))else(flag <- "Ok")
  C <- if( missing(EE))( message("Parameter EE is missing, its required."))else(flag <- "Ok")
  myFlag <- c(A,B,C)
  if( length(myFlag) == 3){
    if( ncol(CC) == nrow(CE)){
      if( ncol(CE) == nrow(EE )){
        return(flag <- TRUE)
      }else{
        message("The number of columns of CE is different from the number of rows of EE")
        return(flag <- FALSE)
      }
    }else{
      message("The number of columns of CC is different from the number of rows of CE.")
      return(flag <- FALSE)
    }
  }else{
    return(flag <- FALSE)
  }
}

resultIgraph <- function(data){
  nRow <- dim(data)[1]
  nCol <- dim(data)[2]
  experts <- dim(data)[3]
  resultCent <- matrix(nrow = experts, ncol = nRow)
  colnames(resultCent) <- colnames(data[,,1])

  for( i in seq_len(experts)){
    distance <- 1 - data[,,i]
    grafo <- graph.adjacency(distance, weighted = T, mode = "directed")
    cent <- betweenness(grafo, directed = T)
    resultCent[i, ] <- cent
  }
  return(resultCent)
}

#
bs_median_function <- function(data, statistic, R, parallel, ncpus){
  output <- tryCatch(
    bs_median <- boot(data = data, statistic = statistic, R = R, parallel = parallel, ncpus= ncpus), error = function(e) NULL)
  #por el momento no tiene warning)
  return(output)
}
#
conpl_function <- function(resultCent, R, conf,  parallel, ncpus){
  cent_table<-data.frame(var=colnames(resultCent),Median=0,LCI=0,UCI=0, Method=0)#
  for(j in seq_len(nrow(cent_table))){
    bs_median <- bs_median_function(data=resultCent[,j], statistic=sample_mean_pl, R=R, parallel=parallel, ncpus=ncpus)
    if(is.null(bs_median)){
      bs_median         <- bs_median_function(resultCent[,j], statistic= sample_mean, R=R, parallel=parallel, ncpus=ncpus)
      cent_ci           <- boot.ci(bs_median,conf = conf,type = "bca")
      cent_table[j,2:4] <- cbind(bs_median$t0,cent_ci$bca[4],cent_ci$bca[5])
      cent_table[j, 5]  <- "mean"
    }else{
      cent_ci           <- boot.ci(bs_median,conf = conf,type = "bca")
      cent_table[j,2:4] <- cbind(bs_median$t0,cent_ci$bca[4],cent_ci$bca[5])
      cent_table[j, 5]  <- "conpl"
    }
  }
  return(cent_table)
}

conlnorm_function <- function(resultCent, R, conf,  parallel, ncpus){
  #browser()
  cent_table<-data.frame(var=colnames(resultCent),Median=0,LCI=0,UCI=0, Method=0)#
  for(j in seq_len(nrow(cent_table))){
    bs_median <- bs_median_function(data=resultCent[,j], statistic=sample_mean_ln, R=R, parallel=parallel, ncpus=ncpus)
    if(is.null(bs_median)){
      bs_median         <- bs_median_function(resultCent[,j], statistic= sample_mean, R=R, parallel=parallel, ncpus=ncpus)
      cent_ci           <- boot.ci(bs_median,conf = conf,type = "bca")
      cent_table[j,2:4] <- cbind(bs_median$t0,cent_ci$bca[4],cent_ci$bca[5])
      cent_table[j, 5]  <- "mean"
    }else{
      cent_ci           <- boot.ci(bs_median,conf = conf,type = "bca")
      cent_table[j,2:4] <- cbind(bs_median$t0,cent_ci$bca[4],cent_ci$bca[5])
      cent_table[j, 5]  <- "conlnorm"
    }
  }
  return(cent_table)
}

resultBoot_median <- function(resultCent, parallel, reps, ncpus, conf){
  cent_table<-data.frame(var=colnames(resultCent),Median=0,LCI=0,UCI=0, Method = 0)#
  for(j in seq_len(nrow(cent_table))){
    bs_median <- boot(data = resultCent[,j], statistic = sample_mean, R = reps, parallel = parallel, ncpus= ncpus)
    cent_ci<-boot.ci(bs_median,conf = conf,type = "bca")
    cent_table[j,2:4]<-cbind(bs_median$t0,cent_ci$bca[4],cent_ci$bca[5])
    cent_table[j, 5]<- "mean"
  }
  return(cent_table)
}
#
bootCent <- function(CE, Model = "mean" ,parallel ,reps, ncpus, conf){
  if( missing(CE)){
    message("Parameter CC is missing, its required.")
    return(NULL)
  }else{
    #Verificar el funcionamiento
    CE <- if( is.list(CE) == TRUE)  listDataframeToArray3D(CE) else CE
    if( nrow(CE) != ncol(CE)){
      message("Only for square matrices.")
      return(NULL)
    }
    output_resultIgraph <- resultIgraph(CE)
    Model <- ifelse( length(Model) != 1, "mean", Model )
    parallel <- ifelse( length(parallel) != 1, "no", parallel)

    if(Model == "conpl"){
      pl <- conpl_function(resultCent=output_resultIgraph, parallel=parallel, R=reps, ncpus=ncpus, conf=conf)
      return(pl)
    }
    else if(Model == "conlnorm"){
      ln <- conlnorm_function(resultCent=output_resultIgraph, parallel=parallel, R=reps, ncpus=ncpus, conf=conf)
      return(ln)
    }
    else if(Model == "mean"){
      normal_mean <- resultBoot_median(output_resultIgraph, parallel, reps, ncpus, conf)
      return(normal_mean)
    }
  }
}
bootCent_rect <- function(CC,CE,EE, Model = "mean", parallel ,reps, ncpus , conf){
  flag <- validations(CC = CC, CE = CE, EE = EE)
  if( flag == TRUE){
    # Verificar el funcionamiento
    CC <- if( is.list(CC) == TRUE)  listDataframeToArray3D(CC) else CC
    CE <- if( is.list(CE) == TRUE)  listDataframeToArray3D(CE) else CE
    EE <- if( is.list(EE) == TRUE)  listDataframeToArray3D(EE) else EE
    data_matrix <- transformMatrix(CC,CE,EE)
    #
    output_resultIgraph <- resultIgraph(data_matrix)
    Model <- ifelse( length(Model) != 1, "mean", Model )
    print(Model) # Eliminar esta linea despues
    parallel <- ifelse( length(parallel) != 1, "no", parallel)
    print(parallel)

    if(Model == "conpl"){
      pl <- conpl_function(resultCent=output_resultIgraph, parallel=parallel, R=reps, ncpus=ncpus, conf=conf)
      return(pl)
    }
    else if(Model == "conlnorm"){
      ln <- conlnorm_function(resultCent=output_resultIgraph, parallel=parallel, R=reps, ncpus=ncpus, conf=conf)
      return(ln)
    }
    else if(Model == "mean"){
      normal_mean <- resultBoot_median(output_resultIgraph, parallel, reps, ncpus, conf)
      return(normal_mean)
    }
  }else{
    return(NULL)
  }
}
#' @title Significant Centrality For Complete Graphs
#' @name centrality.sq
#' @aliases centrality.sq
#'
#' @description Perform the calculation of the median betweenness centrality with R bootstrap replicas for bipartite graphs.
#'
#' @param CE Three-dimensional matrix, where each submatrix along the z-axis is a square and reflective incidence matrix, or a \code{list} of \code{data.frames} containing square and reflective incidence matrices.
#' @param Model Bootstrap with one of the following statistics: \code{"conpl","conlnorm","mean"}. By default \code{Model = "mean"}.
#' @param parallel The type of parallel operation that is used (if applicable). The options are \code{"multicore", "snow" and "no"}. By default \code{parallel = "no"}.
#' @param reps Number of bootstrap replicas. By default \code{reps = 10,000}.
#' @param ncpus \code{Integer}. Number of processes that are used in the parallel implementation. By default \code{ncpus = 1}.
#' @param conf \code{Real}. Indicates the confidence levels of the required intervals. By default \code{conf = 0.95}.
#'
#' @details
#'The Model parameter makes use of the PoweRlaw package. For \code{“conpl”} the median of a power distribution is calculated according to Newman, M. E. (2005)
#', or \code{"conlnorm"} can be used according to Gillespie CS (2015). In the event that either of the two statistical methods fails,
#' the error will be reported and the mean centrality will be calculated.
#'
#'The parallel and ncpus options are not available on Windows operating systems.
#'
#' @return Returns a data.frame containing the following:
#' \item{Var}{  Variable name.}
#' \item{Median}{ Calculated median.}
#' \item{LCI}{  Lower confidence interval.}
#' \item{UCI}{  Upper confidence interval.}
#' \item{Method}{ Statistical method used.}
#'
#'
#' @references
#'
#' Canty A, Ripley BD (2021). boot: Bootstrap R (S-Plus) Functions. R package version 1.3-28.
#'
#' Csardi G, Nepusz T (2006). "The igraph software package for complex network research." InterJournal, Complex Systems, 1695
#'
#' Gillespie CS (2015). "Fitting Heavy Tailed Distributions: The poweRlaw Package." Journal of Statistical Software, 64(2), 1-16.
#'
#' Newman, M. E. (2005). Power laws, Pareto distributions and Zipf's law. Contemporary physics, 46(5), 323-351.
#'
#' @examples
#' # a <- centrality.rect( CC = m$AA, CE = AB, EE = BB, Model = "conpl")
#' @export

centrality.sq <- function(CE, Model = c("conpl", "conlnorm","mean") ,parallel=c("multicore","snow","no") ,reps = 10000, ncpus = 1, conf =0.95){
  output <- bootCent(CE, Model, parallel, reps, ncpus, conf)
  return(output)
}



#' @title Significant Centrality For Bipartite Graphs
#' @name centrality.rect
#' @aliases centrality.rect
#'
#' @description Perform the calculation of the median betweenness centrality with R bootstrap replicas.
#'
#' @param CC Three-dimensional matrix, where each submatrix along the z-axis is a square and reflective incidence matrix, or a \code{list} of \code{data.frames} containing square and reflective incidence matrices.
#' @param CE Three-dimensional matrix, where each submatrix along the z-axis is a reflective rectangular incidence matrix, or a \code{list} of \code{data.frames} containing reflective and rectangular incidence matrices.
#' @param EE Three-dimensional matrix, where each submatrix along the z-axis is a square and reflective incidence matrix, or a \code{list} of \code{data.frames} containing square and reflective incidence matrices.
#' @param Model Bootstrap with one of the following statistics: \code{"conpl","conlnorm","mean"}. By default \code{Model = "mean"}.
#' @param parallel The type of parallel operation that is used (if applicable). The options are \code{"multicore", "snow" and "no"}. By default \code{parallel = "no"}.
#' @param reps Number of bootstrap replicas. By default \code{reps = 10,000}.
#' @param ncpus \code{Integer}. Number of processes that are used in the parallel implementation. By default \code{ncpus = 1}.
#' @param conf \code{Real}. Indicates the confidence levels of the required intervals. By default \code{conf = 0.95}.
#'
#' @details
#'The Model parameter makes use of the PoweRlaw package. For \code{“conpl”} the median of a power distribution is calculated according to Newman, M. E. (2005)
#', or \code{"conlnorm"} can be used according to Gillespie CS (2015). In the event that either of the two statistical methods fails,
#' the error will be reported and the mean centrality will be calculated.
#'
#'The parallel and ncpus options are not available on Windows operating systems.
#'
#' @return Returns a data.frame containing the following:
#' \item{Var}{  Variable name.}
#' \item{Median}{ Calculated median.}
#' \item{LCI}{  Lower confidence interval.}
#' \item{UCI}{  Upper confidence interval.}
#' \item{Method}{ Statistical method used.}
#'
#'
#' @references
#'
#' Canty A, Ripley BD (2021). boot: Bootstrap R (S-Plus) Functions. R package version 1.3-28.
#'
#' Csardi G, Nepusz T (2006). "The igraph software package for complex network research." InterJournal, Complex Systems, 1695
#'
#' Gillespie CS (2015). "Fitting Heavy Tailed Distributions: The poweRlaw Package." Journal of Statistical Software, 64(2), 1-16.
#'
#' Newman, M. E. (2005). Power laws, Pareto distributions and Zipf's law. Contemporary physics, 46(5), 323-351.
#'
#' @examples
#' # b <- centrality.rect( CC= m$AA, CE = m$AB, EE= m$EE, Model = "conpl")
#' @export
centrality.rect <- function(CC,CE,EE, Model = c("conpl", "conlnorm","mean") ,parallel=c("multicore","snow","no") ,reps = 10000, ncpus = 1, conf = 0.95){
  output <- bootCent_rect(CC,CE,EE, Model, parallel, reps, ncpus, conf)
  return(output)
}

