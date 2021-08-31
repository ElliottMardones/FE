#'
#' @import wBoot
#' @import boot

wrapper.de.sq <- function(CC, reps=10000, pValue = 0.05, thr = 0.5, delete = FALSE){
  if( missing(CC)){
    message("Parameter CC is missing, its required.")
    return(NULL)
  }else{
    if(is.list(CC)){
      CC <- listDataframeToArray3D(CC)
    }
    bootCC <- data.frame(From = character(), To = character(), Mean = numeric(), UCI= numeric(), p.value = numeric())
    numdim  <- (dim(CC)[1] * dim(CC)[2] - dim(CC)[2] ) * dim(CC)[3]
    rownamesData <- rownames(CC[,,1])
    colnamesData <- colnames(CC[,,1])
    vector_Value  <- numeric()
    for(x in seq_len(dim(CC)[1]) ){
      for( y in seq_len(dim(CC)[2]) ){
        if( x != y){
          vector_Value <- CC[x,y,]
          valuesFromArrays <- (as.numeric(vector_Value))
          Data.CI<-boot.one.bca(valuesFromArrays,mean,null.hyp = thr, alternative = "less",R=reps)
          bootCC <- rbind( bootCC, data.frame(From = rownamesData[x],
                                              To = colnamesData[y],
                                              Mean = mean(valuesFromArrays),
                                              UCI= Data.CI$Confidence.limits[1],
                                              p.value = Data.CI$p.value))

        }
      }
    }
    if( delete){
      borrar <-  which(bootCC$p.value < pValue | ( bootCC$Mean < thr & is.nan(bootCC$p.value)))
      temp <- bootCC[borrar, ]
      for( ii in seq_len(nrow(temp)) ){
        From <- temp[ii, 1]
        To <- temp[ii, 2]
        From <- which(rownamesData == From)
        To <- which(colnamesData == To)
        CC[From, To, ] <- 0
      }
      #browser()
      if(length(borrar > 0)){
        message("deleting data...")
        bootCC <- bootCC[-borrar, ]
        rownames(bootCC) <- seq_len(nrow(bootCC))
        return(list(Data=CC,directEffects=bootCC ))
      }else{
        message("There is no data to delete...")
        rownames(bootCC) <- seq_len(nrow(bootCC))
        return(list(Data=CC,directEffects=bootCC ))
      }

    }else{
      rownames(bootCC) <- seq_len(nrow(bootCC))
      return(list(directEffects=bootCC ))
    }
  }
}

wrapper.de.rect <- function( CC, CE, EE, reps = 10000, pValue = 0.05, thr = 0.5,  delete = FALSE){
  flag <- validations(CC = CC, CE = CE, EE=EE)
  if( flag == TRUE){
    CCdata <- de.sq(CC, reps, pValue, thr,  delete)
    CEdata <- de.sq(CE, reps, pValue, thr, delete)
    EEdata <- de.sq(EE, reps, pValue, thr, delete)
    directEffects_CC <- CCdata$directEffects
    directEffects_CE <- CEdata$directEffects
    directEffects_EE <- EEdata$directEffects
    #browser()
    output_AllDirecEffects <- rbind(directEffects_CC, directEffects_CE)
    output_AllDirecEffects <- rbind(output_AllDirecEffects,directEffects_EE)
    if(delete == TRUE){
      return(list(CC = CCdata$Data,
                  CE = CEdata$Data,
                  EE = EEdata$Data,
                  directEffects = output_AllDirecEffects))
    }else{
      return(directEffects = output_AllDirecEffects )
    }
  }
}


#' @title Significant Direct Effects For Complete Graphs
#' @name de.sq
#' @aliases de.sq
#'
#' @description Performs the calculation of the mean incidence, left one-sided confidence interval and p-value
#'  with multiple experts. The function allows to eliminate the edges whose mean incidences are not significant
#'  at the set p-value.
#'
#' @param CC Three-dimensional matrix, where each submatrix along the z-axis is a square and reflective incidence matrix,
#'  or a \code{list} of \code{data.frames} containing square and reflective incidence matrices.
#' @param reps Number of bootstrap replicas. By default \code{reps = 10,000}.
#' @param pValue Real between [0,1]. Significance threshold of the left one-sided mean t-test. By default \code{pValue = 0.05}.
#' @param thr \code{Real} between [0,1]. Defines the degree of significant truth. By default \code{thr = 0.5}.
#' @param delete \code{Logical}. For significant = TRUE, remove the rows and columns whose mean occurrences are not significant at the set p-value.
#'  By default \code{significant = FALSE.}

#'
#' @details The function implements "boot.one.bca" from the wBoot package to get the UCI and the p-value.
#'
#' @return
#' For \code{delete = TRUE}, it returns the three-dimensional matrix entered, but with the non-significant rows and columns
#' removed and a "directEffects" \code{data.frame} containing the following values:
#' \item{From}{Origin of the incident}
#' \item{To}{Destination of the incident}
#' \item{Mean}{Average incidence}
#' \item{UCI}{Upper Confidence Interval}
#' \item{p-value}{The p-value for the t test}
#' For \code{delete = FALSE}, only the \code{data.frame} "directEffects" is returned.
#'
#'
#' @examples
#' # saveData <- de.sq(CC =m$AA, reps = 1000, pValue = 0.05, thr = 0.5, delete = TRUE)
#' @export
de.sq <- function(CC, reps=10000, pValue = 0.05, thr = 0.5, delete = FALSE){
  output <- wrapper.de.sq(CC, reps, pValue, thr, delete)
  return(output)
}


#' @title Significant Direct Effects For Bipartite Graphs
#' @name de.rect
#' @aliases de.rect
#' @description Performs the calculation of the mean incidence, left one-sided confidence interval and p-value
#' with multiple experts for rectangular matrices. The function allows to eliminate the edges whose mean
#' incidences are not significant at the set p-value.
#'
#' @param CC Three-dimensional matrix, where each submatrix along the z-axis is a square and reflective incidence matrix,
#'  or a \code{list} of \code{data.frames} containing square and reflective incidence matrices.
#' @param CE Three-dimensional matrix, where each submatrix along the z-axis is a reflective rectangular incidence matrix,
#' or a \code{list} of \code{data.frames} containing reflective and rectangular incidence matrices.
#' @param EE Three-dimensional matrix, where each submatrix along the z-axis is a square and reflective incidence matrix,
#'  or a \code{list} of \code{data.frames} containing square and reflective incidence matrices.
#' @param reps Number of bootstrap replicas. By default \code{reps = 10,000}.
#' @param pValue Real between [0,1]. Significance threshold of the left one-sided mean t-test. By default \code{pValue = 0.05}.
#' @param thr \code{Real} between [0,1]. Defines the degree of significant truth. By default \code{thr = 0.5}.
#' @param delete \code{Logical}. For significant = TRUE, remove the rows and columns whose mean occurrences are not significant at the set p-value.
#'  By default \code{significant = FALSE.}

#'
#' @details The function implements "boot.one.bca" from the wBoot package to get the UCI and the p-value.
#'
#' @return
#' For \code{delete = TRUE}, it returns the three-dimensional matrix entered, but with the non-significant rows and columns
#' removed and a "directEffects" \code{data.frame} containing the following values:
#' \item{From}{Origin of the incident}
#' \item{To}{Destination of the incident}
#' \item{Mean}{Average incidence}
#' \item{UCI}{Upper Confidence Interval}
#' \item{p-value}{The p-value for the t test}
#' For delete = FALSE, only the \code{data.frame} "directEffects" is returned.
#'
#' @export
de.rect <- function( CC, CE, EE, reps = 10000, pValue = 0.05, thr = 0.5,  delete = FALSE){
  output <- wrapper.de.rect( CC, CE, EE, reps, pValue, thr, delete)
  return(output)
}
