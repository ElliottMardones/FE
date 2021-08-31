#'
#' @import wBoot
#' @import boot
wrapper.BootMargin<-function(CC,  reps=10000,  pValue = 0.05, thr = 0.5, delete=FALSE){
  if( missing(CC)){
    message("Parameter CC is missing, its required.")
    return(NULL)
  }
  if( nrow(CC) == ncol(CC)){
    nn<-nrow(CC)
    promFilas<-data.frame(Var=colnames(CC[,,1]),Mean=0,UCI=0,p.value=0)
    promColumnas<-data.frame(Var=colnames(CC[,,1]),Mean=0,UCI=0,p.value=0)
    #
    saveRowData  <- data.frame(Var=character() ,Mean=numeric(),UCI=numeric(),p.value=numeric())
    saveColData  <- data.frame(Var=character(),Mean=numeric(),UCI=numeric(),p.value=numeric())
    #
    for(i in 1:nn){
      filasExp<-t(CC[i,,])[,-i] # FILA
      fila<-rowMeans(filasExp)
      fila.CI<-boot.one.bca(fila, mean, null.hyp = thr , alternative = "less", R=reps)
      fila.Mean<-fila.CI$Mean
      fila.UCI<-fila.CI$Confidence.limits[1]
      fila.p_value<-fila.CI$p.value
      #####################################
      columnaExp<-t(CC[,i,])[,-i] # COLUMNA
      columna<-rowMeans(columnaExp)
      columna.CI<-boot.one.bca(columna, mean, null.hyp = thr, alternative =  "less", R=reps)
      columna.Mean<-columna.CI$Mean
      columna.UCI<-columna.CI$Confidence.limits[1]
      columna.p_value<-columna.CI$p.value
      promFilas[i,2:4]    <-cbind(fila.Mean, fila.UCI, fila.p_value)
      promColumnas[i,2:4] <-cbind(columna.Mean ,columna.UCI, columna.p_value)
    }
    if( delete ){
      #
      saveRowData <-(promFilas) # Eliminar esto para eliminar las 3ra salida.
      saveColData <- (promColumnas) # Eliminar esto para eliminar la 4ta salida.
      #
      delete_in_rows <-  which(promFilas$p.value < pValue | ( promFilas$Mean < thr & is.nan(promFilas$p.value)))
      delete_in_cols <- which(promColumnas$p.value < pValue | (promColumnas$Mean < thr & is.nan(promColumnas$p.value)))
      promFilas <- promFilas[-delete_in_rows, ]
      promColumnas <- promColumnas[-delete_in_cols, ]
      #
      rn_mR<- rownames(promFilas)
      rn_mC <- rownames(promColumnas)
      marginRow <- which( rn_mR %in% rn_mC)
      marginCol <- which( rn_mC %in% rn_mR)

      promFilas_output <- promFilas[ marginRow, ]
      promColumnas_output <- promColumnas[ marginCol, ]
      rownames(promFilas_output) <- seq_len(nrow(promFilas_output))
      rownames(promColumnas_output) <- seq_len(nrow(promColumnas_output))
      #browser()
      if( length(marginRow) ==0 ){
        message("All data has been deleted...")
        return(NULL)
      }
      new_CC <- CC[promFilas_output$Var, promColumnas_output$Var, ]
      return(list(Data=new_CC,
                  MarginRow =  promFilas_output,
                  MarginCol = promColumnas_output))

    }
    return(list(MarginRow = promFilas,
                MarginCol = promColumnas))

  }else{
    message("Only for complete graphs")
  }
}


#' @title Significant Causes and Effects For Complete Graphs
#' @name bootMargin
#' @aliases bootMargin
#'
#' @description Performs the calculation of the mean incidence for each cause and each effect for all experts,
#'  the left one-sided confidence interval and the p-value. The function makes it possible to eliminate causes
#'  and effects whose average incidence is not significant at the established p-value.
#'
#' @param CC Three-dimensional matrix, where each submatrix along the z-axis is a square and reflective incidence matrix,
#'  or a \code{list} of \code{data.frames} containing square and reflective incidence matrices.
#' @param reps Number of bootstrap replicas. By default \code{reps = 10,000}.
#' @param pValue Real between [0,1]. Significance threshold of the left one-sided mean t-test. By default \code{pValue = 0.05}.
#' @param thr \code{Real} between [0,1]. Defines the degree of significant truth. By default \code{thr = 0.5}.
#' @param delete \code{Logical}. For significant = TRUE, remove the rows and columns whose mean occurrences are not significant at the set p-value.
#'  By default \code{significant = FALSE}.

#'
#' @details The function implements "boot.one.bca" from the wBoot package to get the UCI and the p-value.
#'
#' @return
#' For \code{delete = TRUE}, it returns the three-dimensional matrix entered, but with the non-significant rows and columns
#' removed and a "directEffects" \code{data.frame} containing the following values:
#' \item{Var}{Variable name}
#' \item{Mean}{Calculated mean}
#' \item{ICU}{Upper Confidence Interval}
#' \item{p-value}{The p-value for the t test}
#'
#'
#' @examples
#' # Example of use
#' # saveData <- bootMargin(CC = m$AA, reps = 1000, pValue = 0.05, thr = 0.5, delete = TRUE)
#' # To see the new data set:
#' ##    saveData$Data
#' # To see the results obtained by row:
#' ##    savaData$MarginRow
#' # To see the results obtained by col:
#' ##    savaData$MarginCol
#' @export
bootMargin <-function(CC,  reps=10000,  pValue = 0.05, thr = 0.5, delete=FALSE){
  output <- wrapper.BootMargin(CC,  reps,  pValue, thr, delete)
  return(output)
}

