
CountriesList  <- c("Afganistan","Albania","Alemania","Algeria","Andorra","Angola","Anguilla","Antigua y Barbuda","Antillas Holandesas"
                    ,"Arabia Saudita","Argentina"       ,"Armenia","Aruba","Australia","Austria","Azerbaiyan","Bahamas","Bahrein"
                    ,"Bangladesh","Barbados","Belgica","Belice","Benin","Bermudas","Bielorrusia","Bolivia"
                    ,"Bosnia y Herzegovina","Botsuana","Brasil","Brunei","Bulgaria","Burkina Faso","Burundi","Butan",
                    "Cabo Verde","Camboya","Camerun","Canada","Chad","Chile","China","Chipre","Colombia","Comores",
                    "Congo","Cook, Islas","Corea del Norte","Corea del Sur","Costa de Marfil","Costa Rica","Croacia","Dinamarca","Ecuador",
                    "Egipto","El Salvador","Emiratos Arabes Unidos","Eritrea","Eslovaquia","Eslovenia", "Espa?a","Estados Unidos","Estonia","Etiopia","Feroe, Islas","Filipinas","Finlandia","Fiyi","Francia")
CountriesList <- toupper(CountriesList)
length(CountriesList)
CountriesList <-chartr("??????", "AEIOUN", toupper(CountriesList))
CreateCountries <- function(CountriesList, nFilas, nColumnas, nexp ){
  datos <- array(rep(0,nFilas*nColumnas*nexp), dim = c(nFilas,nColumnas,nexp))
  for(i in 1:nFilas){
    for(j in 1:nColumnas){
      if(i!=j){
        media_ij<-runif(1,0,1)
        sd_ij<-rgamma(1,1,10)
        rnd<-runif(nexp,0,1)
        ijk<-qnorm(rnd,mean = media_ij,sd = sd_ij)
        ijk[which(ijk < 0)] <- 0
        ijk[which(ijk > 1)] <- 1
        ijk2 <- as.numeric(formatC( ijk, format="f", digits = 17))
        ijk3 <- as.character(round(ijk2, 1))

        datos[i,j, ] <-   as.numeric(ijk3)
      } else{
        datos[i,j, ] <- 1
      }
    }
  }
  return(datos)
}
#' @title createDataSet
#'
#' @param rows Number of rows
#' @param columns Number of columns
#' @param experts Number of experts
#'
#' @return returns a list containing 3 three-dimensional matrices
#' @export
#'
createDataSet <- function( rows, columns, experts){
  AA <-AA <- CreateCountries(CountriesList,rows,columns,experts)
  rownames(AA) <- CountriesList[1:nrow(AA)]
  colnames(AA) <- CountriesList[1:ncol(AA)]

  AB <- CreateCountries(CountriesList, rows, (columns-1), experts)
  rownames(AB) <- CountriesList[1:(nrow(AB))]
  colnames(AB) <- CountriesList[(rows+2):(ncol(AB)+(rows+1))]


  BB <- CreateCountries(CountriesList, (rows-1),(columns-1), experts)
  rownames(BB) <- CountriesList[(rows+2):(nrow(BB)+(rows+1))]
  colnames(BB) <- CountriesList[(rows+2):(ncol(BB)+(rows+1))]
  return(list(AA=AA,AB=AB,BB=BB))
}
# m <-createDataSet(15,15,15)
# m$AA
# m$AB
# m$BB
