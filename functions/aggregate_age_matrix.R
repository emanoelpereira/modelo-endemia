#' aggregate_age_matrix.
#' 
#' A function to operate the average contacts with weighting between the age strucuture of a population.
#'
#' @param mat.cont N by N matrix with contacts between age bins 'mat.cont'.
#' @param aggregate_indices A list of vector of indexes of age bins that should be aggregated 'aggregate_indices'.
#' @param age_structure A vector of size N with population distribution to act as weights in a weighted average 'age_structure'.
#'
#' @return A new matrix with size \code{length(aggregate_index)} by \code{length(aggregate_index)} with average contacts weighted by population structure
#' @export
#'
#' @examples
aggregate_age_matrix <- function(mat.cont, aggregate_indices, age_structure){
  if (missing(age_structure))
  {
    age_structure = rep(c(1),each=dim(mat.cont)[1])
  }
  age_structure <- as.numeric(age_structure)
  if(dim(mat.cont)[1] != dim(mat.cont)[2]){
    stop("mat.cont is not a square matrix")
  }
  if(length(age_structure) < dim(mat.cont)[1]){
    stop("age_structure is smaller than the linear dimension of mat.cont")
  }
  if(length(age_structure) > dim(mat.cont)[1]){
    warning("age_structure is bigger than linear size of mat.cont, aggregating the last values in one")
    age_structure[dim(mat.cont)[1]] =  sum(age_structure[dim(mat.cont)[1]:length(age_structure)])
  }
  new_mat.cont = matrix(0,nrow=length(aggregate_indices),ncol=length(aggregate_indices))
  for (i in 1:length(aggregate_indices)){
    for(j in 1:length(aggregate_indices)){
      sum_age_struc = 0;
      for(k in aggregate_indices[[i]]){
        for(l in aggregate_indices[[j]]){
          new_mat.cont[i,j] = new_mat.cont[i,j] + age_structure[k]*mat.cont[k,l]
        }
        sum_age_struc = sum_age_struc + age_structure[k]
      }
      new_mat.cont[i,j] = new_mat.cont[i,j]/sum_age_struc
    }
  }
  return(new_mat.cont)
}
