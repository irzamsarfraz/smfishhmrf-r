#' SeqFISHplus dataset
#'
#' Data from SeqFISH experiment on SS cortex. This is a dataset with 523 cells and the expression of about 500 spatial genes
#' 
#' @docType data
#'
#' @usage data(seqfishplus)
#'
#' @format A list containing the following fields: y, nei, blocks, damp, mu, sigma
#' \describe{
#' \item{y}{gene expression matrix}
#' \item{nei}{cell adjacency matrix}
#' \item{blocks}{vertex (or cell) update order; a list of vertex colors; cells marked with the same color are updated at once}
#' \item{damp}{dampening constants (length k, the number of clusters)}
#' \item{mu}{initialization (means). Means is a (i,k) matrix}
#' \item{sigma}{initialization (sigmas). Sigmas is a (i,j,k) 3D matrix. k is cluster id. (i,j) is covariance matrix}
#' }
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{Eng2019}{smfishHmrf}
#' @examples
#' data(seqfishplus)
"seqfishplus"
