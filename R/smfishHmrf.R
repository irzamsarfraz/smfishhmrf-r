#' smfishHmrf: A package for running hidden markov random field on smFISH and other spatial transcriptomic datasets
#' 
#' A package for running hidden markov random field \insertCite{Zhu2018}{smfishHmrf} on smFISH and other spatial transcriptomic datasets.
#'
#' @section Input:
#' The inputs of HMRF are the following:
#' - Gene expression matrix
#' - Cell neighborhood matrix
#' - Initial centroids of clusters
#' - Number of clusters
#' - beta
#' 
#' smfishHmrf has been tested to work on seqFISH, MERFISH, starMAP, 10X Visium and other datasets. See Giotto \insertCite{Dries701680}{smfishHmrf} for examples of such datasets and to learn about the technologies. 
#' smfishHmrf is a general algorithm, and should probably work with other data types.
#'
#' @section Running:
#' The first step is to calculate initial centroids on the gene expression matrix given k (the number of clusters). The function **smfishHmrf.generate.centroid.it** is used for this purpose. 
#'
#' The next step is to run the HMRF algorithm given the expression matrix, and cell neighborhood matrix. The function **smfishHmrf.hmrfem.multi.it.min** is used for this purpose.
#'
#' @section Variations:
#' You might notice several variations of the functions:
#' - `smfishHmrf.hmrfem.multi.it.min`: supports multiple betas; supports file names as inputs. **This is the recommended function.**
#' - `smfishHmrf.hmrfem.multi.it`: supports multiple betas; supports R data structures as inputs.
#' - `smfishHmrf.hmrfem.multi`: supports a single beta; supports R data structures as inputs.
#' - Note: beta is the smoothness parameter of HMRF
#'
#' Also:
#' - `smfishHmrf.generate.centroid.it`: supports file names as inputs. **This is the recommended function** 
#' - `smfishHmrf.generate.centroid`: supports R matrices as inputs. Assumes input files have been read into R matrices.
#' - `smfishHmrf.generate.centroid.use.exist`: loads existing centroids. Assumes that centroids have been generated previously and saved to disk.
#'
#' @docType package
#' @name smfishHmrf
#' @references
#' \insertRef{Zhu2018}{smfishHmrf}
#'
#' \insertRef{Dries701680}{smfishHmrf}
#'
NULL
