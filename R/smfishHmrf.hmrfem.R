#' @title smfishHmrf.hmrfem
#' @name smfishHmrf.hmrfem
#' @description HMRF for single variate normal distribution
#' @param y gene expression matrix
#' @param neighbors adjacency matrix between cells
#' @param numnei an array of number of neighbors per cell
#' @param blocks a list of cell colors for deciding the order of cell update
#' @param beta the beta to try (smoothness parameter)
#' @param mu an array of cluster mean
#' @param sigma an array of cluster standard deviation
#' @param err the error that is allowed between successive iterations
#' @param maxit maximum number of iterations
#' @param verbose TRUE or FALSE
#' @return A list of prob, new mu, and new sigma after iterations finish
#' @keywords internal
smfishHmrf.hmrfem <- function(y, neighbors, numnei, blocks, beta=0.5, mu, sigma, err=1e-4, maxit=20, verbose){
    checkErrors(mu=mu, sigma=sigma, err=err)
    if (length(err) < 2) 
		err <- rep(err, length.out = 2)
    k <- length(mu)
    nvert <- length(y)
    maxnei <- ncol(neighbors)
    nblocks <- length(blocks)
    neighbors <- structure(as.integer(neighbors), dim = dim(neighbors))

    yunique <- sort(unique(y))
    n.yunique <- length(yunique)
    nvert <- length(y)
    ymatch <- match(y, yunique)

    muold <- rep(0,k)
    sigmaold <- rep(0,k)
    niter <- 10
    indices <- initialIndices(y, nvert, mu, sigma, k)

    it <- 0

    repeat{
        it <- it + 1
        muold <- mu
        sigmaold <- sigma
        for(j in 1:niter){
            den <- getDen(yunique, n.yunique, ymatch, mu, sigma)

            #indices <- updateIndicesHMRFEM(neighbors, numnei, maxnei, blocks, nblocks, beta, k, indices, check, den)
            indices <- updateIndicesHMRFEM(neighbors, numnei, maxnei, blocks, nblocks, beta, k, indices, den)
        }
        prob  <- updateProbabilities(neighbors, numnei, indices, den, k, beta)
        mu <- updateMeans(prob, y)
        sigma <- updateStdevs(prob, y, mu, nvert, k)

        flag <- checkStopVerbose(muold, mu, sigmaold, sigma, err, it, maxit, verbose)
        if(flag==1) break
    }

    list(prob=prob, mu=mu, sigma=sigma)
}

#' @name smfishHmrf.generate.centroid.use.exist
#' @title Use existing cluster centroids
#' @description This function assumes that cluster centroids have already been generated from previously applying kmeans on the dataset. The results should have been saved. It will load cluster centroids from existing clustering result files. The results should be found in `input_dir` directory. The function looks for the following two kmeans result files:
#' 1. `{input_dir}`/k_`{par_k}`/f`{name}`.gene.ALL.centroid.txt
#' 1. `{input_dir}`/k_`{par_k}`/f`{name}`.gene.ALL.kmeans.txt
#' 
#' where `{}` refers to the value of parameters down below
#' @param name name of this run
#' @param input_dir input directory
#' @param par_k number of clusters
#' @return A kmeans object which is a list with centers and cluster fields
#' @importFrom utils read.table
#' @importFrom fs path
#' @export
#' @examples
#' kmeans_results = system.file("extdata", package="smfishHmrf")
#' kk = smfishHmrf.generate.centroid.use.exist(name="test", input_dir=kmeans_results, par_k=9)
smfishHmrf.generate.centroid.use.exist <-function(name="test", input_dir=".", par_k){
	indir <- fs::path(input_dir, paste0("k_", par_k))
	centroid_file <- fs::path(indir, paste0("f", name, ".gene.ALL.centroid.txt"))
	kmeans_file <- fs::path(indir, paste0("f", name, ".gene.ALL.kmeans.txt"))
	centroid<-utils::read.table(centroid_file, header=F, row.names=1)
	centroid<-as.matrix(centroid)	
	kmeans<-utils::read.table(kmeans_file, header=F, row.names=1)
	kmeans<-as.matrix(kmeans)
	kk<-NULL
	kk$centers<-centroid
	kk$cluster<-kmeans
	kk
}

#' @name smfishHmrf.generate.centroid.it
#' @title Generate cluster centroids, where input is a file.
#' @description This function generates cluster centroids from applying kmeans. It accepts an expression matrix file as input. 
#' @param expr_file expression matrix file. The expression file should be a space-separated file. The rows are genes. The columns are cells. There is no header row. The first column is a gene index (ranges from 1 to the number of genes). Note the first column is not gene name.
#' @param par_k number of clusters
#' @param par_seed random generator seed (-1 if no fixing). Change the par_seed to vary the initialization.
#' @param nstart number of starts (kmeans). It is recommended to set nstart to at least 100 (preferrably 1000).
#' @param name name of this run
#' @param output_dir output directory; where to store the kmeans results
#' @return A kmeans object which is a list with centers and cluster fields
#' @importFrom fs path
#' @details
#' Note that after running kmeans step, the function also automatically saves the kmeans results to the `output_dir` directory. The results will be found in two files:
#' 1. `{output_dir}`/k_`{par_k}`/f`{name}`.gene.ALL.centroid.txt
#' 1. `{output_dir}`/k_`{par_k}`/f`{name}`.gene.ALL.kmeans.txt
#'
#' where `{}` refers to the value of parameters.
#' @examples 
#' mem_file = system.file("extdata", "ftest.expression.txt", package="smfishHmrf")
#' kk = smfishHmrf.generate.centroid.it(mem_file, par_k=9, par_seed=100, 
#'     nstart=100, name="test", output_dir=tempdir())
#' @export
smfishHmrf.generate.centroid.it <- function(expr_file, par_k, par_seed=-1, nstart, name="test", output_dir="."){
	y<-smfishHmrf.read.expression(expr_file)
	kk<-smfishHmrf.generate.centroid(y, par_k, par_seed, nstart)
	smfishHmrf.generate.centroid.save(kk, name=name, output_dir=output_dir)
	kk
}

#' @title smfishHmrf.read.expression
#' @name smfishHmrf.read.expression
#' @description Reads expression
#' @param expr_file expression matrix file. The expression file should be a space-separated file. The rows are genes. The columns are cells. There is no header row. The first column is a gene index (ranges from 1 to the number of genes). Note the first column is not gene name. 
#' @importFrom utils read.table
#' @return A matrix with gene expression matrix
#' @export
#' @examples 
#' expr_file = system.file("extdata", "ftest.expression.txt", package="smfishHmrf")
#' y<-smfishHmrf.read.expression(expr_file)
smfishHmrf.read.expression <- function(expr_file){
	y<-utils::read.table(expr_file, header=F, row.names=1)
	y<-as.matrix(y)
	y
}

#' @name smfishHmrf.generate.centroid.save
#' @title Save the cluster centroids to file
#' @description This function is run after the kmeans step. It takes a kmeans object (containing the kmeans result) as input and save the cluster centroids to file.
#' 
#' Note that the location of saving and the file names are decided by the following rule:
#' 1. `{output_dir}`/k_`{par_k}`/f`{name}`.gene.ALL.centroid.txt
#' 1. `{output_dir}`/k_`{par_k}`/f`{name}`.gene.ALL.kmeans.txt
#' 
#' where `{}` refers to the value of parameters.
#' @param kk kmeans object
#' @param name name of the run
#' @param output_dir output directory; where to save the results
#' @importFrom utils write.table
#' @export
#' @examples
#' expr_file = system.file("extdata", "ftest.expression.txt", package="smfishHmrf")
#' y<-smfishHmrf.read.expression(expr_file)
#' kk = smfishHmrf.generate.centroid(y, par_k=9, par_seed=100, nstart=100)
#' smfishHmrf.generate.centroid.save(kk, name="test", output_dir=tempdir())
#kk is kmeans clustering results
smfishHmrf.generate.centroid.save <- function(kk, name="test", output_dir="."){
	par_k <- dim(kk$centers)[1]
	savedir <- fs::path(output_dir, paste0("k_", par_k))
	dir.create(savedir, showWarnings=TRUE)
	centroid_file <- fs::path(savedir, paste0("f", name, ".gene.ALL.centroid.txt"))
	kmeans_file <- fs::path(savedir, paste0("f", name, ".gene.ALL.kmeans.txt"))
	utils::write.table(kk$cluster, file=kmeans_file, sep=" ", quote=F, col.names=F, row.names=T)
	utils::write.table(kk$centers, file=centroid_file, sep=" ", quote=F, col.names=F, row.names=T)
}

#' @name smfishHmrf.generate.centroid
#' @title Generate cluster centroids, where input is given as a matrix.
#' @description This function assumes that the input gene expression matrix file has been already loaded into a matrix. The function accepts a matrix and applies kmeans clustering to generate cluster centroids.
#' @param y expression matrix
#' @param par_k number of clusters
#' @param par_seed random generator seed (to fix it set it to above 0, or -1 if no fixing). Change the par_seed to vary the initialization.
#' @param nstart number of starts (kmeans parameter). It is recommended to set nstart to at least 100 (preferrably 1000).
#' @return A kmeans list with centers and cluster fields
#' @importFrom stats kmeans
#' @export
#' @examples
#' data(seqfishplus)
#' kk<-smfishHmrf.generate.centroid(seqfishplus$y, par_k=9, par_seed=100, nstart=100)
smfishHmrf.generate.centroid <- function(y, par_k, par_seed=-1, nstart){
	#par_k <- commandArgs(trailingOnly = TRUE)[1]
	#par_seed <- commandArgs(trailingOnly = TRUE)[2]
	#nstart <- commandArgs(trailingOnly = TRUE)[3]
	#mem_file <- commandArgs(trailingOnly = TRUE)[4]
	#centroid_file <- commandArgs(trailingOnly = TRUE)[5]
	#kmeans_file <- commandArgs(trailingOnly = TRUE)[6]
	#par_k <- as.integer(par_k)
	#par_seed <- as.integer(par_seed)
	#nstart <- as.integer(nstart)
	if(par_seed!=-1 & par_seed>0){
		set.seed(par_seed)
	}
	#y<-read.table(mem_file, header=F, row.names=1)
	y<-as.matrix(y)
	k<-par_k
	m<-dim(y)[2]
	kk<-stats::kmeans(y, k, nstart=nstart, iter.max=100)
	#if(write_file==TRUE){
	#	write.table(kk$cluster, file=kmeans_file, sep=" ", quote=F, col.names=F, row.names=T)
	#	write.table(kk$centers, file=centroid_file, sep=" ", quote=F, col.names=F, row.names=T)
	#}
	kk$centers<-as.matrix(kk$centers)
	kk$cluster<-as.matrix(kk$cluster)
	kk
}

#' @name smfishHmrf.hmrfem.multi.it.min
#' @title Perform HMRF on multivariate normal distributions. Accepts file names as inputs. Accepts multiple betas.
#' @description This function performs HMRF \insertCite{Zhu2018}{smfishHmrf} for multi variate normal distributions. It takes minimum required inputs (inputs being file names). There are a couple of files required:
#' 1. a file containing expression matrix
#' 1. a file containing cell neighborhood matrix 
#' 1. a file containing node (or cell) color. This is used for updating cells during HMRF iterations.
#' 
#' HMRF needs users to specify the initializations of parameters (mu and sigma). It is recommended to use the kmeans centroids as initializations (specified by `kk` parameter). Note: kmeans should be run prior to this function.
#'
#' 
#' @section Data preprocessing:
#' It assumes that the expression values follow a multivariate gaussian distribution. We generally recommend using **log2 transformed counts further normalized by z-scores (in both x- and y- dimensions)**. Double z-scoring this way helps to remove the inherent bias of zscoring just one dimension (as the results might present a bias towards cell counts).
#' @section Betas:
#' Beta is the smoothness parameter in HMRF. The higher the beta, the more the HMRF borrows information from the neighbors. This function runs HMRF across a range of betas. To decide which beta range, here are some guideline:
#' * if the number of genes is from 10 to 50, the recommended range is 0 to 10 at beta increment of 0.5.
#' * if the number of genes is below 50, the recommended range is 0 to 15 at beta increment of 1. 
#' * if the number of genes is between 50 to 100, the range is 0 to 50 at beta increment of 2. 
#' * if the number of genes is between 100 and 500, the range is 0 to 100 at beta increment of 5.
#'
#' Within the range of betas, we recommend selecting the best beta by the Bayes information criterion. This requires first performing randomization of spatial positions to generate the null distribution of log-likelihood scores for randomly distributed cells for the same range of betas. Then find the beta where the difference between the observed and the null log-likelihood is maximized. 
#' 
#' @section Variations:
#' - `smfishHmrf.hmrfem.multi.it.min` (this function): supports multiple betas; supports file names as inputs. Recommended.
#' - `smfishHmrf.hmrfem.multi.it`: supports multiple betas; supports R data structures as inputs.
#' - `smfishHmrf.hmrfem.multi`: supports a single beta; supports R data structures as inputs.
#'
#' @param mem_file expression file. The expression file should be a space-separated file. The rows are genes. The columns are cells. There is no header row. The first column is a gene index (ranges from 1 to the number of genes). Note the first column is not gene name. See section **Data preprocessing** for which form of expression works best.
#' @param nei_file file containing cell neighborhood matrix. This should be a space-separated file. The rows are cells. The columns are neighbors. There is no header row. The first column is the cell index (1 to number of cells). Each row lists the indices of neighbor cells. The dimension of the cell neighborhood matrix is (num_cell, max_num_neighbors). If a cell does not have enough neighbors, the remaining entries of that row is padded with -1. The R package Giotto <http://spatialgiotto.com> \insertCite{Dries701680}{smfishHmrf} contains a number of functions for generating the cell neighborhood network.
#' @param block_file file containing cell colors (which determines cell update order). The order of updating the state probabilities of each cell can matter the result. Cells (or nodes) and their immediate neighbors are not updated at the same time. This is akin to the vertex coloring problem. This file contains the color of each cell such that no two neighbor cells have the same color. The file is 2-column, space-separated. Column 1 is cell ID, and column 2 is the cell color (integer starting at 1). The python utility get_vertex_color.py <https://bitbucket.org/qzhudfci/smfishhmrf-py/src/master/get_vertex_color.py> (requires smfishHmrf-py package <https://pypi.org/project/smfishHmrf/>) can generate this file.
#' @param kk kmeans results (object returned by kmeans). Kmeans (one of functions smfishHmrf.generate.centroid.it or smfishHmrf.generate.centroid) should be run before this function.
#' @param par_k number of clusters
#' @param name name for this run (eg test)
#' @param output_dir output directory
#' @param tolerance tolerance
#' @param beta,beta_increment,beta_num_iter 3 values specifying the range of betas to try: the initial beta, the beta increment, and the number of betas. Beta is the smoothness parameter. Example: `beta`=0, `beta_increment`=2, `beta_num_iter`=6 means to try betas: 0, 2, 4, 6, 8, 10. See section **Betas** for more information.
#' @importFrom stats cov dist
#' @importFrom utils read.table
#' @export
#' @examples 
#' mem_file = system.file("extdata", "ftest.expression.txt", package="smfishHmrf")
#' nei_file = system.file("extdata", "ftest.adjacency.txt", package="smfishHmrf")
#' block_file = system.file("extdata", "ftest.blocks.txt", package="smfishHmrf")
#' par_k = 9
#' name = "test"
#' output_dir = tempdir()
#'     
#' \dontrun{
#' kk = smfishHmrf.generate.centroid.it(mem_file, par_k, par_seed=100, 
#' nstart=100, name=name, output_dir=output_dir)
#' }
#'
#' # alternatively, if you already have run kmeans before, you can load it directly
#' kmeans_results = system.file("extdata", package="smfishHmrf")
#' kk = smfishHmrf.generate.centroid.use.exist(name=name, input_dir=kmeans_results, par_k)
#'
#' smfishHmrf.hmrfem.multi.it.min(mem_file, nei_file, block_file, kk, par_k, 
#' name=name, output_dir=output_dir, tolerance=1e-5, 
#' beta=28, beta_increment=2, beta_num_iter=1)
#'     
#' \dontrun{
#' # alternatively, to test a larger set of beta's
#' smfishHmrf.hmrfem.multi.it.min(mem_file, nei_file, block_file, kk, par_k,
#' name=name, output_dir=output_dir, tolerance=1e-5, 
#' beta=0, beta_increment=2, beta_num_iter=20)
#' }
#'  
#' @references
#' \insertRef{Zhu2018}{smfishHmrf}
#'
#' \insertRef{Eng2019}{smfishHmrf}
#'
#' \insertRef{Dries701680}{smfishHmrf}
#'
smfishHmrf.hmrfem.multi.it.min <- function(mem_file, nei_file, block_file, kk, par_k, 
name="test", output_dir=".", tolerance=1e-5, beta=0, beta_increment=1, beta_num_iter=10){
	#file reading
	y<-utils::read.table(mem_file, header=F, row.names=1)
	y<-as.matrix(y)
	nei<-utils::read.table(nei_file, header=F, row.names=1)
	colnames(nei)<-NULL
	rownames(nei)<-NULL
	nei<-as.matrix(nei)
	blocks<-utils::read.table(block_file, header=F, row.names=1)
	blocks<-c(t(blocks))
	maxblock <- max(blocks)
	blocks<-lapply(1:maxblock, function(x) which(blocks == x))
	numnei<-apply(nei, 1, function(x) sum(x!=-1))
	#centroid<-read.table(centroid_file, header=F, row.names=1)
	#centroid<-as.matrix(centroid)
	centroid<-kk$centers
	#parameter setting
	k<-par_k
	m<-dim(y)[2]
	sigma <-array(0, c(m,m,k))
	for(i in 1:k){
		sigma[, ,i] <- cov(y)
		print(rcond(sigma[,,i]))
	}
	mu<-array(0, c(m,k))
	kk2<-centroid
	for(i in 1:k){
		mu[,i] <- kk2[i,]
	}
	numcell<-dim(y)[1]
	kk_dist<-array(0, c(numcell, k))
	for(i in 1:numcell){
		for(j in 1:k){
			kk_dist[i,j] <- stats::dist(rbind(y[i,], mu[,j]), method="euclidean")
		}
	}
	clust_mem<-apply(kk_dist, 1, function(x) which(x==min(x))) 
	lclust<-lapply(1:k, function(x) which(clust_mem == x))
	damp<-array(0, c(k))
	for(i in 1:k){
		sigma[, , i] <- stats::cov(y[lclust[[i]], ])
		#default tolerance is 1e-60
		di<-findDampFactor(sigma[,,i], factor=1.05, d_cutoff=tolerance, startValue=0.0001)
		if(is.null(di)){
			damp[i] = 0
		}else{
			damp[i] = di
		}
	}

	smfishHmrf.hmrfem.multi.it(name, output_dir, k, y, nei, beta=beta, beta_increment=beta_increment, beta_num_iter=beta_num_iter, 
	numnei, blocks, mu, sigma, damp)
}

#' @name smfishHmrf.hmrfem.multi.save
#' @title Save the HMRF result
#' @description This function assumes that HMRF has been run via smfishHmrf.hmrfem.multi, smfishHmrf.hmrfem.multi.it or smfishHmrf.hmrfem.multi.it.min function. It assumes the results have been generated. This function saves the results of each beta to the output directory. It will return void.
#' @param name name for this run (eg test)
#' @param outdir output directory
#' @param beta beta to save
#' @param tc.hmrfem the result of running of hmrfem on single beta (from smfishHmrf.hmrfem.multi)
#' @param k number of clusters
#' @importFrom utils write.table
#' @export
#' @examples
#' data(seqfishplus)
#' s <- seqfishplus
#' tc.hmrfem<-smfishHmrf.hmrfem.multi(s$y, s$nei, s$numnei, s$blocks, beta=28,
#'     mu=s$mu, sigma=s$sigma, err=1e-7, maxit=50, verbose=TRUE, dampFactor=s$damp,
#'     tolerance=1e-5)
#' smfishHmrf.hmrfem.multi.save(name="test", outdir=tempdir(), beta=28, tc.hmrfem=tc.hmrfem, k=9)
#needs y, nei, beta, numnei, blocks, mu, sigma, damp
smfishHmrf.hmrfem.multi.save <- function(name, outdir, beta, tc.hmrfem, k){
	savedir <- fs::path(outdir, paste0("k_", k))
	dir.create(savedir, showWarnings=TRUE)
	out_file <- fs::path(savedir, sprintf("f%s.%.1f.prob.txt", name, beta)) #hmrfem probability
	out_file_2 <- fs::path(savedir, sprintf("f%s.%.1f.centroid.txt", name, beta)) #hmrfem centroids
	out_file_3 <- fs::path(savedir, sprintf("f%s.%.1f.hmrf.covariance.txt", name, beta)) #hmrfem covariance
	out_file_unnorm <- gsub("prob", "unnormprob", out_file)

	utils::write.table(tc.hmrfem$prob, file=out_file, sep=" ", quote=F, col.names=F, row.names=T)
	utils::write.table(tc.hmrfem$unnormprob, file=out_file_unnorm, sep=" ", quote=F, col.names=F, row.names=T)
	utils::write.table(t(tc.hmrfem$mu), file=out_file_2, sep=" ", quote=F, col.names=F, row.names=T)
	utils::write.table(tc.hmrfem$sigma[,,1], file=out_file_3, sep=" ", quote=F, col.names=F, row.names=T)
	for(i in 2:k){
		utils::write.table(tc.hmrfem$sigma[,,i], file=out_file_3, sep=" ", quote=F, col.names=F, row.names=T, append=T)
	}
}

#' @name smfishHmrf.hmrfem.multi.it
#' @title Perform HMRF for multivariate normal distribution. Accepts R data structures as inputs. Accepts multiple betas.
#' @description This function performs HMRF model \insertCite{Zhu2018}{smfishHmrf} on inputs which are directly R data structures. Different from smfishHmrf.hmrfem.multi, this function iterates over multiple betas rather than a single beta. Different from smfishHmrf.hmrfem.multi.it.min, this function accepts R data structures (i.e. parameters y, nei, blocks) as inputs to the function rather than accepting file names. This function will save the results of HMRF to the output directory. It will return void. 
#'
#' This function exists for legacy and compatibility reason. User should use **smfishHmrf.hmrfem.multi.it.min** function.
#' @param name name for this run (eg test)
#' @param outdir output directory
#' @param k number of clusters
#' @param y gene expression matrix
#' @param nei adjacency matrix between cells
#' @param beta initial beta
#' @param beta_increment beta increment
#' @param beta_num_iter number of betas to try
#' @param numnei a vector containing number of neighbors per cell
#' @param blocks a list of cell colors for deciding the order of cell update
#' @param mu a 2D matrix (i,j) of cluster mean (initialization) 
#' @param sigma a 3D matrix (i,j,k) where (i,j) is the covariance of cluster k (initialization)
#' @param damp a list of dampening factors (length = k)
#'
#' @section More information:
#' Arguments mu and sigma refer to the cluster centroids from running kmeans algorithm. 
#' They serve as initialization of HMRF.
#' Users should refer to **smfishHmrf.hmrfem.multi.it.min** for more information about function parameters and the requirements.
#'
#' @references
#' \insertRef{Zhu2018}{smfishHmrf}
#'
#' @export
#' @examples 
#' y<-as.matrix(read.table(system.file("extdata", "ftest.expression.txt", 
#'     package="smfishHmrf"), header=FALSE, row.names=1))
#' nei<-as.matrix(read.table(system.file("extdata", "ftest.adjacency.txt", 
#'     package="smfishHmrf"), header=FALSE, row.names=1))
#' colnames(nei)<-NULL; rownames(nei)<-NULL
#' blocks<-c(t(read.table(system.file("extdata", "ftest.blocks.txt", 
#'     package="smfishHmrf"), header=FALSE, row.names=1)))
#' blocks<-lapply(1:max(blocks), function(x) which(blocks == x))
#' numnei<-apply(nei, 1, function(x) sum(x!=-1))
#' k<-9
#' kmeans_results = system.file("extdata", package="smfishHmrf")
#' kk = smfishHmrf.generate.centroid.use.exist(name="test", input_dir=kmeans_results, k)
#' numcell<-dim(y)[1]; m<-dim(y)[2]
#' mu<-t(kk$centers) #should be dimension (m,k)
#' lclust<-lapply(1:k, function(x) which(kk$cluster == x))
#' damp<-array(0, c(k)); sigma<-array(0, c(m,m,k))
#' for(i in 1:k){
#'     sigma[, , i] <- cov(y[lclust[[i]], ])
#'     di<-findDampFactor(sigma[,,i], factor=1.05, d_cutoff=1e-5, startValue=0.0001)
#'     damp[i]<-ifelse(is.null(di), 0, di)
#' } 
#' smfishHmrf.hmrfem.multi.it(name="test", outdir=tempdir(), k=k, y=y, nei=nei, 
#'     beta=28, beta_increment=2, beta_num_iter=1, numnei=numnei, blocks=blocks, 
#'     mu=mu, sigma=sigma, damp=damp)
#'
#' \dontrun{
#' # alternatively, to test a larger set of betas:
#' smfishHmrf.hmrfem.multi.it(name="test", outdir=tempdir(), k=k, y=y, nei=nei, 
#'     beta=0, beta_increment=2, beta_num_iter=20, numnei=numnei, blocks=blocks, 
#'     mu=mu, sigma=sigma, damp=damp)
#' }
#'
#'
smfishHmrf.hmrfem.multi.it <- function(name, outdir, k, y, nei, beta=0, beta_increment=1, beta_num_iter=10, 
numnei, blocks, mu, sigma, damp){
	beta_current <- beta
	res <- c()
	for(bx in 1:beta_num_iter){
		print(sprintf("Doing beta=%.3f", beta_current))
		tc.hmrfem<-smfishHmrf.hmrfem.multi(y=y, neighbors=nei, beta=beta_current, numnei=numnei, 
		blocks=blocks, mu=mu, sigma=sigma, verbose=T, err=1e-7, maxit=50, dampFactor=damp)
		smfishHmrf.hmrfem.multi.save(name, outdir, beta_current, tc.hmrfem, k)
		#do_one(name, outdir, k, y, nei, beta_current, numnei, blocks, mu, sigma, damp)
		t_key <- sprintf("k=%d b=%.2f", k, beta_current)
		res[[t_key]] <- tc.hmrfem
		beta_current <- beta_current + beta_increment
	}
	res
}

#' @name smfishHmrf.hmrfem.multi
#' @title Performs HMRF for multivariate normal distribution. Accepts R data structures as inputs. Accepts a single beta.
#' @description This function performs HMRF \insertCite{Zhu2018}{smfishHmrf} on multivariate normal distributions. Different from other variations, this function accepts R data structures directly as inputs, and only accepts a single value of beta.
#'
#' This function exists for legacy and compatibility reason. User should use **smfishHmrf.hmrfem.multi.it.min** function.
#' @param y gene expression matrix
#' @param neighbors adjacency matrix between cells
#' @param numnei a vector containing number of neighbors per cell
#' @param blocks a list of cell colors for deciding the order of cell update
#' @param beta the beta to try (smoothness parameter)
#' @param mu a 2D matrix (i,j) of cluster mean (initialization) 
#' @param sigma a 3D matrix (i,j,k) where (i,j) is the covariance of cluster k (initialization)
#' @param err the error that is allowed between successive iterations
#' @param maxit maximum number of iterations
#' @param dampFactor the dampening factor
#' @param forceDetectDamp will auto detect a dampening factor instead of using the specified one
#' @param tolerance applicable when forceDetectDamp is set to TRUE
#' @param verbose TRUE or FALSE
#' @return A list of prob, new mu, new sigma, unnormalized prob after iterations finish
#' @importFrom pracma std eps
#' @useDynLib smfishHmrf, .registration=TRUE
#'
#' @section More information:
#' Arguments mu and sigma refer to the cluster centroids from running kmeans algorithm. 
#' They serve as initialization of HMRF.
#' Users should refer to **smfishHmrf.hmrfem.multi.it.min** for more information about function parameters and the requirements.
#'
#' @references
#' \insertRef{Zhu2018}{smfishHmrf}
#'
#' @export
#' @examples
#'     data(seqfishplus)
#'     s <- seqfishplus
#'     res<-smfishHmrf.hmrfem.multi(s$y, s$nei, s$numnei, s$blocks, beta=28, 
#'     mu=s$mu, sigma=s$sigma, err=1e-7, maxit=50, verbose=TRUE, dampFactor=s$damp, 
#'     tolerance=1e-5)
smfishHmrf.hmrfem.multi <- function(y, neighbors, numnei, blocks, beta=0.5, mu, sigma, err=1e-7, maxit=50, verbose, dampFactor=NULL, forceDetectDamp=FALSE, tolerance=1e-60){
    checkErrorsMulti(mu=mu, sigma=sigma, err=err)
    if (length(err) < 2)
		err <- rep(err, length.out = 2)

	k <- dim(mu)[2] #mu is a m * k matrix
	nvert <- dim(y)[1] #number of data points, n
	numdim <- dim(y)[2] #number of dimensions, m
    maxnei <- ncol(neighbors)
    nblocks <- length(blocks)
    neighbors <- structure(as.integer(neighbors), dim = dim(neighbors))
    muold <- array(0, dim(mu)) # m * k matrix
    sigmaold <- array(0, dim(sigma)) # m * m * k matrix 
    niter <- 10

    indices <- initialIndicesMulti(y, mu, sigma, k, dampFactor=dampFactor, forceDetectDamp=forceDetectDamp)

    it <- 0
	eps <- pracma::eps()
	very_small <- mean(pracma::std(y)) * (1/k) * 0.0001

    repeat{
        it <- it + 1
        muold <- mu
        sigmaold <- sigma
		den = matrix(0, nrow=nvert, ncol=k)
		for(j in 1:k){
			determinant <- det(sigma[,,j])
			invsigma <- array(0, dim(sigma[,,j]))
			#choice 1======================
			if(forceDetectDamp==FALSE & !is.null(dampFactor)){
				sigma[,,j]<-sigma[,,j] + diag(numdim) * dampFactor[j]
				invsigma <- solve(sigma[,,j])
				determinant <- det(sigma[,,j])
			}
			else if(forceDetectDamp==TRUE){
				damp<-findDampFactor(sigma[,,j], factor=1.05, d_cutoff=tolerance, startValue=0.0001)
				if(!is.null(damp)){
					sigma[,,j]<-sigma[,,j] + diag(numdim) * damp
					invsigma <- solve(sigma[,,j])
					determinant <- det(sigma[,,j])
					#print(damp)
				}
			}
			dist <- y - matrix(rep(t(mu[,j]),nvert),nrow=nvert,ncol=numdim,byrow=T)
			exponent <- -0.5 * (rowSums((dist %*% invsigma) * dist))
			const <- (numdim/2)*log(2*pi) + (1/2)*log(determinant)
			exp_diff <- exponent - const
			lim_exp_diff <- sapply(exp_diff, function(x) min(700, max(x, -700)))
			den[,j] <- exp(lim_exp_diff)
		}
        for(j in 1:niter){
            indices <- updateIndicesHMRFEM(neighbors, numnei, maxnei, blocks, nblocks, beta, k, indices, den)
        }
		unnorm_prob <- updateUnnormProbabilitiesMulti(neighbors, numnei, indices, den, k, beta)
        prob  <- updateProbabilitiesMulti(neighbors, numnei, indices, den, k, beta)
        mu <- updateMeansMulti(prob, y)
        sigma <- updateCovariancesMulti(prob, y, mu, nvert, k)

        flag <- checkStopVerboseMulti(muold, mu, sigmaold, sigma, err, it, maxit, verbose)

        if(flag==1) break
    }

	class_id <- max.col(prob, "first")

	t_tot <- 0
	t_prob <- rep(0, dim(unnorm_prob)[2])
	for(j in 1:dim(unnorm_prob)[1]){
		t_tot_density <- 0
		for(k in 1:dim(unnorm_prob)[2]){
			t_tot_density <- t_tot_density + unnorm_prob[j,k]
		}
		t_this_prob <- rep(0, dim(unnorm_prob)[2])
		for(k in 1:dim(unnorm_prob)[2]){
			t_this_prob[k] <- unnorm_prob[j,k] / t_tot_density
			t_prob[k] = t_prob[k] + t_this_prob[k]
		}
		t_tot <- t_tot + 1
	}
	for(k in 1:dim(unnorm_prob)[2]){
		t_prob[k] <- t_prob[k] / t_tot
	}
	t_likelihood <- 0
	for(j in 1:dim(unnorm_prob)[1]){
		t_this_like <- 0
		for(k in 1:dim(unnorm_prob)[2]){
			t_this_like <- t_this_like + t_prob[k] * unnorm_prob[j,k]
		}
		t_likelihood <- t_likelihood + log(t_this_like)
	}

	#t_likelihood <- t_likelihood / t_tot
    list(prob=prob, mu=mu, sigma=sigma, unnormprob=unnorm_prob, class=class_id, likelihood=t_likelihood)
}

