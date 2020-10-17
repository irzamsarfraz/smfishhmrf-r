
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
    indices <- initialIndices(y, nvert, mu, sigma, k, sub=FALSE)

    it <- 0

    repeat{
        it <- it + 1
        muold <- mu
        sigmaold <- sigma
        for(j in 1:niter){
            den <- getDensity(yunique, n.yunique, ymatch, mu, sigma)
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

smfishHmrf.generate.centroid.use.exist <-function(name="test", input_dir=".", par_k){
	indir <- fs::path(input_dir, paste0("k_", par_k))
	centroid_file <- fs::path(indir, paste0("f", test, ".gene.ALL.centroid.txt"))
	kmeans_file <- fs::path(indir, paste0("f", test, ".gene.ALL.kmeans.txt"))
	centroid<-read.table(centroid_file, header=F, row.names=1)
	centroid<-as.matrix(centroid)	
	kmeans<-read.table(kmeans_file, header=F, row.names=1)
	kmeans<-as.matrix(kmeans)
	kk<-NULL
	kk$centers<-centroid
	kk$cluster<-kmeans
	kk
}

smfishHmrf.generate.centroid.it <- function(expr_file, par_k, par_seed=-1, nstart, name="test", output_dir="."){
	y<-smfishHmrf.read.expression(expr_file)
	kk<-smfishHmrf.generate.centroid(y, par_k, par_seed, nstart)
	savedir <- fs::path(output_dir, paste0("k_", par_k))
	dir.create(savedir, showWarnings=TRUE)
	centroid_file <- fs::path(savedir, paste0("f", test, ".gene.ALL.centroid.txt"))
	kmeans_file <- fs::path(savedir, paste0("f", test, ".gene.ALL.kmeans.txt"))
	smfishHmrf.generate.centroid.save(kk, centroid_file, kmeans_file)
	kk
}

smfishHmrf.read.expression <- function(expr_file){
	y<-read.table(expr_file, header=F, row.names=1)
	y<-as.matrix(y)
	y
}

#kk is kmeans clustering results
smfishHmrf.generate.centroid.save <- function(kk, centroid_file, kmeans_file){
	write.table(kk$cluster, file=kmeans_file, sep=" ", quote=F, col.names=F, row.names=T)
	write.table(kk$centers, file=centroid_file, sep=" ", quote=F, col.names=F, row.names=T)
}

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
	kk<-kmeans(y, k, nstart=nstart, iter.max=100)
	#if(write_file==TRUE){
	#	write.table(kk$cluster, file=kmeans_file, sep=" ", quote=F, col.names=F, row.names=T)
	#	write.table(kk$centers, file=centroid_file, sep=" ", quote=F, col.names=F, row.names=T)
	#}
	kk$centers<-as.matrix(kk$centers)
	kk$cluster<-as.matrix(kk$cluster)
	kk
}

smfishHmrf.hmrfem.multi.it.min <- function(mem_file, nei_file, block_file, kk, par_k, 
name="test", output_dir=".", tolerance=1e-5, beta=0, beta_increment=1, beta_num_iter=10){
	#file reading
	y<-read.table(mem_file, header=F, row.names=1)
	y<-as.matrix(y)
	nei<-read.table(nei_file, header=F, row.names=1)
	colnames(nei)<-NULL
	rownames(nei)<-NULL
	nei<-as.matrix(nei)
	blocks<-read.table(block_file, header=F, row.names=1)
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
			kk_dist[i,j] <- dist(rbind(y[i,], mu[,j]), method="euclidean")
		}
	}
	clust_mem<-apply(kk_dist, 1, function(x) which(x==min(x))) 
	lclust<-lapply(1:k, function(x) which(clust_mem == x))
	damp<-array(0, c(k))
	for(i in 1:k){
		sigma[, , i] <- cov(y[lclust[[i]], ])
		#default tolerance is 1e-60
		di<-findDampFactor(sigma[,,i], factor=1.05, d_cutoff=tolerance, startValue=0.0001)
		if(is.null(di)){
			damp[i] = 0
		}else{
			damp[i] = di
		}
	}

	smfishHmrf.hmrfem.multi.it(name, outdir, k, y, nei, beta=beta, beta_increment=beta_increment, beta_num_iter=beta_num_iter, 
	numnei, blocks, mu, sigma, damp)
}

#needs y, nei, beta, numnei, blocks, mu, sigma, damp
smfishHmrf.hmrfem.multi.save <- function(name, outdir, beta, tc.hmrfem, k){
	out_file <- sprintf("%s/%s.%.1f.prob.txt", outdir, name, par_beta) #hmrfem probability
	out_file_2 <- sprintf("%s/%s.%.1f.centroid.txt", outdir, name, par_beta) #hmrfem centroids
	out_file_3 <- sprintf("%s/%s.%.1f.hmrf.covariance.txt", outdir, name, par_beta) #hmrfem covariance
	out_file_unnorm <- gsub("prob", "unnormprob", out_file)

	write.table(tc.hmrfem$prob, file=out_file, sep=" ", quote=F, col.names=F, row.names=T)
	write.table(tc.hmrfem$unnormprob, file=out_file_unnorm, sep=" ", quote=F, col.names=F, row.names=T)
	write.table(t(tc.hmrfem$mu), file=out_file_2, sep=" ", quote=F, col.names=F, row.names=T)
	write.table(tc.hmrfem$sigma[,,1], file=out_file_3, sep=" ", quote=F, col.names=F, row.names=T)
	for(i in 2:k){
		write.table(tc.hmrfem$sigma[,,i], file=out_file_3, sep=" ", quote=F, col.names=F, row.names=T, append=T)
	}
}

smfishHmrf.hmrfem.multi.it <- function(name, outdir, k, y, nei, beta=0, beta_increment=1, beta_num_iter=10, 
numnei, blocks, mu, sigma, damp){
	beta_current <- beta
	for(bx in 1:beta_num_iter){
		print(sprintf("Doing beta=%.3f", beta_current))
		tc.hmrfem<-smfishHmrf.hmrfem.multi(y=y, neighbors=nei, beta=beta_current, numnei=numnei, 
		blocks=blocks, mu=mu, sigma=sigma, verbose=T, err=1e-7, maxit=50, dampFactor=damp)
		smfishHmrf.hmrfem.multi.save(name, outdir, beta_current, tc.hmrfem, k)
		#do_one(name, outdir, k, y, nei, beta_current, numnei, blocks, mu, sigma, damp)
		beta_current <- beta_current + beta_increment
	}
}

smfishHmrf.hmrfem.multi <- function(y, neighbors, numnei, blocks, beta=0.5, mu, sigma, err=1e-4, maxit=20, verbose, dampFactor=NULL, forceDetectDamp=FALSE){
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
	eps <- eps()
	very_small <- mean(std(y)) * (1/k) * 0.0001

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
				damp<-findDampFactor(sigma[,,j], factor=1.05, d_cutoff=1e-60, startValue=0.0001)
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

		#need to check this!
        flag <- checkStopVerboseMulti(muold, mu, sigmaold, sigma, err, it, maxit, verbose)

        if(flag==1) break
    }
    list(prob=prob, mu=mu, sigma=sigma, unnormprob=unnorm_prob)
}

