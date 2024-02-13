#' Compute a correlation matrix an return a summary value of the values
#'
#' @param x the dataset on which the correlation should be computed
#' @param method the correlation method to apply (default: "spearman")
#' @param sf the summary function (default: median)
#' @param square TRUE if the correlation values should be squared
#' (default: FALSE)
#' @param byrow if TRUE (default) compute the correlation between rows
#'
#' @return The result of the sf function applied on the correlation
#' lower triangular matrix, excluding diagonal
#'
#' @export
#'
corSummary <- function(
    x, method="spearman", sf=median, square=FALSE, byrow=TRUE
){
  stopifnot(ncol(x)>0, nrow(x)>0)
  if(byrow){
    x <- t(x)
  }
  cv <- as.dist(cor(x, method=method))
  if(square){
    cv <- cv^2
  }
  toRet <- sf(cv)
  return(toRet)
}


##############################################################@
#' Differential co-expression
#'
#' identify differentially coexpression modules comparing two conditions
#' using a permutation test across samples
#'
#' @param sampleTable data.frame with two columns, sample_id and cohort
#' @param contrasts character vector specifying contrasts (e.g. Disease-Control)
#' @param moduleList BEIDlist object containing the modules
#' @param expression expression matrix containing the genes used to construct the modules
#' @param corMeth a character string indicating which correlation coefficient is to
#' be computed
#' @param perm  positive integer corresponding to the number of permutations to be done in
#' order to assess the significance of the results (default: 100)
#' @param mc.cores the number of cores to use during the permutation procedure
#' (default: 5)
#'
#'
#' @return The results of the differential coexpression
#' permutation test. tibble with colnames: Module, R2 cohort 1, R2 cohort 2,
#' Difference_Median, pvalue
#'
#' @export
#'
dcTest <- function(
    sampleTable,
    contrasts,
    moduleList,
    expression,
    corMeth = "spearman",
    perm = 100,
    mc.cores = 5
){

  print("checks")
  stopifnot(is.BEIDList(moduleList))
  stopifnot(all(unlist(moduleList) %in% rownames(expression)))
  stopifnot(all(sampleTable$sample_id %in% colnames(expression)))

  contrasts <- unlist(strsplit(contrasts, split = "-"))

  ## load clusters
  print("load clusters")
  scope <- c("sampleTable", "expression","corMeth",
             "moduleList", "contrasts")
  cl <- parallel::makeCluster(mc.cores)
  on.exit({
    parallel::stopCluster(cl)
  })
  parallel::clusterExport(cl, scope, envir=environment())

  ## randomization
  print("randomization")
  rDC <- do.call(cbind,
                 parallel::parLapply(cl, 1:perm, function(i) {
                   randInd <- sampleTable
                   randInd$cohort <- sample(randInd$cohort)
                   rdiv0 <- cor(t(expression[,
                                             randInd$sample_id[randInd$cohort == contrasts[2]]]),
                                method = corMeth)
                   rdiv1 <- cor(t(expression[,
                                             randInd$sample_id[randInd$cohort == contrasts[1]]]),
                                method = corMeth)
                   toAdd <- unlist(lapply(
                     moduleList,
                     function(x){median(as.dist(rdiv1[unlist(x),unlist(x)]^2)) -
                         median(as.dist(rdiv0[unlist(x),unlist(x)]^2))
                     }))
                   rm(rdiv1, rdiv0)
                   return(toAdd)
                 }))

  ##
  print("cormat")
  corMat0 <- cor(t(expression[,
                              sampleTable$sample_id[sampleTable$cohort == contrasts[2]]]),
                 method=corMeth)
  corMat1 <- cor(t(expression[,
                              sampleTable$sample_id[sampleTable$cohort == contrasts[1]]]),
                 method=corMeth)

  r20 <- sapply(moduleList, function(x){
    median(as.dist(corMat0[x, x]^2))
  })
  r2 <- sapply(moduleList, function(x){
    median(as.dist(corMat1[x,x]^2))
  })
  DC <-   r2 - r20

  ##
  print("permutation")
  ## When R² < 0 -> test whether sign. more values show lower R²
  t <- setdiff(names(DC[DC < 0]), NA)
  pval <- apply(cbind(DC[t],rDC[t,]),1,
                function(x){
                  sum(x <= x[1])/length(x)
                })
  ## When R² > 0 -> test whether sign. more values show higher R²
  t <- setdiff(names(DC[DC > 0]), NA)
  pval <- c(pval, apply(cbind(DC[t],rDC[t,]),1,
                        function(x){
                          sum(x >= x[1])/length(x)
                        }))

  nm1 <- paste("R2", contrasts[2], sep = "_")
  nm2 <- paste("R2", contrasts[1], sep = "_")
  toRet <- tibble::tibble("Module" = names(DC),
                          !!nm1 := r20[names(DC)],
                          !!nm2 := r2[names(DC)],
                          "Difference_Median" = round(DC,digits = 3),
                          "pvalue" = pval[names(DC)])
  return(toRet)
}

#' Grouping list items according to shared elements
#'
#'
#' @param x a list of item associated elements
#' @param sharing average percentage of elements shared by the items
#' (default: 0.5).
#' @param ... parameters for the [lhclust] function.
#'
#' @return A numeric vector providing group index of each term.
#'
#' @export
#'
lgrouping <- function(
    x,
    sharing=0.5,
    ...
){
  hc <- lhclust(x, ...)
  groups <- cutree(hc, h=1-sharing)
  attr(groups, "hc") <- hc
  return(groups)
}

#' Monte Carlo significance of a dataset summary value by feature list
#'
#'
#' @param l the list of feature sets
#' @param d a dataset with rownames corresponding to all the features in l
#' @param f a function taking a dataset (subset of d) as an input and returning
#' a single numerical summary value
#' @param alternative the alternative to be tested: should the actual value
#' be "greater" (default) or "less" than the values based on permuations
#' @param replace should be TRUE if l is redundant and FALSE (default) if l
#' is a dichotomy
#' @param perm a positive integer corresponding to the number
#' of permutations to be done in order to assess the significance of the
#' results (default: 100)
#' @param mc.cores the number of cores to use during the permutation procedure
#' (default: 4)
#' @param flibs a character vector with names of necessary libraries
#' (default: NULL)
#' @param ... additional parameters to pass to f
#'
#' @return A tibble with the following fields:
#' - **name**: the name of the feature set
#' - **n**: the length of the feature set
#' - **value**: the value of the statistics of interest
#' - **p.value**: the significance of the value based on permutations
#' - **FDR**: Benjamini-Hochberg False Discovery Rate
#'
#' @export
#'
lmcTest <- function(
    l, d, f,
    alternative=c("greater", "less"),
    replace=FALSE,
    perm=100,
    mc.cores=4,
    flibs=NULL,
    ...
){

  ############################################################################@
  ## Checks ----
  alternative <- match.arg(alternative)
  greater <- alternative=="greater"
  perm <- as.integer(perm)
  stopifnot(length(perm)==1, !is.na(perm), perm>0)
  mc.cores <- as.integer(mc.cores)
  stopifnot(length(mc.cores)==1, !is.na(mc.cores), mc.cores>0)

  stopifnot(all(unlist(l) %in% rownames(d)))

  ############################################################################@
  ## Compute actual values ----
  mval <- unlist(lapply(
    l,
    function(x){
      toRet <- f(d[x,,drop=FALSE], ...)
      stopifnot(length(toRet)==1)
      return(toRet)
    }
  ))

  ############################################################################@
  ## Compute permuation values ----
  ulibs <- flibs
  ln <- unlist(lapply(l, length))
  mnvec <- do.call(c, lapply(
    names(ln),
    function(n) rep(n, ln[[n]])
  ))
  scope <- c("l", "d", "f", "replace", "mnvec", "ulibs")
  cl <- makeCluster(mc.cores, type="PSOCK")
  on.exit({
    stopCluster(cl)
  })
  clusterExport(cl, scope, envir=environment())
  if(length(ulibs)>0){
    clusterEvalQ(cl, lapply(ulibs, library, character.only = TRUE))
  }
  rval <- do.call(cbind, parLapply(
    cl,
    1:perm,
    function(i){
      toRet <- unlist(lapply(
        split(sample(rownames(d), length(mnvec), replace=replace), mnvec),
        function(x){
          f(d[x, ,drop=FALSE], ...)
        }
      ))
      return(toRet)
    }
  ))

  ############################################################################@
  ## Compute p-values ----
  if(greater){
    pval <- apply(
      cbind(mval, rval[names(mval),]),
      1,
      function(x){
        (sum(x[-1]>=x[1])+1)/length(x)
      }
    )
  }else{
    pval <- apply(
      cbind(mval, rval[names(mval),]),
      1,
      function(x){
        (sum(x[-1]<=x[1])+1)/length(x)
      }
    )
  }
  tval <- tibble(
    "l"=names(mval),
    "n"=ln,
    "value"=mval,
    "p.value"=pval,
    "FDR"=p.adjust(pval, method="BH")
  )
  return(tval)

}


#' Enrichment between 2 lists of vectors of ID
#'
#' @author Patrice Godard (\email{patrice.godard@@ucb.com})
#'
#' @param query a named list of sets of IDs of interest
#' @param reference a named list of sets of IDs of reference
#' @param omega the ID universe for enrichment analysis. If NULL (default),
#' all the IDs in query or in reference
#' @param mc.cores number of cores used to parallelize the analysis (default: 4)
#'
#' @return A tibble with the following fields:
#' - **q**: the set of IDs of interest
#' - **qsize**: the number of IDs in q considered for the enrichment
#' - **r**: the reference term
#' - **rsize**: the number of IDs in r considered for the enrichment
#' - **i**: the length of the intersection between q and r
#' - **P-Value**: the P-Value returned by the hypergeometric test
#' - **FDR**: False Discovery Rate (Benjamini & Hochberg (1995)  method)
#'
#' @export
#'
#'
qrlEnrich <- function(query, reference, omega=NULL, mc.cores=4){
  if(length(omega)==0){
    omega <- union(
      unique(unlist(query)),
      unique(unlist(reference))
    )
  }
  query <- lapply(query, intersect, omega)
  reference <- lapply(reference, intersect, omega)

  query <- lapply(
    names(query),
    function(l){
      toRet <- query[[l]]
      attr(toRet, "name") <- l
      return(toRet)
    }
  )
  cl <- makeCluster(mc.cores)
  on.exit({
    stopCluster(cl)
  })
  # f <- lenrich
  # scope <- "f"
  # clusterExport(cl, scope, envir=environment())
  clusterEvalQ(cl, {
    test <- try(library(TBTools), silent=TRUE)
    n <- 1
    while(inherits(test, "try-error") & n < 5){
      test <- try(library(TBTools), silent=TRUE)
      n <- n+1
    }
    if(inherits(test, "try-error")){
      stop(test)
    }
  })
  toRet <- parLapply(
    cl,
    query,
    lenrich,
    reference, omega
  )

  toRet <- do.call(rbind, toRet)

  toRet <- toRet %>%
    as_tibble %>%
    select("q", "qsize", "r", "rsize", "i", "P-Value", "FDR")
  return(toRet)
}

#' Enrichment of a list of vectors of ID in a set of ID
#'
#' @param x a vector of ID
#' @param reference a named list of sets of IDs of reference
#' @param omega the ID universe for enrichment analysis.
#'
lenrich <- function(x, reference, omega){
  lsize <- length(x)
  toRet <- do.call(rbind, lapply(
    reference,
    function(y){
      il <- length(intersect(x, y))
      return(data.frame(
        rsize=length(y),
        i=il,
        "P-Value"=phyper(
          q=il-1,
          m=length(x),
          n=length(omega)-length(x),
          k=length(y),
          lower.tail=F
        ),
        stringsAsFactors=FALSE,
        check.names=FALSE
      ))
    }
  ))
  toRet$FDR <- p.adjust(toRet$`P-Value`, method="BH")
  toRet$r <- rownames(toRet)
  toRet$q <- attr(x, "name")
  toRet$qsize <- lsize
  toRet <- as_tibble(toRet)
  return(toRet)
}


msCutree <- function(tree, k=NULL, h=NULL, minsize, d=NULL, corMeth="spearman"){
  clusters <- cutree(tree, k=k, h=h)
  if(minsize<=1){
    return(clusters)
  }
  clustSize <- table(clusters)
  if(min(clustSize) >= minsize){
    return(clusters)
  }
  if(sum(clustSize >= minsize) < 2){
    toRet <- rep(1, length(clusters))
    names(toRet) <- names(clusters)
    return(toRet)
  }
  toKeep <- names(clustSize)[which(clustSize >= minsize)]
  toMerge <- names(clustSize)[which(clustSize < minsize)]
  if(length(toMerge)==0){
    return(clusters)
  }
  if(is.null(d)){
    stop("d should be provided for merging clusters")
  }
  eg <- getEigenValues(clusters, d)
  egCor <- cor(t(eg), method=corMeth)
  egCor <- egCor[toKeep, toMerge, drop=FALSE]^2
  toRet <- clusters
  for(i in toMerge){
    selClust <- toKeep[which.max(egCor[,i])]
    toRet[which(toRet==as.numeric(i))] <- selClust
  }
  return(toRet)
}

treeCutQual <- function(
    tree,
    corMat,
    k.min=1,
    k.max=200,
    minsize=1,
    d,
    corMeth="spearman",
    BPPARAM=SerialParam()
){
  toRet <- do.call(rbind, bplapply(
    k.min:k.max,
    function(k){
      modules <- msCutree(tree, k=k, minsize=minsize, d=d, corMeth=corMeth)
      qc <- clQual(modules, corMat=corMat)
      toAdd <- data.frame(
        k=k,
        n=nrow(qc),
        "Size median"=median(qc$size),
        "R2 weighted median"=sum(qc$r2med * qc$size)/sum(qc$size),
        "R2 median"=median(qc$r2med),
        check.names=F
      )
      return(toAdd)
    },
    BPPARAM=BPPARAM
  ))
  return(toRet)
}

plotInfl <- function(infl){
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  par(mfrow=c(3,1))
  par(mar=c(0.5, 4.1, 0.5, 1.1))
  plot(
    infl[,1:2],
    main="",
    xlab="",
    ylab=colnames(infl[2]),
    xaxt="n"
  )
  points(infl[,1], infl$py, col="red", type="l")
  abline(v=attr(infl, "xlim"), col="blue", lty=2)
  text(
    x=attr(infl, "xlim"),
    y=min(infl[,2]),
    attr(infl, "xlim"),
    col="blue",
    pos=4
  )
  ##
  par(mar=c(4.1, 4.1, 0.5, 1.1))
  plot(
    infl[,c(colnames(infl)[1], "dy")],
    main="",
    xlab=colnames(infl)[1],
    ylab="dy/dx"
  )
  points(infl[,1], infl$pdy, col="red", type="l")
  abline(v=attr(infl, "xlim"), col="blue", lty=2)
  ##
  par(mar=c(0.5, 4.1, 0.5, 1.1))
  plot(
    infl[,c(colnames(infl)[1], "ddy")],
    main="",
    xlab="",
    ylab="d(dy/dx)/dx",
    xaxt="n"
  )
  points(infl[,1], infl$pddy, col="red", type="l")
  abline(v=attr(infl, "xlim"), col="blue", lty=2)
  abline(h=0, col="grey", lty=3, lwd=3)
}

getSubModules <- function(corMat, logFC){
  stopifnot(!is.null(names(logFC)))

  modTree <- hclust(as.dist(1-corMat))
  # subMod <- cutree(modTree, h=1)
  # if(length(unique(subMod))>2){
  #   subMod <- cutree(modTree, k=2)
  # }
  subMod <- cutree(modTree, k=2)
  imCor <- median(
    corMat[names(subMod)[which(subMod==1)], names(subMod)[which(subMod==2)]],
    na.rm=TRUE
  )
  if(imCor > 0){
    subMod <- list("1"=rownames(corMat), "2"=c())
  }else{
    subMod <- split(names(subMod), subMod)
  }
  subModDeg <- lapply(
    subMod,
    function(x){
      logFC[intersect(names(logFC), x)]
    }
  )
  smdAvg <- unlist(lapply(subModDeg, mean, na.rm=TRUE))
  smdAvg <- ifelse(is.na(smdAvg), 0, smdAvg)
  if(!is.na(smdAvg[2])){
    if(smdAvg[1] <= smdAvg[2]){
      smdNames <- c("u", "o")
    }else{
      smdNames <- c("o", "u")
    }
  }else{
    if(smdAvg[1] < 0){
      smdNames <- c("u", "o")
    }else{
      smdNames <- c("o", "u")
    }
  }
  names(subMod) <- names(subModDeg) <- smdNames
  return(subMod[c("o","u")])
}

clQual <- function(clusters, corMat){

  ##################################
  ## Internal functions
  ##################################

  ## Clusters to list
  clToList <- function(clusters){
    clusters <- data.frame(
      a=names(clusters),
      b=clusters,
      stringsAsFactors=F
    )
    clusters <- by(clusters, clusters$b, function(d) d$a)
    class(clusters) <- "list"
    return(clusters)
  }

  ## Quality by cluster
  cQual <- function(cluster, corMat){
    n <- length(cluster)
    if(n < 2){
      r2med <- 1
      r2mad <- 0
    }else{
      if(n > 1000){
        scluster <- sample(cluster, 1000, replace=F)
      }else{
        scluster <- cluster
      }
      ## Cluster quality
      corMat <- corMat[scluster, scluster]
      r2 <- as.dist(corMat^2)
      r2med <- median(r2)
      r2mad <- mad(r2)
    }
    return(data.frame(size=n, r2med=r2med, r2mad=r2mad))
  }

  ##################################
  ##################################

  clusters <- clToList(clusters)
  return(do.call(rbind, lapply(
    clusters,
    cQual,
    corMat=corMat
  )))

}



