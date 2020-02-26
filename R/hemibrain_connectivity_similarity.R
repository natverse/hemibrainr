#' Calculate a connectivity similarity score between two connectivity profiles
#'
#' @description Calculate a similarity score between connectivity matrices that penalises small differences between low and does not heavily penalise large differences between high weights. Algorithm from Jarrell et al. 2012.
#' @param x a vector/matrix of connectivities, where each entry in the vector or each column in the matrix is a different target/input neuron/cell type
#' @param y a different vector/matrix of connectivities
#' @param m an n x m adjacency matrix
#' @param c1 determines how negatively we want to punish a case such as the one above. Default C1 is chosen so that 1 and 5 are weakly dissimilar.
#' @param c2 determines the point where the similarity of the two numbers switches from negative to positive. Default C2 is chosen so that 10 and 100 synapses are weakly similar.
#' @param normalise perform a min-max normalisation on the similarity scores as in Schlegel et al. 2015
#' @param diag for connectivity_similarity_distance. Logical value indicating whether the diagonal of the distance matrix should be printed by print.dist.
#' @param upper for connectivity_similarity_distance. Logical value indicating whether the upper triangle of the distance matrix should be printed by print.dist.
#' @references Jarrell TA, Wang Y, Bloniarz AE, Brittin CA, Xu M, Thomson JN, Albertson DG, Hall DH, Emmons SW (2012) "The connectome of a decision-making neural network." Science (80- ) 337: 437–444.
#' @name hemibrain_connectivity_similarity
#' @export
#' @rdname hemibrain_connectivity_similarity
hemibrain_connectivity_similarity <- function(x,y, c1 = 0.5, c2 = 0.18, normalise = TRUE) UseMethod("hemibrain_connectivity_similarity")

#' @export
#' @rdname hemibrain_connectivity_similarity
hemibrain_connectivity_similarity.numeric <- function(x,y, c1 = 0.5, c2 = 0.18, normalise = TRUE){
  s = sapply(seq_along(x), function(i) min(x[i],y[i]) - c1*max(x[i],y[i])*exp(-c2*min(x[i],y[i])))
  score = sum(s)
  if(normalise){
    max.score = max(sapply(seq_along(x), function(i) max(x[i],y[i]) - c1*max(x[i],y[i])*exp(-c2*max(x[i],y[i]))))
    min.score = max(sapply(seq_along(x), function(i) - c1*max(x[i],y[i])))
    score = score/max.score(score-min.score)/(max.score - min.score)
    if(is.infinite(score)){
      score = 0
    }
  }
  score
}

#' @export
#' @rdname hemibrain_connectivity_similarity
hemibrain_connectivity_similarity.matrix<- function(x,y, c1 = 0.5, c2 = 0.18, normalise = TRUE){
  c = c()
  for(i in 1:ncol(x)){
    c = c(c,hemibrain_connectivity_similarity.numeric(x[,i],y[,i], normalise = normalise))
  }
  score = sum(c)
  if(normalise){
    score = score/ncol(x)
  }
  score
}

#' @export
#' @rdname hemibrain_connectivity_similarity
hemibrain_connectivity_similarity_distance <-function(m,c1 = 0.5, c2 = 0.18,normalise = FALSE, diag = FALSE, upper = FALSE){
  if(!is.matrix(m)){
    stop("m is not a matrix")
  }
  M = matrix(0,ncol(m),ncol(m))
  rownames(M) = colnames(M) = colnames(m)
  for(e in colnames(m)){
    for(i in colnames(m)){
      M[e,i] = hemibrain_connectivity_similarity.numeric(x = m[,e], y= m[,i],c1=c1,c2=c2,normalise = normalise)
    }
  }
  stats::dist(M,diag = diag,upper = upper,method = "euclidean")
}

#' @export
#' @rdname hemibrain_connectivity_similarity
hemibrain_connectivity_similarity_matrix <-function(m, c1 = 0.5, c2 = 0.18, normalise = FALSE){
  if(!is.matrix(m)&!is.data.frame(m)){
    stop("m is not a matrix")
  }
  M <- foreach::foreach(e = colnames(m), .combine = cbind) %dopar% {
    S <- foreach::foreach(i = colnames(m), .combine = rbind) %do% {
      hcs <- hemibrain_connectivity_similarity.numeric(x = m[,e], y= m[,i],c1=c1,c2=c2,normalise = normalise)
    }
  }
  rownames(M) = colnames(M) = colnames(m)
  M
}

# hidden
`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`
