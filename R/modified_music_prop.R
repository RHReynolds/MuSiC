#' MuSiC Deconvolution
#'
#' This function is to calculate the MuSiC deconvolution proportions
#'
#' @param bulk.eset ExpressionSet for bulk data
#' @param sc.eset ExpressionSet for single cell data
#' @param markers Vector or list of gene names, default as NULL. If NULL, use
#'   all genes that provided by both bulk and single cell dataset.
#' @param clusters Character, the phenoData of single cell dataset used as
#'   clusters;
#' @param samples Character,the phenoData of single cell dataset used as
#'   samples;
#' @param select.ct Vector of cell types, default as NULL. If NULL, then use all
#'   cell types provided by single cell dataset;
#' @param ct.cov Logical. If TRUE, use the covariance across cell types;
#' @param verbose Logical, default as TRUE.
#' @param iter.max Numeric, maximum iteration number
#' @param nu Regulation parameter, take care of weight when taking reciprocal
#' @param eps Thredshold of convergence
#' @param centered Logic, substract avg of Y and D
#' @param normalize Logic, divide Y and D by their standard deviation
#' @param celltype.batching Logical. If TRUE, no longer require one
#'   expressionSet for the entire dataset. Instead, dataset can be subsetted by
#'   cell type into a number of expressionSets. I.e. Each expressionSet contains
#'   all samples and all genes commonly expressed (non-zero expression) across
#'   all samples, with only a few cell types from the entire dataset. This works
#'   well for large datasets that have > 50,000 cells/nuclei in total.
#' @param celltype.eSet.dir Character. File path to directory containing all
#'   cell-type subsetted expressionSets for an entire dataset. Each
#'   expressionSet should be saved as an .Rds object.
#'
#' @return a list with elements: * Estimates of MuSiC * Estimates of NNLS *
#'   Weight of MuSiC * r.squared of MuSiC * Variance of MuSiC estimates
#' @seealso \code{\link{music_basis}}
#' @export

music_prop_modified = function(bulk.eset, sc.eset = NULL, markers = NULL, clusters, samples, select.ct = NULL, ct.cov = FALSE, verbose = TRUE,
                      iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE, 
                      celltype.batching = FALSE, celltype.eSet.dir = NULL, ... ){
  
  # Additional libraries needed
  library(purrr)
  library(qdapTools)
  library(tidyverse)
  library(stringr)
  library(xbioc)
  
  # Remove genes with non-zero expression in bulk data
  bulk.gene = rownames(bulk.eset)[rowMeans(exprs(bulk.eset)) != 0]
  bulk.eset = bulk.eset[bulk.gene, , drop = FALSE]
  
  if(is.null(markers)){
    sc.markers = bulk.gene
  }else{
    sc.markers = intersect(bulk.gene, unlist(markers))
  }
  
  # Cell type batching argument
  if(celltype.batching == FALSE){
    
    sc.basis = music_basis_modified(sc.eset, non.zero = TRUE, markers = sc.markers, clusters = clusters, samples = samples, select.ct = select.ct, ct.cov = ct.cov, verbose = verbose)
    
  } else{
    
    if(verbose){message(paste('Getting MuSiC basis from batched dataset.'))}
    
    if(is.null(celltype.eSet.dir)){
      stop("No directory provided for cell-type batching.")
    }
    
    eSet_file_paths <- list.files(path = celltype.eSet.dir, full.names = TRUE, pattern = ".Rds")
    
    if(length(eSet_file_paths) < 1){
      stop("No files with extension .Rds found. Was the correct directory provided?")
    }
    
    # For each file, load dataset and perform music_basis()
    for(i in 1:length(eSet_file_paths)){
      
      if(verbose){message(paste('Loading expressionSet from:', eSet_file_paths[i]))}
      
      sc.eSet <- readRDS(eSet_file_paths[i])
      
      # Set non-zero to FALSE as removal of genes with non-zero expression occurs across all samples in default mode
      # Thus snRNA-seq must be filtered prior to entry into music_basis()
      basis_subset <- music_basis_modified(sc.eSet, non.zero = FALSE, markers = sc.markers, clusters = clusters, samples = samples, select.ct = select.ct, ct.cov = ct.cov, verbose = verbose)
      
      if(i == 1){
        
        sc.basis <- basis_subset
        
      } else{
        
        # Need to add conditional statement for if only one cell type in sc.eSet, as might occur when subsetting dataset into pairs of cell types
        # music_basis_modified will return a numeric vector, thus must convert to matrix.
        # Two exceptions: (1) M.S., which must remain a named vector and (2) S is returned as a matrix
        if(length(pVar(sc.eSet, clusters) %>% unique()) == 1){
          
          for(j in 1:length(basis_subset)){
            
            # Exception 1
            if(names(basis_subset[j]) == "M.S"){
              
              sc.basis[[j]] <- c(sc.basis[[j]], basis_subset[[j]])
              
            } else {
              
              # Exception 2
              if(class(basis_subset[[j]]) == "matrix"){
                
                ct_name <- pVar(sc.eSet, clusters) %>% unique() %>% as.character()
                
                # Reorder genes in basis_subset by order in sc.basis
                basis_subset[[j]] <- basis_subset[[j]][match(rownames(sc.basis[[j]]), rownames(basis_subset[[j]])),] %>% 
                  as.matrix()
                
                colnames(basis_subset[[j]]) <- ct_name
                
                sc.basis[[j]] <- sc.basis[[j]] %>%
                  cbind(basis_subset[[j]])
                
              } else {
                
                ct_name <- pVar(sc.eSet, clusters) %>% unique() %>% as.character()
                
                # Reorder genes in basis_subset by order in sc.basis
                basis_subset[[j]] <- basis_subset[[j]][match(rownames(sc.basis[[j]]), names(basis_subset[[j]]))] %>% 
                  as.matrix()
                
                colnames(basis_subset[[j]]) <- ct_name
                
                sc.basis[[j]] <- sc.basis[[j]] %>%
                  cbind(basis_subset[[j]])
                
              }
              
            }
            
          }
          
        } else{
          
          for(j in 1:length(basis_subset)){
            
            if(names(basis_subset[j]) == "M.S"){
              
              sc.basis[[j]] <- c(sc.basis[[j]], basis_subset[[j]])
              
            } else {
              
              # Reorder genes in basis_subset by order in sc.basis
              basis_subset[[j]] <- basis_subset[[j]][match(rownames(sc.basis[[j]]), rownames(basis_subset[[j]])),] 
              
              sc.basis[[j]] <- sc.basis[[j]] %>%
                cbind(basis_subset[[j]])
              
            }
            
          }
          
        }
        
      }
      
    }
    
    if(verbose){message(paste('Full MuSiC basis loaded'))}
    
  }
  
  cm.gene = intersect( rownames(sc.basis$Disgn.mtx), bulk.gene )
  
  if(is.null(markers)){
    if(length(cm.gene)< 0.2*min(length(bulk.gene), nrow(sc.basis$Disgn.mtx)) )
      stop("Too few common genes!")
  }else{
    if(length(cm.gene)< 0.2*length(unlist(markers)))
      stop("Too few common genes!")
  }
  if(verbose){message(paste('Used', length(cm.gene), 'common genes...'))}
  
  m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx)); m.bulk = match(cm.gene, bulk.gene)
  D1 = sc.basis$Disgn.mtx[m.sc, ]; M.S = colMeans(sc.basis$S, na.rm = T);
  Yjg = relative.ab(exprs(bulk.eset)[m.bulk, ]); N.bulk = ncol(bulk.eset);
  if(ct.cov){
    Sigma.ct = sc.basis$Sigma.ct[, m.sc];
    
    if(sum(Yjg[, i] == 0) > 0){
      D1.temp = D1[Yjg[, i]!=0, ];
      Yjg.temp = Yjg[Yjg[, i]!=0, i];
      Sigma.ct.temp = Sigma.ct[, Yjg[,i]!=0];
      if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...') )
    }else{
      D1.temp = D1;
      Yjg.temp = Yjg[, i];
      Sigma.ct.temp = Sigma.ct;
      if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...'))
    }
    
    lm.D1.weighted = music.iter.ct(Yjg.temp, D1.temp, M.S, Sigma.ct.temp, iter.max = iter.max,
                                   nu = nu, eps = eps, centered = centered, normalize = normalize)
    Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
    Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
    weight.gene.temp = rep(NA, nrow(Yjg)); weight.gene.temp[Yjg[,i]!=0] = lm.D1.weighted$weight.gene;
    Weight.gene = cbind(Weight.gene, weight.gene.temp)
    r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
    Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
  }else{
    Sigma = sc.basis$Sigma[m.sc, ];
    
    valid.ct = (colSums(is.na(Sigma)) == 0)&(colSums(is.na(D1)) == 0)&(!is.na(M.S))
    
    if(sum(valid.ct)<=1){
      stop("Not enough valid cell types!")
    }
    
    if(verbose){message(paste('Used', sum(valid.ct), 'cell types in deconvolution...' ))}
    
    D1 = D1[, valid.ct]; M.S = M.S[valid.ct]; Sigma = Sigma[, valid.ct];
    
    Est.prop.allgene = NULL
    Est.prop.weighted = NULL
    Weight.gene = NULL
    r.squared.full = NULL
    Var.prop = NULL
    for(i in 1:N.bulk){
      if(sum(Yjg[, i] == 0) > 0){
        D1.temp = D1[Yjg[, i]!=0, ];
        Yjg.temp = Yjg[Yjg[, i]!=0, i];
        Sigma.temp = Sigma[Yjg[,i]!=0, ];
        if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...') )
      }else{
        D1.temp = D1;
        Yjg.temp = Yjg[, i];
        Sigma.temp = Sigma;
        if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...'))
      }
      
      lm.D1.weighted = music.iter(Yjg.temp, D1.temp, M.S, Sigma.temp, iter.max = iter.max,
                                  nu = nu, eps = eps, centered = centered, normalize = normalize)
      Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp = rep(NA, nrow(Yjg)); weight.gene.temp[Yjg[,i]!=0] = lm.D1.weighted$weight.gene;
      Weight.gene = cbind(Weight.gene, weight.gene.temp)
      r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
  }
  colnames(Est.prop.weighted) = colnames(D1)
  rownames(Est.prop.weighted) = colnames(Yjg)
  colnames(Est.prop.allgene) = colnames(D1)
  rownames(Est.prop.allgene) = colnames(Yjg)
  names(r.squared.full) = colnames(Yjg)
  colnames(Weight.gene) = colnames(Yjg)
  rownames(Weight.gene) = cm.gene
  colnames(Var.prop) = colnames(D1)
  rownames(Var.prop) = colnames(Yjg)
  
  return(list(Est.prop.weighted = Est.prop.weighted, Est.prop.allgene = Est.prop.allgene,
              Weight.gene = Weight.gene, r.squared.full = r.squared.full, Var.prop = Var.prop))
}

#' Prepare Design matrix and Cross-subject Variance for MuSiC Deconvolution
#'
#' This function is used for generating cell type specific cross-subject mean and variance for each gene. Cell type specific library size is also calcualted.
#'
#' @param x ExpressionSet, single cell dataset
#' @param non.zero logical, default as TRUE. If true, remove all gene with zero expression.
#' @param markers vector or list of gene names. Default as NULL. If NULL, then use all genes provided.
#' @param clusters character, the phenoData used as clusters;
#' @param samples character,the phenoData used as samples;
#' @param select.ct vector of cell types. Default as NULL. If NULL, then use all cell types provided.
#' @param ct.cov logical. If TRUE, use the covariance across cell types.
#' @param verbose logical, default as TRUE.
#' @return a list of
#'     * gene by cell type matrix of Design matrix
#'     * subject by celltype matrix of Library size
#'     * vector of average library size for each cell type
#'     * gene by celltype matrix of average relative abudance
#'     * gene by celltype matrix of cross-subject variation
#'
#' @export
music_basis_modified = function(x, non.zero = TRUE, markers = NULL, clusters, samples, select.ct = NULL, ct.cov = FALSE, verbose = TRUE){
  
  if(!is.null(select.ct)){
    s.ct = sampleNames(x)[as.character(pVar(x, clusters)) %in% select.ct]
    x <- x[, s.ct, drop = FALSE]
  }
  
  # This would need to be performed across all subjects and cell types within a single cell dataset to ensure common genes were the same
  # This is the only operation that is not cell-type specific -- but could be batched by gene or possible even performed before eSets were 
  # put through music_basis()
  if(non.zero){  ## eliminate non expressed genes
    nz.gene = rownames(x)[( rowSums(exprs(x)) != 0 )]
    x <- x[nz.gene, , drop = FALSE]
  }
  
  # # If markers used (e.g. all genes in bulk data), filter single-cell data for intersecting genes
  # # This is in lieu of original filtering that occurred following calculation of Sigma.ct, Sigma, M.theta, D
  # if (!is.null(markers)){
  #   ids <- intersect(unlist(markers), rownames(x))
  #   m.ids = match(ids, rownames(x))
  #   x <- x[m.ids, ]
  # }
  
  # Extract cell type and sample values from phenotype data in expressionSet 
  clusters <- as.character(pVar(x, clusters))
  samples <- as.character(pVar(x, samples))
  
  # Use BiocGenerics sapply, which applies function over list-like object
  # For each cell type, 'ct', calculate proportion of expression assigned to one gene
  # I.e. divide each gene's expression by expression of all genes
  # Then average proportion of a gene's expression across all samples
  # Here it is possible to split the dataset by cell types
  M.theta <- sapply(unique(clusters), function(ct){
    my.rowMeans(sapply(unique(samples), function(sid){
      y = exprs(x)[,clusters %in% ct & samples %in% sid, drop = FALSE]
      rowSums(y)/sum(y)
    }), na.rm = TRUE)
  })
  if(verbose){message("Creating Relative Abudance Matrix...")}
  
  # Creating covariance/variance matrix
  if(ct.cov){
    # Only relevant if covariance matrix required
    # Works in a similar manner to variance sapply i.e. functions applied to each cell type independently of another
    nGenes = nrow(x);
    n.ct = length(unique(clusters));
    nSubs = length(unique(samples))
    
    Theta <- sapply(unique(clusters), function(ct){
      sapply(unique(samples), function(sid){
        y = exprs(x)[,clusters %in% ct & samples %in% sid, drop = FALSE]
        return( rowSums(y)/sum(y) )
      })
    })
    if(!is.null(select.ct)){
      m.ct = match(select.ct, colnames(Theta))
      Theta = Theta[, m.ct]
    }
    
    Sigma.ct = sapply(1:nGenes, function(g){
      sigma.temp = Theta[nGenes*(0:(nSubs - 1)) + g, ];
      Cov.temp = cov(sigma.temp)
      Cov.temp1 = cov(sigma.temp[rowSums(is.na(Theta[nGenes*(0:(nSubs - 1)) + 1, ])) == 0, ])
      Cov.temp[which(colSums(is.na(sigma.temp))>0), ] = Cov.temp1[which(colSums(is.na(sigma.temp))>0), ]
      Cov.temp[, which(colSums(is.na(sigma.temp))>0)] = Cov.temp1[, which(colSums(is.na(sigma.temp))>0)]
      return(Cov.temp)
    })
    colnames(Sigma.ct) = rownames(x);
    
    # Code to intersect with marker genes e.g. all bulk genes.
    # Positioning here means Sigma.ct calculated prior to intersection with bulk genes
    if (!is.null(markers)){
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      Sigma.ct <- Sigma.ct[ , m.ids]
    }
    if(verbose){message("Creating Covariance Matrix...")}
  }else{
    # Use BiocGenerics sapply, which applies function over list-like object
    # For each cell type, 'ct', calculate proportion of expression assigned to one gene
    # I.e. divide each gene's expression by expression of all genes
    # Then calculate variance of a gene's expression across all samples
    # Here it is possible to split the dataset by cell types
    Sigma <- sapply(unique(clusters), function(ct){
      apply(sapply(unique(samples), function(sid){
        y = exprs(x)[,clusters %in% ct & samples %in% sid, drop = FALSE]
        rowSums(y)/sum(y)
      }), 1, var, na.rm = TRUE)
    })
    if(!is.null(select.ct)){
      m.ct = match(select.ct, colnames(Sigma))
      Sigma = Sigma[, m.ct]
    }
    
    # Code to intersect with marker genes e.g. all bulk genes.
    # Positioning here means Sigma calculated prior to intersection with bulk genes
    if (!is.null(markers)){
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      Sigma <- Sigma[m.ids, ]
    }
    if(verbose){message("Creating Variance Matrix...")}
  }
  
  # Use BiocGenerics sapply, which applies function over list-like object
  # For each cell type, 'ct', sum expression of all genes within a sample and divide by number of cells of cell type 'ct' within sample 
  # I.e. library size across each cell type
  # Here it is possible to split the dataset by cell types
  S <- sapply(unique(clusters), function(ct){
    my.rowMeans(sapply(unique(samples), function(sid){
      y = exprs(x)[, clusters %in% ct & samples %in% sid, drop = FALSE]
      sum(y)/ncol(y)
    }), na.rm = TRUE)
  })
  if(verbose){message("Creating Library Size Matrix...")}
  
  # For each cell type, calculate average library size across all samples
  S[S == 0] = NA
  M.S = colMeans(S, na.rm = TRUE)
  #S.ra = relative.ab(S, by.col = FALSE)
  #S.ra[S.ra == 0] = NA
  #S[S == 0] = NA
  #M.S = mean(S, na.rm = TRUE)*ncol(S)*colMeans(S.ra, na.rm = TRUE)
  
  # Design matrix
  # For each cell type, multiply each gene's proportion of expression (M.theta) by average library size (M.S)
  D <- t(t(M.theta)*M.S)
  
  if(!is.null(select.ct)){
    m.ct = match(select.ct, colnames(D))
    D = D[, m.ct]
    S = S[, m.ct]
    M.S = M.S[m.ct]
    M.theta = M.theta[, m.ct]
  }
  
  # Code to intersect with marker genes e.g. all bulk genes.
  # Positioning here means D and M.theta calculated prior to intersection with bulk genes
  if (!is.null(markers)){
    ids <- intersect(unlist(markers), rownames(x))
    m.ids = match(ids, rownames(x))
    D <- D[m.ids, ]
    M.theta <- M.theta[m.ids, ]
  }
  
  # Return design matrix D, together with S, M.S, M.theta, and Sigma matrices
  if(ct.cov){
    return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, Sigma.ct = Sigma.ct))
  }else{
    return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, Sigma = Sigma))
  }
}
