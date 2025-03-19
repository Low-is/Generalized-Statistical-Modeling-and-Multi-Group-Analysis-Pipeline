#### Extract meta data for RNA-seq files using ARCHS4 ####
### [Optional]: Can also extract meta data for RNA-seq studies using this
#### For connecting python with R ####
### library(reticulate)
### use_python("/env/bin/python3", required = TRUE)
### py_config()
### py_module_available("archs4py")
### a4 <- import("archs4py")
### file = "/app/data/ARCHS4/human_gene_v2.5.h5"
### a4$ls(file)

# pData_RNA <- lapply(rna_studies, function( geo ) {
#   file = "/app/data/ARCHS4/human_gene_v2.5.h5"
#   a4$meta$series(file, geo)
# })

#### Generating expression matrices ####
generate_exprs_mtx <- function(DNA = NULL, RNA = NULL, dna_studies = list(), rna_mtx = list(), pData_RNA = list()) {
  
  #### Handling duplicated genes ####
  handling_duplicates <- function(data_table) {
    
    if (is.list(data_table$Genes)) {
      data_table[, Genes := vapply(Genes, function(x) {
        # Check if the list is empty or NULL
        if (is.null(x) || length(x) == 0 || all(is.na(x))) {
          return(NA_character_)  # Return NA for empty or NULL elements
        } else if (is.list(x)) {
          return(as.character(unlist(x)))  # Unlist nested lists
        } else {
          return(as.character(x))  # Directly return non-list elements
        }
      }, FUN.VALUE = character(1))] # Ensures results is a character vector
    }
    
    data_table <- na.omit(data_table[(Genes != "" & !is.na(Genes)), ])
    
    if (any(duplicated(data_table$Genes))) {
      exprs_data.table <- data_table[, lapply(.SD, mean), by = Genes, .SDcols =! 'Genes']
    }
    
    exprs_mtx <- as.matrix(exprs_data.table[, -1])
    rownames(exprs_mtx) <- exprs_data.table$Genes
    
    return(exprs_mtx)
  }
  
  #### For DNA microarray data sets ####
  if(!is.null(DNA) && DNA == TRUE) {
    
    #### Checking if dna_files is a list ####
    if(!is.list(dna_studies)) stop("dna_studies needs to be a list")
    
    #### Checking if all elements in the list are character strings ####
    if(!all(sapply(dna_studies, is.character))) stop("Elements in dna_studies needs to be a list of character strings\nExample: list(dataset1 = 'GSEXXXX', dataset2 = 'GSEXXXX') ")
    
    ### Extracting list of DNA microarray expression matrices ####
    list_of_geo_data <- lapply(dna_studies, function( file ) {
      tryCatch({
        geo_data <- getGEO(GEO = file, AnnotGPL = TRUE)[[1]]
        if(is.null(geo_data) || length(geo_data) == 0) stop(paste("Failed to retrieve GEO dataset:", file))
        return(geo_data)
      }, error = function(e) {
        warning(paste("Error retrieving", file, ":", e$message))
        return(NULL) ### Return NULL in case of failure
      })
    })
    
    #### Extracting expression matrices ####
    list_of_dna_mtx <- lapply(list_of_geo_data, function( geo ) {
      if(!is.null(geo)) {
        return(exprs(geo))
      } else {
        return(NULL)
      }
    })
    
    #### Extracting feature data ####
    list_of_fData <- lapply(list_of_geo_data, function( geo ){
      if(!is.null(geo)) {
        return(fData(geo))
      } else {
        return(NULL)
      }
    })
    
    #### Annotating probe IDs ####
    for(i in 1:length(list_of_fData)) {
      fData <- list_of_fData[[i]]
      gene_cols <- c("Gene symbol", "ILMN_Gene", "GB_ACC")
      matched_col <- intersect(gene_cols, colnames(fData))
      
      if(length(matched_col) > 0) {
        selected_col <- matched_col[1]  ### Use the first matched column
        
        #### Ensure probe IDs are character strings ####
        probe_ids <- as.character(fData[[selected_col]])
        
        #### Check if conversion is needed ####
        if(!(selected_col %in% c("Gene symbol", "ILMN_Gene"))) {
          keytype <- ifelse(selected_col == "GB_ACC", "ACCNUM")
          
          #### Convert probe IDs to gene symbols ####
          fData$Genes <- mapIds(
            x = org.Hs.eg.db,
            keys = probe_ids,
            column = "SYMBOL",
            keytype = keytype,
            multiVals = "first"
          )
        } else {
          #### If already gene symbols, just assign them ####
          fData$Genes <- fData[[selected_col]]
        }
      }
      
      list_of_fData[[i]] <- fData  ### Save updated feature data back to list
    }
    
    #### Converting expression matrices to data.table ####
    list_of_dna_data_table <- list()
    
    for (i in 1:length(list_of_dna_mtx)) {
      mtx <- list_of_dna_mtx[[i]]
      fData <- list_of_fData[[i]]
      
      genes_vector <- as.character(fData$Genes)
      
      # Create data.table for each entry
      df <- data.frame(Genes = genes_vector, mtx)
      
      data_table <- as.data.table(df)
      
      # Append the data_table to the list
      list_of_dna_data_table[[i]] <- data_table
    }
    
    list_of_dna_mtx <- lapply(list_of_dna_data_table, function ( tbl ) {
      mtx <- handling_duplicates( tbl )
      return(mtx)
    })
    
    ##### Quantile normalization #####
    list_of_dna_mtx <- lapply(list_of_dna_mtx, function( mtx ) {
      gene_names <- rownames( mtx )
      
      exprs_mtx <- apply(mtx, 2, as.numeric)
      
      log_mtx <- log(exprs_mtx + 1)
      
      if(any(is.nan(log_mtx)) == TRUE) {
        log_mtx[is.nan(log_mtx)] <- 0
      }
      norm_mtx <- normalizeBetweenArrays(log_mtx, method = "quantile")
      
      #### Re-assigning gene names to the normalized matrix ####
      rownames(norm_mtx) <- gene_names
      return(norm_mtx)
    })
    
    #### Saving boxplots of normalized data ####
    for(i in 1:length(list_of_dna_mtx)) {
      pdf(paste0("boxplot_", names(dna_studies)[i], ".pdf"), width = 8, height = 8)
      boxplot <- boxplot(list_of_dna_mtx[[i]], outline = FALSE,
                         col = "skyblue", main = paste0("Log transformation/Quantile", "\n", names(dna_studies)[i]))
      dev.off()
    }
    
    return(list_of_dna_mtx)
  
  }
  
  #### For RNA-Seq data sets ####
  if(!is.null(RNA) && RNA == TRUE) {
    #### Checking if rna_mtx is a list ####
    if(!is.list(rna_mtx)) stop("rna_mtx needs to be a list")
    
    #### Checking if pData_RNA has "condition" column ####
    if(!is.list(pData_RNA)) stop("pData_RNA needs to be a list")
    
    #### Checking there is a column named "condition"
    if(!all(sapply(pData_RNA, function(x) "condition" %in% colnames(x)))) {
      stop("Each element in pData_RNA must have a 'condition' column")
    }
    
    #### Creating empty list ####
    list_of_norm_RNA_mtx <- list()
    list_of_dds <- list()
    
    for(i in 1:length(rna_mtx)) {
      mtx <- rna_mtx[[i]]
      pData <- pData_RNA[[i]]
      
      common_samples <- intersect(colnames(mtx), rownames(pData))
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      pData_fil <- pData[common_samples, , drop = FALSE]
      
      dds_1 <- DESeqDataSetFromMatrix(round(mtx_fil), colData =  pData_fil, design = ~ condition)
      dds_2 <- DESeq( dds_1 )
      
      norm_counts <- counts(dds_2, normalized = TRUE)
      
      list_of_dds[[i]] <- dds_2
      list_of_norm_RNA_mtx[[i]] <- norm_counts
    }
    
    #### Saving boxplots of normalized data ####
    for(i in 1:length(list_of_norm_RNA_mtx)) {
      pdf(paste0("boxplot_", names(rna_mtx)[i], ".pdf"), width = 8, height = 8)
      boxplot <- boxplot(list_of_norm_RNA_mtx[[i]], outline = FALSE,
                         col = "skyblue", main = paste0("Median of Ratios", "\n", names(rna_mtx)[i]))
      dev.off()
    }
    
    return( list(dds_objs = list_of_dds,
                 norm_counts = list_of_norm_RNA_mtx) )
  }
}

####
####

#### Finding common genes across studies #####
find_common_genes <- function(DNA = FALSE, RNA = FALSE, DNA_mtx = list(), RNA_mtx = list()) {
  if(DNA && !RNA) {
    return(Reduce(intersect, lapply(DNA_mtx, rownames)))
  }
  if(RNA && !DNA) {
    if(length(RNA_mtx) == 2) stop("Use RNA_mtx$norm_counts")
    return(Reduce(intersect, lapply(RNA_mtx, rownames)))
  }
  if(DNA && RNA) {
    list_of_matrices <- c(DNA_mtx, RNA_mtx)
    return(Reduce(intersect, lapply(list_of_matrices, rownames)))
  }
  return(NULL)
}

####
####

#### Meta-analysis function ####
meta_results <- function(list_of_studies) {

list_of_es <- lapply(list_of_studies, function( study ) {
  effect.sizes(study)
  })
  
  summary <- combine.effect.sizes(list_of_es)
  
  summary$g <- summary$g[!apply(is.na(summary$g) | is.infinite(summary$g), 1, any), ]
  summary$se.g <- summary$se.g[!apply(is.na(summary$se.g) | is.infinite(summary$se.g), 1, any), ]
  g_genes <- rownames(summary$g)
  
  summary$pooled.estimates$genes <- rownames(summary$pooled.estimates)
  summary$pooled.estimates <- summary$pooled.estimates %>%
    dplyr::filter(genes %in% g_genes)
  
  summary$pooled.estimates <- summary$pooled.estimates %>%
    dplyr::filter(!apply(is.na(.) | is.infinite(as.matrix(.)), 1, any))
  
  summary$pooled.estimates <- as.data.frame(lapply(summary$pooled.estimates, function(x) ifelse(x == 0, 1e-200, x)))
  rownames(summary$pooled.estimates) <- g_genes
  
  g <- summary$g
  se.g <- summary$se.g
  
  pool    <- summary$pooled.estimates[, "summary" ]
  names(pool) <- rownames(g)
  
  se.pool <- summary$pooled.estimates[, "se.summary" ]
  names(se.pool) <- rownames(g)
  
  x.label <- "Standardized Mean Difference (log2 scale)"
  
  consistent_genes <- c()
  
  for (gene in rownames(g)) {
    g_gene <- g[gene, ]
    if ((all(g_gene > 0) & pool[gene] > 0.5) | (all(g_gene < 0) & pool[gene] < -0.5)) {
      consistent_genes <- c(consistent_genes, gene)
    }
  }
  
  for (gene in consistent_genes) {
    g_gene <- g[gene, ]
    se.g_gene <- se.g[gene, ]
    study_names <- gsub("_g", "", names(g_gene))
    
    pdf(file = paste0(gene, ".pdf"), height = 3.5, width = 3.5)
    par(cex = 0.65)
    metaplot(g_gene, se.g_gene,
             labels = study_names,
             summn = pool[gene],
             sumse = se.pool[gene],
             sumnn = 1/se.pool[gene]^2,
             summlabel = "Summary effect",
             xlab = x.label,
             ylab = "",
             main = bquote(italic(.(gene))),
             colors = meta.colors(box = "violetred", lines = "plum", summary = "mediumpurple", 
                                  text = "black", axes = "black", zero = "black"),  
             boxsize = 1,
             lty.random = 1,  
             lwd.random = 2,   
             zero = 0,  
             col.zero = "black", 
             lty.zero = 3)
    dev.off()
  }
  
  return(list(summary = summary,
              genes = consistent_genes))
}

####
####

#### Function to store meta-analysis results ####
generate_list_for_meta_analysis <- function(DNA = FALSE, RNA = FALSE, dna_mtx = list(), rna_mtx = list(), pData = list(), study = list()) {
  
  #### Creating empty list to store meta-analysis results ####
  list_of_studies <- list()
  
  #### Preparing data for meta-analysis ####
  if(DNA && !RNA) {
    common_genes <- find_common_genes(DNA = TRUE, DNA_mtx = dna_mtx)
    for(i in 1:length(study)) {
      
      name <- names(study)[i]
      mtx <- dna_mtx[[i]]
      p <- pData[[i]]
      
      common_samples <- intersect(colnames(mtx), rownames(p))
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      p_fil <- p[common_samples, , drop = FALSE]
      
      #### Dynamically creating list ####
      list_of_studies[[name]] <- list(
        expr = mtx_fil[common_genes, ],
        pheno = p_fil$condition,
        keys = common_genes,
        class = as.numeric(ifelse(p_fil$condition == levels(p_fil$condition)[1], 0, 1))
      )
    }
  }
  
  #### Preparing data for meta-analysis ####
  if(RNA && !DNA) {
    common_genes <- find_common_genes(RNA = TRUE, RNA_mtx = RNA_matrices$norm_counts)
    for(i in 1:length(study)) {
      name <- names(study)[i]
      mtx <- rna_mtx[[i]]
      p <- pData[[i]]
      
      common_samples <- intersect(colnames(mtx), rownames(p))
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      p_fil <- p[common_samples, , drop = FALSE]
     
      #### Dynamically creating list ####
      list_of_studies[[name]] <- list(
        expr = mtx_fil[common_genes, ],
        pheno = p_fil$condition,
        keys = common_genes,
        class = as.numeric(ifelse(p_fil$condition == levels(p_fil$condition)[1], 0, 1))
      )
    }
  }
  
  
  if(DNA && RNA) {
    common_genes <- find_common_genes(DNA = TRUE, DNA_mtx = DNA_matrices, RNA = TRUE, RNA_mtx = RNA_matrices$norm_counts)
    for(i in 1:length(study)) {
      
      name <- names(study)[i]
      list_of_mtx <- c(DNA_mtx, RNA_mtx)
      mtx <- list_of_mtx[[i]]
      p <- pData[[i]]
      
      common_samples <- intersect(colnames(mtx), rownames(p))
      mtx_fil <- mtx[, common_samples, drop = FALSE]
      p_fil <- p[common_samples, , drop = FALSE]
      
      #### Dynamically creating list ####
      list_of_studies[[name]] <- list(
        expr = mtx_fil[common_genes, ],
        pheno = p_fil$condition,
        keys = common_genes,
        class = as.numeric(ifelse(p_fil$condition == levels(p_fil$condition)[1], 0, 1))
      )
    }
  }
  
  return(list(studies = list_of_studies,
              results = meta_results(list_of_studies)))
}

####
####

#### Return pvalues and qvalues ####
get_p_q_values <- function(meta_results) {
  
  cleaned_studies <- list(
    studies = lapply(meta_results$studies, function(study) {
      
      valid_rows <- complete.cases(study$class, study$pheno)
      expr_clean <- study$expr[, valid_rows]
      
      valid_genes <- apply(expr_clean, 1, function(gene_values) {
        all(!is.na(gene_values) & !is.nan(gene_values) & !is.infinite(gene_values))
      })
      expr_clean <- expr_clean[valid_genes, ]
      
      study_cleaned <- list(
        expr = expr_clean,
        pheno = study$pheno[valid_rows],
        keys = study$keys[valid_genes],
        class = study$class[valid_rows]
      )
      return(study_cleaned)
    }),
    
    results = meta_results$results
  )
  
  list.of.sigs <- lapply(cleaned_studies$studies, function( study ) {
    fisher_pval <- ttest.Pvalues(study)
    genes <- cleaned_studies$results$genes
    fisher_pval <- fisher_pval %>% dplyr::filter(keys%in%genes)
  })
  
  output.Fisher <- as.data.frame(sum.of.logs(list.of.sigs)) # For fisher p-value
  adj_pval <- adjust.fisher(output.Fisher=output.Fisher)[, c("F.Qval.up", "F.Qval.down")]
  
  genes <- cleaned_studies$results$genes
  pool <- cleaned_studies$results$summary$pooled.estimates[genes, "summary"]
  se.pool <- cleaned_studies$results$summary$pooled.estimates[genes, "se.summary" ]
  
  es_pval <- numeric()
  es_qval <- numeric()
  
  for(i in 1:nrow(output.Fisher)) {
    if(output.Fisher$F.pval.up[i] <= 0.05) {
      es_pval <- c(es_pval, output.Fisher$F.pval.up[i])
    }
    if(output.Fisher$F.pval.down[i] <= 0.05) {
      es_pval <- c(es_pval, output.Fisher$F.pval.down[i])
    }
  }
  
  for(i in 1:nrow(adj_pval)) {
    if(adj_pval$F.Qval.up[i] <= 0.05) {
      es_qval <- c(es_qval, adj_pval$F.Qval.up[i])
    }
    if(adj_pval$F.Qval.down[i] <= 0.05) {
      es_qval <- c(es_qval, adj_pval$F.Qval.down[i])
    }
  }
  
  sig <- data.frame(Genes = genes,
                    Pooled_es = pool,
                    Pooled_es_se = se.pool,
                    es_pval = es_pval,
                    es_qval = es_qval)
  return(sig)
}
