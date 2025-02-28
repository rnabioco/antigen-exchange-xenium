
# Utility ----

#' Save Seurat object and meta.data
#' 
#' @param ob_in Seurat object to save.
#' @param prfx Prefix to use for saved files. If set to NULL, the name of the
#' object is used.
#' @param ob_dir Directory to save files. The default is the current project
#' directory.
#' @export
save_objs <- function(ob_in, prfx = NULL, ob_dir = "") {
  
  if (is.null(prfx)) {
    prfx <- deparse(substitute(ob_in))
  }
  
  ob_in %>%
    qs::qsave(file.path(ob_dir, str_c(prfx, ".qs")))
  
  if ("Seurat" %in% class(ob_in)) {
    ob_in@meta.data %>%
      tibble::as_tibble(rownames = "cell_id") %>%
      vroom::vroom_write(file.path(ob_dir, str_c(prfx, ".tsv.gz")))
  }
}

#' Load slimmed down Seurat object
#' 
#' @param file File path to qs object
#' @param diet Return slimmed down object
load_obj <- function(file, diet = TRUE) {
  obj <- qread(file)
  
  if (diet) {
    obj <- obj %>%
      DietSeurat(
        assays = "RNA",
        counts = FALSE
      )
  }
  
  gc()
  
  obj
}

#' Export counts and meta.data tables
#' 
#' @param sobj_in Seurat object
#' @param assays Assays to include in counts matrix
#' @param feat_type Feature type for specified assay, provide a feature type
#' for each assay
#' @param gene_prefix Prefix to add to gene names for counts matrix
#' @param columns meta.data columns to export
#' @param out_dir Output directory
#' @param file_prefix Prefix to add to output files
#' @return Counts and meta.data tables
#' @export
export_matrices <- function(sobj_in, assays = "RNA", feat_type = "Gene Expression",
                            gene_prefix = "", columns, out_dir, file_prefix = "") {
  
  if (length(gene_prefix) == 1) {
    gene_prefix <- rep(gene_prefix, length(assays))
  }
  
  # Format count matrices
  assays <- purrr::set_names(feat_type, assays)
  
  mat <- assays %>%
    purrr::imap(~ Seurat::GetAssayData(sobj_in, assay = .y, "counts"))
  
  gene_names <- mat %>%
    purrr::map(rownames) %>%
    c(recursive = TRUE, use.names = FALSE)
  
  feat_type <- mat %>%
    purrr::imap(~ rep(assays[[.y]], nrow(.x))) %>%
    c(recursive = TRUE, use.names = FALSE)
  
  mat <- mat %>%
    purrr::reduce(rbind)
  
  barcodes <- colnames(mat)
  
  # Create HDF5 file
  counts_out <- file.path(out_dir, str_c(file_prefix, "count_matrix.h5"))
  
  h5file <- H5File$new(counts_out, mode = "w")
  
  h5file$create_group("matrix")
  h5file$create_group("matrix/features")
  
  h5file[["matrix/barcodes"]]              <- barcodes
  h5file[["matrix/features/name"]]         <- gene_names
  h5file[["matrix/features/feature_type"]] <- feat_type
  
  # Save matrix data (in sparse format)
  # * convert to COO format (i, j, x)
  triplet <- summary(mat)
  
  h5file[["matrix/data"]]    <- triplet$x                       # Non-zero values
  h5file[["matrix/indices"]] <- triplet$i - 1                   # Adjust indices for 0-based indexing
  h5file[["matrix/indptr"]]  <- c(0, cumsum(table(triplet$j)))  # Compressed column pointer
  h5file[["matrix/shape"]]   <- c(nrow(mat), ncol(mat))         # Matrix dimensions
  
  h5file$close()
  
  # Write meta.data table
  meta_out <- file.path(out_dir, str_c(file_prefix, "metadata.tsv.gz"))
  
  sobj_in@meta.data %>%
    tibble::as_tibble(rownames = "cell_id") %>%
    dplyr::select(any_of(columns)) %>%
    readr::write_tsv(meta_out)
}

#' Calculate a pseudo count for a given vector
#' 
#' @param x Vector of values to use for calculation
#' @param frac The pseudo count is calculated by multiplying the smallest
#' non-zero value by frac.
#' @return Values with pseudo count added
#' @export
add_pseudo <- function(x, frac = 0.5) {
  
  pseudo <- min(x[x > 0]) * frac
  
  x + pseudo
}

#' Run hypergeometric test
#' 
#' Arguments match those used for dhyper()
#' 
#' @param x number of white balls drawn
#' @param k number of total balls drawn
#' @param m number of white balls in urn
#' @param n number of black balls in urn
#' @param alt alternative hypothesis, 'greater' tests whether more white balls
#' were drawn than expected
#' @export
.calc_fisher <- function(x, k, m, n, alt = "two.sided") {
  tot <- m + n
  k   <- k - x
  m   <- m - x
  n   <- n - k
  
  # Example contingency table
  # the sum of the matrix should equal the total number of cells
  # 23  244  | 267
  # 51  3235 | 3286
  #
  # 74  3479 | 3553
  
  mat <- c(x, k, m, n) %>%
    matrix(nrow = 2)
  
  if (sum(mat) != tot) {
    stop(
      "To create contingency table, the following must be TRUE: ",
      "x + (k - x) + (m - x) + (n - k + x) == m + n"
    )
  }
  
  res <- mat %>%
    fisher.test(alternative = alt)
  
  res$p.value
}

#' Format p values for labels
#' 
#' modified from djvdj
.format_pvalue <- function(p, digits = 1, cutoffs = NULL, show_decimal = 0.1) {
  
  if (p == 0) {
    p <- str_c("italic(p) < 1*x*10^-16")
    
    return(p)
  }
  
  # Set p label based on vector of cutoffs
  if (!is.finite(p)) return(as.character(NA))
  
  if (!is.null(cutoffs)) {
    if (any(duplicated(cutoffs))) {
      cli::cli_abort("Cutoff values for p_label must be unique.")
    }
    
    # Set default labels when not provided by user
    if (is.null(names(cutoffs))) {
      cutoffs <- sort(cutoffs, decreasing = TRUE)
      
      names(cutoffs) <- purrr::imap_chr(
        cutoffs, ~ paste0(rep("*", .y), collapse = "")
      )
    }
    
    cutoffs <- sort(cutoffs)
    p_label <- as.character(NA)
    
    for (val in names(cutoffs)) {
      if (p < cutoffs[val]) {
        p_label <- val
        
        break()
      }
    }
    
    # Treat "value" as a keyword that will allow user to display actual
    # p-value for a certain cutoff
    # All custom labels need to be wrapped in quotes for parsing
    if (!identical(p_label, "value")) {
      if (!is.na(p_label)) p_label <- paste0("\'", p_label, "\'")
      
      return(p_label)
    }
  }
  
  # Format p-value label
  # label_scientific will round 0.095 to 0.1 when digits = 1
  if (round(p, digits + 1) >= show_decimal) return(as.character(round(p, 1)))
  
  p <- scales::label_scientific(digits = digits)(p)
  
  ex <- str_extract_all(p, "[+\\-][0-9]+$")
  
  p <- sub(paste0("\\", ex, "$"), "", p)
  
  ex <- as.numeric(ex)
  ex <- as.character(ex)
  
  p <- sub("e", "*x*10^", p)
  p <- paste0(p, ex)
  
  p
}

# Processing helpers ----

#' Classify cell types based on cluster mean expression
#' 
#' @param so_in Seurat object.
#' @param feats List of features to use for classifying clusters
#' @param filt Expression to use for filtering clusters, e.g. Cd3e < 0.1
#' @param type Cell type label to use for cells identified by filtering
#' expression
#' @param clst_col meta.data column containing cell clusters to use for
#' calculating mean expression.
#' @param type_col meta.data column to add cell type label.
#' @param summary_fn Function to use for summarizing marker gene expression.
#' @return Seurat object containing new cell type classifications.
#' @export
classify_markers <- function(so_in, feats, filt, type_label, clst_col, type_col,
                             summary_fn = mean) {
  
  clsts <- so_in
  
  if (is(so_in, "Seurat")) {
    clsts <- so_in %>%
      FetchData(unique(c(feats, clst_col, type_col)))
  }
  
  num_feats <- clsts %>%
    keep(is.numeric) %>%
    colnames()
  
  num_feats <- feats[feats %in% num_feats]
  chr_feats <- feats[!feats %in% num_feats]
  
  clsts <- clsts %>%
    group_by(!!sym(clst_col)) %>%
    summarize(
      across(all_of(num_feats), summary_fn),
      across(all_of(chr_feats), unique),
      .groups = "drop"
    )
  
  n_clsts <- nrow(clsts)
  
  clsts <- clsts %>%
    filter({{filt}}) %>%
    pull(clst_col) %>%
    as.character()
  
  if (n_distinct(so_in[[clst_col]]) != n_clsts) {
    warning("multiple values of one of the `feats` of type character are present for some clusters")
  }
  
  res <- so_in %>%
    mutate_meta(
      mutate,
      !!sym(type_col) := ifelse(
        !!sym(clst_col) %in% clsts,
        type_label,
        !!sym(type_col)
      )
    )
  
  res
}

#' Wrapper to create Seurat object
#' 
#' @param mat_dir Directory containing matrix generated by Cell Ranger.
#' @param proj_name Project name to include in meta.data table.
#' @param hash_ids Name of cell hashing antibodies included in matrix.
#' @param adt_count_min If CITE-seq was performed, this option will remove
#' antibodies where the sum total counts is less than adt_count_min.
#' @param gene_min Minimum number of detected genes for cell.
#' @param gene_max Maximum number of detected genes for cell.
#' @param mito_max Maximum percentage of mitochondrial reads for cell.
#' @param mt_str String to use for identifying mitochondrial genes.
#' @param rna_assay Name of RNA assay if multiple assays are being added to the
#' object (e.g. if CITE-seq data is included).
#' @param adt_assay Name of ADT assay for Seurat object.
#' @return Seurat object
#' @export
create_sobj <- function(mat_dir, proj_name = "SeuratProject", hash_ids = NULL, adt_count_min = 0,
                        gene_min = 250, gene_max = 5000, mito_max = 20, mt_str = "^mt-",
                        rna_assay = "Gene Expression", adt_assay = "Antibody Capture") {
  
  # Load matrices
  mat_list <- Seurat::Read10X(mat_dir)
  rna_mat  <- mat_list
  
  # Create Seurat object using gene expression data
  if (is_list(mat_list)) {
    rna_mat <- mat_list[[rna_assay]]
  }
  
  res <- rna_mat %>%
    Seurat::CreateSeuratObject(
      project   = proj_name,
      min.cells = 5
    )
  
  # Add antibody capture data to Seurat object
  if (is_list(mat_list)) {
    adt_mat <- mat_list[[adt_assay]]
    
    # Double check that cells match for both assays
    if (!identical(colnames(res), colnames(adt_mat))) {
      adt_mat <- adt_mat[, colnames(res)]
      
      warning("Not all cells are shared between RNA and ADT assays.")
    }
    
    # Remove ADT features that have low total counts and likely failed or
    # were omitted
    n_feats    <- nrow(adt_mat)
    count_sums <- rowSums(as.matrix(adt_mat))
    
    adt_mat <- adt_mat[count_sums >= adt_count_min, ]
    
    if (n_feats != nrow(adt_mat)) {
      warning("Some ADT features were removed due to low counts (<", adt_count_min, ").")
    }
    
    res[["ADT"]] <- Seurat::CreateAssayObject(adt_mat)
  }
  
  # Calculate percentage of mitochondrial reads
  res <- res %>%
    Seurat::PercentageFeatureSet(
      pattern  = mt_str, 
      col.name = "pct_mito"
    )
  
  # Add QC classifications to meta.data
  res <- res %>%
    mutate_meta(
      mutate,
      qc_class = case_when(
        pct_mito     > mito_max ~ "high_mito_reads",
        nFeature_RNA > gene_max ~ "high_gene_count",
        nFeature_RNA < gene_min ~ "low_gene_count",
        TRUE ~ "pass"
      )
    )
  
  res
}

#' Wrapper to normalize and scale Seurat object
#' 
#' @param sobj_in Seurat object.
#' @param rna_assay Name of RNA assay in object.
#' @param adt_assay Name of ADT assay in object.
#' @param cc_scoring Score cell cycle genes using cc.genes included in Seurat.
#' @param regress_vars Variables to regress out when scaling data.
#' @param rna_method Method to use with NormalizeData for RNA assay.
#' @param adt_method Method to use with NormalizeData for ADT assay.
#' @param scale_data Scale data after normalization.
#' @return Seurat object
#' @export
norm_sobj <- function(sobj_in, rna_assay = "RNA", adt_assay = "ADT", cc_scoring = FALSE,
                      regress_vars = NULL, rna_method = "LogNormalize", adt_method = "CLR",
                      scale_data = TRUE) {
  
  # Normalize counts
  res <- sobj_in %>%
    Seurat::NormalizeData(
      assay                = rna_assay,
      normalization.method = rna_method
    )
  
  # Score cell cycle genes
  if (cc_scoring) {
    s.genes <- cc.genes$s.genes %>%
      str_to_title()
    
    g2m.genes <- cc.genes$g2m.genes %>%
      str_to_title()
    
    res <- res %>%
      Seurat::CellCycleScoring(
        s.features   = s.genes,
        g2m.features = g2m.genes
      )
  }
  
  # Scale data
  # By default variable features will be used
  if (scale_data) {
    res <- res %>%
      Seurat::FindVariableFeatures(
        selection.method = "vst",
        nfeatures        = 2000
      ) %>%
      Seurat::ScaleData(vars.to.regress = regress_vars)
  }
  
  # Normalize ADT data
  if (!is.null(adt_assay) && adt_assay %in% names(res)) {
    res <- res %>%
      Seurat::NormalizeData(
        assay                = adt_assay,
        normalization.method = adt_method
      )
    
    if (scale_data) {
      res <- res %>%
        Seurat::ScaleData(assay = adt_assay)
    }
  }
  
  res
}

#' Perform k-means clustering on meta.data variable
#' 
#' @param dat data.frame with single column containing data to use for clustering.
#' @param k Number of clusters.
#' @param dat_clmn Column containing data to cluster
#' @param out_clmn Name of output column containing cell classifications.
#' @param clst_nms Labels to use for cell clusters.
#' @return data.frame containing cell clusters
#' @export
.run_km <- function(dat, k = 2, dat_clmn = "data", out_clmn = "km_cluster",
                    clst_nms = NULL) {
  
  # Data column name
  if (!is.null(colnames(dat))) {
    dat_clmn <- colnames(dat)
  }
  
  # K-means clustering
  res <- dat %>%
    stats::kmeans(centers = k)
  
  # Format results data.frame
  res <- res$cluster %>%
    data.frame()
  
  colnames(res) <- out_clmn
  
  if (!identical(rownames(dat), rownames(res))) {
    stop("Input and results rownames do not match.")
  }
  
  res <- dplyr::bind_cols(res, dat)
  
  # Add cluster names
  if (!is.null(clst_nms)) {
    if (length(clst_nms) != k) {
      stop("Must provide same number of cluster names as k.")
    }
    
    nms <- res %>%
      dplyr::group_by(!!sym(out_clmn)) %>%
      dplyr::summarize(mn = mean(!!sym(dat_clmn)), .groups = "drop") %>%
      dplyr::arrange(mn) %>%
      dplyr::pull(out_clmn)
    
    clst_nms <- purrr::set_names(clst_nms, nms)
    
    res <- res %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(
        !!sym(out_clmn) := clst_nms[as.character(!!sym(out_clmn))]
      ) %>%
      tibble::column_to_rownames()
  }
  
  res
}

#' Cluster meta.data variable using gaussian mixture model
#' 
#' @param dat data.frame with single column containing data to use for clustering.
#' @param k Number of clusters.
#' @param out_clmn Name of output column containing cell classifications.
#' @param clst_nms Labels to use for cell clusters.
#' @param prob Probability cutoff to use for classifying cells.
#' @param quiet Suppress output messages.
#' @return data.frame containing cell clusters
#' @export
.run_gmm <- function(dat, k = 2, out_clmn = "gmm_cluster", clst_nms = c("low", "high"),
                     prob = 0.5, quiet = TRUE) {
  
  if (length(clst_nms) != k) {
    stop("Must provide same number of cluster names as k.")
  }
  
  # Data column name
  dat_clmn <- "data"
  
  if (!is.null(colnames(dat))) {
    dat_clmn <- colnames(dat)
  }
  
  # Fit GMM for ova signal
  quiet_EM <- quietly(~ mixtools::normalmixEM(., k = k))
  
  if (!quiet) {
    quiet_EM <- mixtools::normalmixEM
  }
  
  set.seed(42)
  
  mdl <- dat %>%
    dplyr::pull(dat_clmn) %>%
    quiet_EM()
  
  if (quiet) {
    mdl <- mdl$result
  }
  
  # New column names
  comp_nms <- colnames(mdl$posterior)
  
  if (mdl$mu[1] > mdl$mu[2]) {
    clst_nms <- rev(clst_nms)
  }
  
  post              <- as.data.frame(mdl$posterior)
  colnames(post)    <- clst_nms
  names(comp_nms)   <- clst_nms
  names(mdl$mu)     <- clst_nms
  names(mdl$sigma)  <- clst_nms
  names(mdl$lambda) <- clst_nms
  
  # Format results data.frame
  clmns <- c("mu", "sigma", "lambda")
  
  clmns <- purrr::set_names(
    stringr::str_c(out_clmn, "_", clmns),
    clmns
  )
  
  res <- dplyr::bind_cols(dat, post)
  
  res <- res %>%
    tibble::rownames_to_column() %>%
    dplyr::mutate(
      !!sym(out_clmn) := if_else(
        !!sym(clst_nms[2]) >= prob,
        clst_nms[2],
        clst_nms[1]
      ),
      !!sym(clmns[["mu"]])     := mdl$mu[!!sym(out_clmn)],
      !!sym(clmns[["sigma"]])  := mdl$sigma[!!sym(out_clmn)],
      !!sym(clmns[["lambda"]]) := mdl$lambda[!!sym(out_clmn)],
      .before                   = !!sym(dat_clmn)
    ) %>%
    dplyr::select(-all_of(clst_nms)) %>%
    tibble::column_to_rownames()
  
  # Check that results match input data
  if (!identical(rownames(dat), rownames(res))) {
    stop("Input and results rownames do not match.")
  }
  
  res
}

#' Cluster meta.data variable
#' 
#' @param sobj_in Seurat object.
#' @param data_column meta.data column containing data to use for clustering.
#' @param k Number of clusters.
#' @param grp_column meta.data column containing cell labels to use for
#' dividing data. Clusters will be identified independently for each group.
#' @param filt Cell group present in grp_column to use for filtering cells before clustering.
#' All other cells will be labeled "other".
#' @param data_slot Slot to pull data from.
#' @param clust_column Name of meta.data column to output cell classifications.
#' @param clust_names Labels to use for cell clusters.
#' @param return_sobj Return a Seurat object. If FALSE a data.frame is
#' returned.
#' @param method Method to use for clustering, can be either "km" or "gmm".
#' @return Seurat object with cell classifications added to meta.data.
#' @export
cluster_signal <- function(sobj_in, data_column, k = 2, grp_column = NULL,
                           filt = NULL, data_slot = "counts",
                           clust_column = "clust",
                           clust_names = c("low", "high"),
                           return_sobj = TRUE, method = "gmm") {

  # Select method
  .funs <- list(
    km  = .run_km,
    gmm = .run_gmm
  )
  
  if (!method %in% names(.funs)) {
    stop("Must select one of the following methods: ", str_c(names(.funs), collapse = ", "))
  }
  
  .fun <- .funs[[method]]
  
  # Filter Seurat object
  so_flt <- sobj_in
  
  if (!is.null(filt)) {
    so_flt <- so_flt %>%
      subset(!!sym(grp_column) == filt)
  }
  
  # Split meta.data by grp_column
  so_flt <- list(so_flt)
  
  if (!is.null(grp_column)) {
    so_flt <- so_flt[[1]] %>%
      Seurat::SplitObject(grp_column)
  }
  
  # Cluster signal
  res <- so_flt %>%
    imap_dfr(~ {
      .x <- .x %>%
        Seurat::FetchData(data_column, slot = data_slot) %>%
        
        .fun(
          k        = k,
          out_clmn = clust_column,
          clst_nms = clust_names
        ) %>%
        
        tibble::rownames_to_column()
      
      if (!is.null(grp_column)) {
        .x <- .x %>%
          dplyr::mutate(!!sym(grp_column) := .y, .before = !!sym(clust_column))
      }  
      
      .x %>%  
        tibble::column_to_rownames()
    })
  
  # Return data.frame
  if (!return_sobj) {
    return(res)
  }
  
  # Add clusters to meta.data for input object
  res <- res %>%
    dplyr::select(-all_of(c(data_column, grp_column)))
  
  res <- sobj_in %>%
    Seurat::AddMetaData(res)
  
  # Add "other" label for cells not included in the comparison
  res <- res %>%
    djvdj::mutate_meta(mutate, !!sym(clust_column) := replace_na(!!sym(clust_column), "other"))
  
  res
}
