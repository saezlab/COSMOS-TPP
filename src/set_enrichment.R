#' Set enrichment for eprot or pprot containers
#' @param obj object of class eprot or pprot
#' @param ... additional parameters
#' @export
set_enrichment <- function(obj, ...) {
  UseMethod("set_enrichment", obj)
}

#' Set enrichment for eprot containers
#'
#' This methods performs an enrichment test of differential testing results from the
#' eprot container against sets of genes or proteins.
#'
#' @param obj object of class eprot
#' @param ontologies ontologies either a character vector or a list of lists. If given
#'                   as character vector must contain valid \code{ontology} arguments for
#'                   \code{\link{get_ontology}}.
#' @param gene_mapping two-column data frame or matrix with mappings between Uniprot
#'                     protein accession and gene IDs. This is needed if the Ontologies
#'                     are given as sets of gene IDs. If no mapping is provided
#'                     the ontologies are assumed to be provided on protein
#'                     level. Default to NULL.
#' @param alpha the alpha used for calling significant proteins (using the adjusted p), default to 0.05.
#' @param log2fc_cutoff the log2 fold-change cutoff for calling significant proteins, default to log2(1.5).
#' @param alternative the test alternative for the set test, either "greater" (default) or "less".
#' @param min_n_sig_set the minimal number of significant proteins/genes in a set in order for the set
#'                      to be tested.
#' @param padj_method the method for p-value adjustment for the sets of one ontology. Default to "BH".
#' @param ... further arguments to \code{\link{enrich_test}}
#'
#' @return an object of class eprot with the enrichment results assigned to
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom tibble tibble
#' @importFrom dplyr filter pull left_join group_by ungroup summarize select mutate
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ## get mapping from Uniprot accessions to EntrezIDs
#' ## this requires a connection (con) to the Cellzome Oracle database
#' gene_mapping <- czutils::map_uniprot_ids(ids = get_features(ep), con = con, type = "EntrezID")
#'
#' ## perform a set enrichment against Metabase maps (pathways)
#' ## this requires the metabaser package to be installed and an active connection
#' ## to the Metabase server
#' ep <- set_enrichment(ep, ontologies = "mb_maps", gene_mapping = gene_mapping)
#' }
#'
#' @rdname set_enrichment
#' @author Stephan Gade
set_enrichment.eprot <- function(obj,
                                 ontologies,
                                 gene_mapping = NULL,
                                 alpha = 0.05,
                                 log2fc_cutoff = log2(1.5),
                                 alternative = c("greater", "less"),
                                 min_n_sig_set = 2,
                                 padj_method = "BH",
                                 ...) {
  if (is.character(ontologies)) {
    onts <- lapply(ontologies, get_ontology)
    names(onts) <- ontologies
    ontologies <- onts
    rm(onts)
  } else if (is.list(ontologies)) {
    if (is.null(names(ontologies)) || !all(sapply(ontologies, is.list))) {
      stop("Ontologies have to be provided as either character vector or named list of lists!")
    }
  } else {
    stop("Ontologies have to be provided as either character vector or named list of lists!")
  }

  if (!is.null(gene_mapping)) {
    if (!is.data.frame(gene_mapping) && !is.matrix(gene_mapping) && ncol(gene_mapping) != 2) {
      stop("Gene mapping has to be given as a data frame or a matrix with 2 columns.")
    }

    gene_mapping <- as.data.frame(gene_mapping)
    colnames(gene_mapping) <- c("protein", "gene")

    # filter NA entries in gene_mapping
    gene_mapping <- gene_mapping %>%
      filter(!is.na(protein) & !is.na(gene))

    if (nrow(gene_mapping) == 0) {
      stop("The protein->gene mapping does not contain non NA entries!")
    }
  }

  if (!check_slots(obj, "diff_results")) {
    stop("The eprot container does not contain a slot with differential testing results. No set enrichment can be performed!")
  }

  if (check_slots(obj, "enrich_results")) {
    enrich_results <- obj[["enrich_results"]]

    if (!any(names(enrich_results) == "ORA")) {
      enrich_results$ORA <- list()
    }
  } else {
    enrich_results <- list(ORA = list())
  }

  m <- names(ontologies) %in% names(enrich_results$ORA)
  if (any(m)) {
    warning(
      "The eprot container contains already enrichment test results with the ontologies ",
      paste(names(ontologies)[m], collapse = ", "),
      ". The results will be overwritten."
    )
  }

  diff_results <- obj[["diff_results"]]

  for (ont in names(ontologies)) {

    eres <- list()
    for (cont in unique(diff_results$contrast)) {

      data_test <- diff_results %>%
        filter(contrast == cont)

      # filter out proteins with missing adj. p values or log2 fold-changes
      # these cannot be significant and shouldn't be in the universe
      data_test <- data_test %>%
        filter(!is.na(logFC) & !is.na(adj.P.Val))

      universe <- unique(data_test$ID)

      genelist <- data_test %>%
        filter(adj.P.Val < alpha & abs(logFC) > log2fc_cutoff) %>%
        pull(ID) %>%
        unique()

      if (length(genelist) > 0) {

        if (!is.null(gene_mapping)) {

          if (!any(universe %in% gene_mapping[, 1])) {
            stop("None of the proteins for contrast ", cont, " could be found in the gene mapping table. Something is wrong!")
          }

          if (!any(genelist %in% gene_mapping[, 1])) {
            stop("None of the sign. proteins for contrast ", cont, " could be found in the gene mapping table. Something is wrong!")
          }

          universe <- unique(gene_mapping[gene_mapping[, 1] %in% universe, 2])
          genelist <- unique(gene_mapping[gene_mapping[, 1] %in% genelist, 2])

        }

        eres[[cont]] <- enrich_test(
          genelist = genelist,
          universe = universe,
          sets = ontologies[[ont]],
          alternative = alternative,
          min_n_sig_set = min_n_sig_set,
          padj_method = padj_method,
          ...
        )

        if (!is.null(eres[[cont]])) {

          # replace list entries with sign. proteins/genes in results
          # with collapsed list of gene symbols and protein/gene IDs
          sig_sets <- tibble(
            set_name = rep(eres[[cont]]$set_name,
                           sapply(eres[[cont]]$sig_set, length)),
            id = unlist(eres[[cont]]$sig_set)
          )

          if (!is.null(gene_mapping)) {

            sig_sets_proteins <- sig_sets %>%
              left_join(gene_mapping, by = c("id" = "gene")) %>%
              left_join(distinct(data_test[, c("ID", "gene_name")]),
                        by = c("protein" = "ID")) %>%
              mutate(sig_set_proteins = paste0(gene_name, " (", id,
                                               "|", protein, ")")) %>%
              group_by(set_name) %>%
              summarize(sig_set_proteins = paste(sig_set_proteins,
                                                 collapse = ", ")) %>%
              ungroup()

          } else {

            sig_sets_proteins <- sig_sets %>%
              left_join(distinct(data_test[, c("ID", "gene_name")]),
                        by = c("id" = "ID")) %>%
              mutate(sig_set_proteins = paste0(gene_name, " (", id, ")")) %>%
              group_by(set_name) %>%
              summarize(sig_set_proteins = paste(sig_set_proteins,
                                                 collapse = ", ")) %>%
              ungroup()

          }

          eres[[cont]] <- eres[[cont]] %>%
            left_join(sig_sets_proteins, by = "set_name") %>%
            dplyr::select(-sig_set)

        }

      } else {
        warning("No significant proteins for given criteria in contrast ", cont, ".")
      }
    }

    # if all contrasts did not return a valid enrichment result
    # we store an empty table, so we know the enrichment test
    # was performed and no result could be found
    m <- sapply(eres, is.null)
    if (all(m)) {

      eres <- tibble(
        set_name = character(0),
        n_sig = numeric(0),
        n_set = numeric(0),
        n_sig_set = numeric(0),
        n_universe = numeric(0),
        P.Value = numeric(0),
        adj.P.Val = numeric(0),
        sig_set_proteins = character(0),
        contrast = character(0)
      )

    } else {

      eres <- eres[!m]
      eres <- do.call("rbind", eres) %>%
        mutate(contrast = rep(names(eres), sapply(eres, nrow)))

    }

    attr(eres, "alpha") <- alpha
    attr(eres, "log2fc_cutoff") <- log2fc_cutoff
    attr(eres, "padj_method") <- padj_method
    attr(eres, "alternative") <- alternative

    enrich_results$ORA[[ont]] <- eres
  }

  obj[["enrich_results"]] <- enrich_results

  obj
}



#' Set enrichment for pprot containers
#'
#' This methods performs an enrichment test of differential testing results from the
#' pprot container against sets of genes or proteins.
#'
#' @param obj object of class pprot
#' @param ontologies ontologies either a character vector or a list of lists. If given
#'                   as character vector must contain valid \code{ontology} arguments for
#'                   \code{\link{get_ontology}}.
#' @param gene_mapping two-column data frame or matrix with mappings between Uniprot
#'                     protein accession and gene IDs. This is needed if the Ontologies
#'                     are given as sets of gene IDs. If no mapping is provided
#'                     the ontologies are assumed to be provided on protein
#'                     level. Default to NULL.
#' @param alpha the alpha used for calling significant proteins (using the adjusted p), default to 0.05.
#' @param log2fc_cutoff the log2 fold-change cutoff for calling significant proteins, default to log2(1.5).
#' @param alternative the test alternative for the set test, either "greater" (default) or "less".
#' @param min_n_sig_set the minimal number of significant proteins/genes in a set in order for the set
#'                      to be tested.
#' @param padj_method the method for p-value adjustment for the sets of one ontology. Default to "BH".
#' @param ... further arguments to \code{\link{enrich_test}}
#'
#' @return an object of class pprot with the enrichment results assigned to
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom tibble tibble
#' @importFrom dplyr filter pull left_join group_by ungroup summarize select mutate slice arrange
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ## get mapping from Uniprot accessions to EntrezIDs
#' ## this requires a connection (con) to the Cellzome Oracle database
#' gene_mapping <- czutils::map_uniprot_ids(ids = get_features(ep), con = con, type = "EntrezID")
#'
#' ## perform a set enrichment against Metabase maps (pathways)
#' ## this requires the metabaser package to be installed and an active connection
#' ## to the Metabase server
#' ep <- set_enrichment(ep, ontologies = "mb_maps", gene_mapping = gene_mapping)
#' }
#'
#' @rdname set_enrichment
#' @author Stephan Gade, adapted by Mathias Kalxdorf
set_enrichment.pprot <- function(obj,
                                 ontologies,
                                 gene_mapping = NULL,
                                 alpha = 0.05,
                                 log2fc_cutoff = log2(1.5),
                                 alternative = c("greater", "less"),
                                 min_n_sig_set = 2,
                                 padj_method = "BH",
                                 ...) {
  if (is.character(ontologies)) {
    onts <- lapply(ontologies, get_ontology)
    names(onts) <- ontologies
    ontologies <- onts
    rm(onts)
  } else if (is.list(ontologies)) {
    if (is.null(names(ontologies)) || !all(sapply(ontologies, is.list))) {
      stop("Ontologies have to be provided as either character vector or named list of lists!")
    }
  } else {
    stop("Ontologies have to be provided as either character vector or named list of lists!")
  }

  if (!is.null(gene_mapping)) {
    if (!is.data.frame(gene_mapping) && !is.matrix(gene_mapping) && ncol(gene_mapping) != 2) {
      stop("Gene mapping has to be given as a data frame or a matrix with 2 columns.")
    }

    gene_mapping <- as.data.frame(gene_mapping)
    colnames(gene_mapping) <- c("protein", "gene")

    # filter NA entries in gene_mapping
    gene_mapping <- gene_mapping %>%
      filter(!is.na(protein) & !is.na(gene))

    if (nrow(gene_mapping) == 0) {
      stop("The protein->gene mapping does not contain non NA entries!")
    }
  }

  if (!check_slots(obj, "diff_results")) {
    stop("The pprot container does not contain a slot with differential testing results. No set enrichment can be performed!")
  }

  if (check_slots(obj, "enrich_results")) {
    enrich_results <- obj[["enrich_results"]]

    if (!any(names(enrich_results) == "ORA")) {
      enrich_results$ORA <- list()
    }
  } else {
    enrich_results <- list(ORA = list())
  }

  m <- names(ontologies) %in% names(enrich_results$ORA)
  if (any(m)) {
    warning(
      "The pprot container contains already enrichment test results with the ontologies ",
      paste(names(ontologies)[m], collapse = ", "),
      ". The results will be overwritten."
    )
  }

  diff_results <- obj[["diff_results"]]

  for (ont in names(ontologies)) {
    eres <- list()
    for (cont in unique(diff_results$contrast)) {
      data_test <- diff_results %>%
        filter(contrast == cont)

      # aggregate data to phospho-protein level
      # for each psite per protein, select entry with the smallest adj.P.Val
      # filter out sites without valid logFC or adj. p-value
      data_test <- data_test %>%
        filter(!is.na(logFC) & !is.na(adj.P.Val)) %>%
        group_by(representative, gene_name) %>%
        arrange(adj.P.Val) %>%
        slice(1) %>%
        ungroup() %>%
        select(
          logFC, P.Value, adj.P.Val, contrast, ID = representative, gene_name
        )

      # get most significant psite per protein
      # temp2 <- suppressWarnings(aggregate(data_test[, "adj.P.Val"],
      #                                     by = list(ID = data_test$representative, gene_name = data_test$gene_name), FUN = min, na.rm = T))
      # # get corresponding logFC and uncorrected pval
      # temp2$logFC <- data_test$logFC[match(temp2$adj.P.Val, data_test$adj.P.Val)]
      # temp2$P.Value <- data_test$P.Value[match(temp2$adj.P.Val, data_test$adj.P.Val)]

      # # store
      # data_test <- data.frame(
      #   logFC = temp2$logFC,
      #   P.Value = temp2$P.Value,
      #   adj.P.Val = temp2$adj.P.Val,
      #   contrast = cont,
      #   ID = temp2$ID,
      #   gene_name = temp2$gene_name
      # )
      # rownames(data_test) <- make.unique(data_test$gene_name)


      universe <- unique(data_test$ID)

      genelist <- data_test %>%
        filter(adj.P.Val < alpha & abs(logFC) > log2fc_cutoff) %>%
        pull(ID) %>%
        unique()

      if (length(genelist) > 0) {
        if (!is.null(gene_mapping)) {
          if (!any(universe %in% gene_mapping[, 1])) {
            stop("None of the proteins for contrast ", cont, " could be found in the gene mapping table. Something is wrong!")
          }

          if (!any(genelist %in% gene_mapping[, 1])) {
            stop("None of the sign. proteins for contrast ", cont, " could be found in the gene mapping table. Something is wrong!")
          }

          universe <- unique(gene_mapping[gene_mapping[, 1] %in% universe, 2])
          genelist <- unique(gene_mapping[gene_mapping[, 1] %in% genelist, 2])
        }

        eres[[cont]] <- enrich_test(
          genelist = genelist,
          universe = universe,
          sets = ontologies[[ont]],
          alternative = alternative,
          min_n_sig_set = min_n_sig_set,
          padj_method = padj_method,
          ...
        )

        if (!is.null(eres[[cont]])) {

          # replace list entries with sign. proteins/genes in results
          # with collapsed list of gene symbols and protein/gene IDs
          sig_sets <- tibble(
            set_name = rep(eres[[cont]]$set_name, sapply(eres[[cont]]$sig_set, length)),
            id = unlist(eres[[cont]]$sig_set)
          )

          if (!is.null(gene_mapping)) {
            sig_sets_proteins <- sig_sets %>%
              left_join(gene_mapping, by = c("id" = "gene")) %>%
              left_join(distinct(data_test[, c("ID", "gene_name")]),
                        by = c("protein" = "ID")) %>%
              mutate(sig_set_proteins = paste0(gene_name, " (", id,
                                               "|", protein, ")")) %>%
              group_by(set_name) %>%
              summarize(sig_set_proteins = paste(sig_set_proteins,
                                                 collapse = ", ")) %>%
              ungroup()
          } else {
            sig_sets_proteins <- sig_sets %>%
              left_join(distinct(data_test[, c("ID", "gene_name")]),
                        by = c("id" = "ID")) %>%
              mutate(sig_set_proteins = paste0(gene_name, " (", id, ")")) %>%
              group_by(set_name) %>%
              summarize(sig_set_proteins = paste(sig_set_proteins,
                                                 collapse = ", ")) %>%
              ungroup()
          }

          eres[[cont]] <- eres[[cont]] %>%
            left_join(sig_sets_proteins, by = "set_name") %>%
            dplyr::select(-sig_set)

        }

      } else {
        warning("No significant phospho-proteins for given criteria in contrast ", cont, ".")
      }
    }

    # if all contrasts did not return a valid enrichment result
    # we store an empty table, so we know the enrichment test
    # was performed and no result could be found
    m <- sapply(eres, is.null)
    if (all(m)) {

      eres <- tibble(
        set_name = character(0),
        n_sig = numeric(0),
        n_set = numeric(0),
        n_sig_set = numeric(0),
        n_universe = numeric(0),
        P.Value = numeric(0),
        adj.P.Val = numeric(0),
        sig_set_proteins = character(0),
        contrast = character(0)
      )

    } else {

      eres <- do.call("rbind", eres) %>%
        mutate(contrast = rep(names(eres), sapply(eres, nrow)))

    }

    attr(eres, "alpha") <- alpha
    attr(eres, "log2fc_cutoff") <- log2fc_cutoff
    attr(eres, "padj_method") <- padj_method
    attr(eres, "alternative") <- alternative

    enrich_results$ORA[[ont]] <- eres
  }

  obj[["enrich_results"]] <- enrich_results

  obj
}


#' Helper function to perform a set enrichment test.
#'
#' The function calculates an enrichment test for a list of genes against
#' some gene sets using a hyper-geometric test (Fisher's exaxt test).
#'
#' @param genelist vector of IDs (e.g. Entrez gene IDs)
#' @param sets list with gene sets (must be the same IDs as in genelist)
#' @param universe the universe (background list) of IDs. If missing,
#'                 the non-redudant union of all sets is used as universe.
#' @param restrict_universe boolean indicating whether the universe should
#'                          be restricted to only genes also present in the
#'                          the gene sets (TRUE) or used as provided (FALSE).
#'                          Default to FALSE.
#' @param alternative the test alternative for the set test, either
#'                    "greater" (default) or "less".
#' @param min_n_sig_set the minimal number of significant proteins/genes
#'                      in a set in order for the set to be tested. Please
#'                      note, that filtering is applied after the adjusted
#'                      p-values are calculated.
#' @param padj_method the method for p-value adjustment for the sets of
#'                    one ontology. Default to "BH".
#' @param verbose if TRUE print additional messages, default to FALSE
#'
#' @return a tibble with the test results:
#'         - set_name: the name of the set (e.g. name of the pathway),
#'         - n_sig number of IDs (genes/proteins) in the gene list
#'         - n_set number of IDs (gene/proteins) in the set
#'         - n_sig_set the number of IDs from the input list (e.g. sign. genes) in the set
#'         - n_universe number of IDs in the universe (background list)
#'         - p the hyper-geometrical p-value,
#'         - padj the adjusted p-value
#'         - sig_set a list column holding the overlap of IDs in the input list and the set
#'
#'         If the gene list and the given sets do not overlap or none of the
#'         sets contain at least \code{min_n_sig_set} genes from the gene
#'         list the function returns NULL.
#'
#' @importFrom magrittr "%>%"
#' @importFrom tibble tibble
#' @importFrom dplyr filter select mutate
#' @importFrom stats phyper
#' @export
#'
#' @note For efficiency reasons the function uses \code{phyper} instead of
#'       \code{fisher.test}.
#'
#' @rdname set_enrichment
#' @author Stephan Gade
enrich_test <- function(genelist, sets, universe = NULL,
                        restrict_universe = FALSE,
                        alternative = c("greater", "less"),
                        min_n_sig_set = 2,
                        padj_method = "BH",
                        verbose = F) {

  alternative <- match.arg(alternative)

  if (missing(sets)) {
    stop("No gene sets to test aginst specified!")
  }

  if (!is.list(sets)) {
    stop("Gene sets must be a list!")
  }

  if (!all(sapply(sets, is.vector))) {
    stop("Gene sest should be a list of vectors representing genes!")
  }

  if (missing(genelist)) {
    stop("No data is specified!")
  }

  if (!is.vector(genelist)) {
    stop("Data should be represented as vector of genes!")
  }

  if (length(genelist) > 2000) {
    message("Note that your input list contains over 2,000 genes. This is a very large signature.")
  }

  if (is.null(universe)) {
    universe <- unique(unlist(sets))
  }

  if (!any(universe %in% unlist(sets))) {
    stop("None of the genes in the universe is in the sets! Are these compatible?")
  }

  if (restrict_universe) {
    universe <- intersect(universe, unlist(sets))
  }

  if (!any(genelist %in% universe)) {
    stop("None of the genes in the gene list is in the universe! Are these compatible?")
  }

  genelist <- intersect(genelist, universe)

  if (verbose) {
    message("Processsing list of ", length(genelist), " genes")
  }


  if (!any(genelist %in% unlist(sets))) {
    warning("None of the genes in the list overlap with any set!")
    return(NULL)
  }

  g1 <- length(genelist)
  total <- length(universe)

  res <- lapply(
    sets,
    function(set) {

      g2 <- length(intersect(set, universe))
      overlap <- length(intersect(genelist, set))

      p <- switch(alternative,
                  "greater" = phyper(overlap - 1, g2, total - g2, g1,
                                     lower.tail = F),
                  "less" = phyper(overlap, g2, total - g2, g1,
                                  lower.tail = T))

      tibble(
        n_sig = g1, n_set = g2,
        n_sig_set = overlap,
        n_universe = total,
        P.Value = p,
        sig_set = list(intersect(genelist, set))
      )
    }
  )

  res <- do.call("rbind", res) %>%
    mutate(set_name = rep(names(res), sapply(res, nrow)))

  res <- res %>%
    mutate(adj.P.Val = p.adjust(P.Value, method = padj_method)) %>%
    dplyr::select(
      set_name, n_sig, n_set, n_sig_set, n_universe,
      P.Value, adj.P.Val, sig_set
    )

  if (!is.null(min_n_sig_set)) {
    res <- res %>%
      filter(n_sig_set >= min_n_sig_set)
  }

  if (nrow(res) == 0) {
    warning(
      "None of the sets contain at least ", min_n_sig_set,
      " genes from the imput gene list!"
    )
    return(NULL)
  }

  attr(res, "alternative") <- alternative

  res
}
