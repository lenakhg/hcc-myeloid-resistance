#!/bin/bash
#SBATCH --job-name=dgidb_P18excl_65
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=00:20:00
#SBATCH --partition=batch
#SBATCH --output=/ibex/user/alghamlk/HCC_P18_excluded/dgidb_65_%J.out
#SBATCH --error=/ibex/user/alghamlk/HCC_P18_excluded/dgidb_65_%J.err

echo "Job started: $(date)"
cd /ibex/user/alghamlk/HCC_P18_excluded/
module load R/4.5.0/gnu-12.2.0

Rscript - <<'EOF'
.libPaths(c("/home/alghamlk/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({library(httr2); library(jsonlite); library(dplyr)})

cat("=== DGIdb GraphQL API -- unbiased query of all 65 P18-excluded significant genes ===\n\n")

sig <- read.csv("/ibex/user/alghamlk/HCC_P18_excluded/data/deseq2_PRvsPD_P18excl_SIG.csv")
gene_col <- if ("gene" %in% names(sig)) sig$gene else sig[[1]]
genes_of_interest <- as.character(gene_col)
cat("Significant genes queried (all):", length(genes_of_interest), "\n\n")

query <- '
query getInteractions($genes: [String!]) {
  genes(names: $genes) {
    nodes {
      name
      interactions {
        drug { name approved }
        interactionScore
        interactionTypes { type }
        sources { sourceDbName }
      }
    }
  }
}'

resp <- tryCatch({
  request("https://dgidb.org/api/graphql") %>%
    req_body_json(list(query = query, variables = list(genes = genes_of_interest))) %>%
    req_headers("Content-Type" = "application/json") %>%
    req_timeout(120) %>% req_perform()
}, error = function(e) { cat("REQUEST FAILED:", conditionMessage(e), "\n"); NULL })

if (!is.null(resp)) {
  cat("HTTP status:", resp_status(resp), "\n\n")
  result <- resp_body_json(resp)
  rows <- list()
  for (gn in result$data$genes$nodes) {
    if (length(gn$interactions) == 0) next
    for (it in gn$interactions) {
      rows[[length(rows)+1]] <- data.frame(
        gene=gn$name, drug=it$drug$name,
        approved=ifelse(is.null(it$drug$approved),NA,it$drug$approved),
        score=ifelse(is.null(it$interactionScore),NA,it$interactionScore),
        interaction_types=paste(sapply(it$interactionTypes,function(x)x$type),collapse="; "),
        sources=paste(sapply(it$sources,function(x)x$sourceDbName),collapse="; ")) }
  }
  df <- bind_rows(rows)

  # ALL druggable genes (>=1 interaction), sorted
  druggable <- df %>% group_by(gene) %>%
    summarise(n_drugs=n_distinct(drug),
              n_approved=n_distinct(drug[approved==TRUE], na.rm=TRUE),
              .groups="drop") %>% arrange(desc(n_drugs))
  cat("=== ALL druggable genes among the 65 (>=1 DGIdb interaction) ===\n")
  print(as.data.frame(druggable), row.names=FALSE)
  cat("\nTotal druggable genes:", nrow(druggable), "of 65\n")

  # Your 3 nominated axes -- confirm the specific drugs
  cat("\n=== Confirming nominated drugs ===\n")
  for (pat in c("FOSTAMATINIB","GILTERITINIB","NAMODENOSON")) {
    hit <- df %>% filter(grepl(pat, toupper(drug)))
    if (nrow(hit)>0) {
      cat(sprintf("\n%s targets in our gene set:\n", pat))
      print(as.data.frame(hit %>% select(gene, drug, interaction_types, approved) %>% distinct()), row.names=FALSE)
    } else cat(sprintf("\n%s: not returned by DGIdb for these genes\n", pat))
  }

  write.csv(df, "/ibex/user/alghamlk/HCC_P18_excluded/data/dgidb_druggability_P18excl_65.csv", row.names=FALSE)
  write.csv(druggable, "/ibex/user/alghamlk/HCC_P18_excluded/data/dgidb_druggable_genes_P18excl_65.csv", row.names=FALSE)
  cat("\nSaved dgidb_druggability_P18excl_65.csv and dgidb_druggable_genes_P18excl_65.csv\n")
  cat("DGIdb query date:", as.character(Sys.Date()), "\n")
} else cat("\nNo results -- run on login node (needs internet).\n")
cat("\nDone:", as.character(Sys.time()), "\n")
EOF
echo "Job finished: $(date)"
