source("init.R")
library(gtools)
library(GGally)
library(viridis)

## ---- CONST --------
SPECIES <- "human" # or "mouse"

## ---- rnahybrid --------
if (SPECIES == "human") {
  hybrid <-
    vroom("./data-matrix/human_cdna_rnahybrid_matrix_format.tsv") # input
} else {
  hybrid <-
    vroom("./data-matrix/mouse_cdna_rnahybrid_matrix_format.tsv") # input
}

## ---- expression matrix preparing --------
if (SPECIES == "human") {
  srna.exp <- "./data-matrix/human_rnahybrid_input_sRNA.tsv" %>% vroom() # input
  gene.exp <- vroom("./data-matrix/human_gene_tpm_matrix_mean.csv") # input
} else if (SPECIES == "mouse") {
  srna.exp <- "./data-matrix/mouse_rnahybrid_input_sRNA.tsv" %>% vroom() # input
  gene.exp <- vroom("./data-matrix/mouse_gene_tpm_matrix_mean.csv") # input
}

## ---- pattern search --------
patternSearch <- function(
    df, stage, total_stage, foldchange = 5, base_expression = 10, depletion = T) {
  remain <- total_stage %>% lubridate::setdiff(stage)
  totalLen <- length(total_stage)
  stageLen <- length(stage)
  remainLen <- length(remain)

  if (depletion) {
    swap <- stage
    stage <- remain
    remain <- swap
    stageLen <- length(stage)
    remainLen <- length(remain)
  }

  for (i in 1:stageLen) {
    df <- df %>%
      filter(.data[[stage[i]]] >= base_expression)
  }

  # foldchange filtering
  for (i in 1:stageLen) {
    for (j in 1:remainLen) {
      df <- df %>%
        filter(
          (.data[[stage[i]]] / ((.data[[remain[j]]]) + 1)) > foldchange
        )
    }
  }
  return(df)
}

depletion_func <- function(
    gene.exp,
    num = 1,
    total_stage = c("PN5", "E2C", "L2C", "4C", "8C"),
    depletion = T,
    base_expression = 10,
    foldchange = 5) {
  if (num == 1) {
    total_stage %>%
      set_names(total_stage) %>%
      map(
        ~ {
          stage <- .x
          patternSearch(
            gene.exp,
            stage = stage,
            total_stage,
            depletion = depletion,
            base_expression = base_expression,
            foldchange = 5
          )
        }
      )
  } else if (num == 2) {
    input <- combinations(length(total_stage), num, total_stage)
    res <- input %>%
      as.data.frame() %>%
      mutate(V3 = paste(V1, V2, sep = "_"))
    res <- res %>%
      pmap(
        ~ {
          stage <- c(..1, ..2)
          name <- ..3
          patternSearch(
            gene.exp,
            stage = stage,
            total_stage,
            depletion = depletion
          )
        }
      ) %>%
      set_names(res$V3)
  } else if (num == 3) {
    input <- combinations(length(total_stage), num, total_stage)
    res <- input %>%
      as.data.frame() %>%
      mutate(V4 = paste(V1, V2, V3, sep = "_"))
    res %>%
      pmap(
        ~ {
          stage <- c(..1, ..2, ..3)
          name <- ..4
          patternSearch(
            gene.exp,
            stage = stage,
            total_stage,
            depletion = depletion
          )
        }
      ) %>%
      set_names(res$V4)
  } else if (num == 4) {
    input <- combinations(length(total_stage), num, total_stage)
    res <- input %>%
      as.data.frame() %>%
      mutate(V5 = paste(V1, V2, V3, V4, sep = "_"))
    res %>%
      pmap(
        ~ {
          stage <- c(..1, ..2, ..3, ..4)
          name <- ..5
          patternSearch(
            gene.exp,
            stage = stage,
            total_stage,
            depletion = depletion
          )
        }
      ) %>%
      set_names(res$V5)
  }
}


if (SPECIES == "human") {
  # get the expression pattern
  total_stage <- c("1C", "4C", "8C", "Morula")
  gene.1 <- depletion_func(gene.exp, num = 1, total_stage = total_stage)
  gene.2 <- depletion_func(gene.exp, num = 2, total_stage = total_stage)
  gene.3 <- depletion_func(gene.exp, num = 3, total_stage = total_stage)
  srna.1 <- depletion_func(srna.exp, num = 1, total_stage = total_stage, depletion = F)
  srna.2 <- depletion_func(srna.exp, num = 2, total_stage = total_stage, depletion = F)
  srna.3 <- depletion_func(srna.exp, num = 3, total_stage = total_stage, depletion = F)
} else if (SPECIES == "mouse") {
  # get the expression pattern
  gene.1 <- depletion_func(gene.exp, num = 1)
  gene.2 <- depletion_func(gene.exp, num = 2)
  gene.3 <- depletion_func(gene.exp, num = 3)
  gene.4 <- depletion_func(gene.exp, num = 4)
  srna.1 <- depletion_func(sncRNA, num = 1, depletion = F)
  srna.2 <- depletion_func(sncRNA, num = 2, depletion = F)
  srna.3 <- depletion_func(sncRNA, num = 3, depletion = F)
  srna.4 <- depletion_func(sncRNA, num = 4, depletion = F)
}

# hybrid filtering
get_pattern <- function(hybrid, q, t) {
  q <- q %>%
    pull(gene) %>%
    funique()
  t <- t %>%
    pull(gene) %>%
    funique()
  hybrid %>%
    filter(
      query %in% q,
      gene_name %in% t
    )
}
write_data <- function(hybrid, srna, gene) {
  names <- gene %>%
    names() %>%
    cusDiscard("gene")
  for (stage in names) {
    out <- get_pattern(hybrid, srna[[stage]], gene[[stage]])
    out %>% glimpse()
    out %>%
      write_tsv(paste0(
        "./results-figure-data/",
        SPECIES,
        "_pattern/",
        SPECIES,
        "_",
        stage,
        "_pattern.tsv"
      ))
  }
}

# write data
if (SPECIES == "human") {
  write_data(hybrid, srna.1, gene.1)
  write_data(hybrid, srna.2, gene.2)
  write_data(hybrid, srna.3, gene.3)
} else if (SPECIES == "mouse") {
  write_data(hybrid, srna.1, gene.1)
  write_data(hybrid, srna.2, gene.2)
  write_data(hybrid, srna.3, gene.3)
  write_data(hybrid, srna.4, gene.4)
}
