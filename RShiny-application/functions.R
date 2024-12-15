# Libraries ----
library(mailR)
library(rJava)
library(httr)
library(jsonlite)
library(cowplot)
library(nanoparquet)

# The gene_graph function ----
gene_graph <- function(geneName, pData, expr, column = "subject", ylim = c(min(expr),max(expr))) {
  if (!is.numeric(expr[,length(colnames(expr))])) {
    expr <- expr[,-length(colnames(expr))]
  }
  if (length(geneName) == 1) {
    gene <- tibble(expr = unlist(as.vector(expr[paste(geneName), ])), time = factor(pData[, "time"]),
                   ID = factor(pData[, "subject"]), trait = factor(pData[, column]))
    p <- ggplot(gene, aes(x = time, y = expr, group = ID)) +
      geom_line(aes(color = trait)) +
      theme_linedraw() +
      ggtitle(geneName) +
      theme(title = element_text(size = 15, face = "bold"), axis.title = element_text(size = 10)) +
      xlab(label = "Time of sampling") +
      ylab(label = "log2(normalized expression)") +
      lims(y = ylim) +
      labs(color = column)
    return(p)
  } else {
    stop("You can only graph 1 gene at a time")
  }
}
## Probe edition (for microarrays)
gene_graph_p <- function(probeName, pData, expr, column = "subject", ylim = c(min(expr),max(expr))) {
  symbol <- rownames(expr) # Move symbols to a new vector
  rownames(expr) <- expr[,length(colnames(expr))] # Move probes to rows
  expr <- expr[,-length(colnames(expr))] # Remove probe column
  if (length(probeName) == 1) {
    gene <- tibble(expr = unlist(as.vector(expr[paste(probeName), ])), time = factor(pData[, "time"]),
                   ID = factor(pData[, "subject"]), trait = factor(pData[, column]))
    if (any(is.na(gene$expr))) {
      stop("Probe ID not found in this dataset.")
    }
    p <- ggplot(gene, aes(x = time, y = expr, group = ID)) +
      geom_line(aes(color = trait)) +
      theme_linedraw() +
      ggtitle(paste(symbol[which(rownames(expr) == probeName)], ", ", probeName, sep = "")) +
      theme(title = element_text(size = 15, face = "bold"), axis.title = element_text(size = 10)) +
      xlab(label = "Time of sampling") +
      ylab(label = "log2(normalized expression)") +
      lims(y = ylim) +
      labs(color = column)
    return(p)
  } else {
    stop("You can only graph 1 gene at a time")
  }
}