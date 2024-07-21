## Affymetrix HGU133 Plus 2.0 control probe custom code

## Check poly-A probe expression ----
polyA_check <- function(mat, meta, groupby) {
  if (!(identical(colnames(mat), rownames(meta)))) {
    stop("Colnames of expression matrix must be identical to rownames of metadata.")
  }
  if (!groupby %in% names(meta)) {
    stop("The groupby variable must be a column in the metadata.")
  }
  if (length(unique(meta[[groupby]])) < 2) {
    stop("The groupby variable must contain at least two groups.")
  }
  index <- grep(pattern = "AFFX-r2-Bs-.*|thr|lys|phe|dap|Thr|Lys|Phe|Dap", x = rownames(mat))
  polyA <- as_tibble(mat[index,])
  polyA <- polyA %>% 
    pivot_longer(cols = names(polyA),
                 names_to = "sample",
                 values_to = "expr" )
  polyA$group <- factor(meta[polyA$sample, groupby])
  by <- unique(polyA$group)
  if (length(by) == 2) {
    pval <- t.test(x = polyA$expr[which(polyA$group == levels(polyA$group)[1])], y = polyA$expr[which(polyA$group == levels(polyA$group)[2])])$p.value
    p <- ggplot(data = polyA, aes(x = group, y = expr, fill = group)) +
      geom_boxplot() +
      theme_linedraw() +
      geom_label(aes(x = length(by)/2, y = 10, label = paste("p =", round(pval,3)), fill = NULL)) +
      ggtitle("Average polyA probe expression across groups.")
  } else {
    pval <- aov(expr ~ group, data = polyA)
    pval <- summary(pval)[[1]][["Pr(>F)"]][1]
    p <- ggplot(data = polyA, aes(x = group, y = expr, fill = group)) +
      geom_boxplot() +
      theme_linedraw() +
      geom_label(aes(x = length(by)/2, y = 10, label = paste("p =", round(pval,3)), fill = NULL)) +
      ggtitle("Average polyA probe expression across groups.")
  }
  polyA <- as_tibble(mat[index,])
  polyA$probe <- rownames(mat[index,])
  polyA$probe <- tolower(str_extract(polyA$probe, "thr|lys|phe|dap|Thr|Lys|Phe|Dap"))
  polyA$probe <- factor(polyA$probe, levels = c("lys", "phe", "thr", "dap"))
  q <- data.frame(Sample = names(polyA)[1:(length(names(polyA))-1)], Test = rep(0, times = length(names(polyA))-1))
  for (i in 1:(length(names(polyA))-1)) {
    r <- polyA %>% 
      group_by(probe) %>% 
      summarise(mean = mean(.data[[q$Sample[i]]]))
    q$Test <- (r$mean[1] < r$mean[2] & r$mean[2] < r$mean[3] & r$mean[3] < r$mean[4])
  }
  return(list(p,`Are probes ordered properly?` = q))
}

## Check hybridization expression----
hybrid_check <- function(mat, meta, groupby) {
  if (!(identical(colnames(mat), rownames(meta)))) {
    stop("Colnames of expression matrix must be identical to rownames of metadata.")
  }
  if (!groupby %in% names(meta)) {
    stop("The groupby variable must be a column in the metadata.")
  }
  if (length(unique(meta[[groupby]])) < 2) {
    stop("The groupby variable must contain at least two groups.")
  }
  index <- grep(pattern = "cre|Cre|BioC|BioD|BioB|bioC|bioD|bioB", x = rownames(mat))
  polyA <- as_tibble(mat[index,])
  polyA <- polyA %>% 
    pivot_longer(cols = names(polyA),
                 names_to = "sample",
                 values_to = "expr" )
  polyA$group <- factor(meta[polyA$sample, groupby])
  by <- unique(polyA$group)
  if (length(by) == 2) {
    pval <- t.test(x = polyA$expr[which(polyA$group == levels(polyA$group)[1])], y = polyA$expr[which(polyA$group == levels(polyA$group)[2])])$p.value
    p <- ggplot(data = polyA, aes(x = group, y = expr, fill = group)) +
      geom_boxplot() +
      theme_linedraw() +
      geom_label(aes(x = length(by)/2, y = 10, label = paste("p =", round(pval,3)), fill = NULL)) +
      ggtitle("Average hybridization control probe expression across groups.")
  } else {
    pval <- aov(expr ~ group, data = polyA)
    pval <- summary(pval)[[1]][["Pr(>F)"]][1]
    p <- ggplot(data = polyA, aes(x = group, y = expr, fill = group)) +
      geom_boxplot() +
      theme_linedraw() +
      geom_label(aes(x = length(by)/2, y = 10, label = paste("p =", round(pval,3)), fill = NULL)) +
      ggtitle("Average hybridization control probe expression across groups.")
  }
  polyA <- as_tibble(mat[index,])
  polyA$probe <- rownames(mat[index,])
  polyA$probe <- tolower(str_extract(polyA$probe, "cre|Cre|BioC|BioD|BioB|bioC|bioD|bioB"))
  polyA$probe <- factor(polyA$probe, levels = c("biob", "bioc", "biod", "cre"))
  q <- data.frame(Sample = names(polyA)[1:(length(names(polyA))-1)], Test = rep(0, times = length(names(polyA))-1))
  for (i in 1:(length(names(polyA))-1)) {
    r <- polyA %>% 
      group_by(probe) %>% 
      summarise(mean = mean(.data[[q$Sample[i]]]))
    q$Test <- (r$mean[1] < r$mean[2] & r$mean[2] < r$mean[3] & r$mean[3] < r$mean[4])
  }
  return(list(p,`Are probes ordered properly?` = q))
}

## Check positive control expression----  
positive_check <- function(mat, meta, groupby) {
  if (!(identical(colnames(mat), rownames(meta)))) {
    stop("Colnames of expression matrix must be identical to rownames of metadata.")
  }
  if (!groupby %in% names(meta)) {
    stop("The groupby variable must be a column in the metadata.")
  }
  if (length(unique(meta[[groupby]])) < 2) {
    stop("The groupby variable must contain at least two groups.")
  }
  index <- grep(pattern = "200801_x_at|200887_s_at|209969_s_at|212581_x_at|213453_x_at|213867_x_at|217398_x_at|224594_x_at|AFFX-HSAC07/X00351_3_at|AFFX-HSAC07/X00351_5_at|AFFX-HSAC07/X00351_M_at|AFFX-HUMGAPDH/M33197_3_at|AFFX-HUMGAPDH/M33197_5_at|AFFX-HUMGAPDH/M33197_M_at|AFFX-HUMISGF3A/M97935_3_at|AFFX-HUMISGF3A/M97935_5_at|AFFX-HUMISGF3A/M97935_MA_at|AFFX-HUMISGF3A/M97935_MB_at", 
                x = rownames(mat))
  polyA <- as_tibble(mat[index,])
  polyA <- polyA %>% 
    pivot_longer(cols = names(polyA),
                 names_to = "sample",
                 values_to = "expr" )
  polyA$group <- factor(meta[polyA$sample, groupby])
  by <- unique(polyA$group)
  if (length(by) == 2) {
    pval <- t.test(x = polyA$expr[which(polyA$group == levels(polyA$group)[1])], y = polyA$expr[which(polyA$group == levels(polyA$group)[2])])$p.value
    p <- ggplot(data = polyA, aes(x = group, y = expr, fill = group)) +
      geom_boxplot() +
      theme_linedraw() +
      geom_label(aes(x = length(by)/2, y = 10, label = paste("p =", round(pval,3)), fill = NULL)) +
      ggtitle("Average positive control probe expression across groups.")
  } else {
    pval <- aov(expr ~ group, data = polyA)
    pval <- summary(pval)[[1]][["Pr(>F)"]][1]
    p <- ggplot(data = polyA, aes(x = group, y = expr, fill = group)) +
      geom_boxplot() +
      theme_linedraw() +
      geom_label(aes(x = length(by)/2, y = 10, label = paste("p =", round(pval,3)), fill = NULL)) +
      ggtitle("Average positive control probe expression across groups.")
  }
  return(p)
}