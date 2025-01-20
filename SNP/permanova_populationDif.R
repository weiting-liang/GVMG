#dingqiuxia
#!/usr/bin/env Rscript

library(vegan)  # Used for PERMANOVA test
library(Hmisc)  # Used for data imputation

# Read input
args <- commandArgs(trailingOnly = TRUE)
d1 <- as.matrix(read.table(args[1], header = TRUE, row.names = 1, check.names=FALSE,sep = "\t"))
phe <- read.table(args[2], header = TRUE, row.names=1,sep = "\t")

# Replace empty values with NA
phe[phe == ""] <- NA

# Fill NA in Reads_number with median
phe$Reads_number <- as.numeric(impute(phe$Reads_number, median))

#rownames(d1) <- trimws(rownames(d1))
intern <- intersect(rownames(phe), rownames(d1))
phe <- phe[intern, , drop = FALSE]
d1 <- as.matrix(d1)[intern, intern]

# Ensure correct variable types
phe$Source_group <- as.factor(phe$Source_group)  # Ensure Source_group is a factor
continuous_vars <- setdiff(colnames(phe), "Source_group")
phe[continuous_vars] <- lapply(phe[continuous_vars], as.numeric)  # Ensure continuous variables are numeric

# Initialize result list
result_list <- list()

# Iterate over all variables, each as the target variable
variables <- colnames(phe)
for (target_var in variables) {
  covariates <- setdiff(variables, target_var)
  if ("Age" %in% covariates) {
    covariates <- setdiff(covariates, "Age")  # Remove Age
  }
  # Construct formula, excluding factor variables with only one level
non_na_indices <- complete.cases(phe[, c(covariates, target_var)])
  filtered_phe <- phe[non_na_indices, , drop = FALSE]
  
valid_covariates <- covariates[sapply(covariates, function(var) {
    !is.factor(filtered_phe[[var]]) || length(unique(filtered_phe[[var]])) > 1
  })]

  
  # Filter non-NA samples
  filtered_phe <- filtered_phe[, c(valid_covariates, target_var)]
  
 d2<-d1[rownames(filtered_phe),rownames(filtered_phe)]
 
 formula <- as.formula(
    paste("d2 ~", paste(valid_covariates, collapse = " + "), "+", target_var)
  )
  
  # Check sample size
  sample_size <- nrow(filtered_phe)

# Special condition: If the target variable is Source_group, ensure sample size > 5 for at least two groups
  if (target_var == "Source_group") {
    group_counts <- table(phe[[target_var]])
    valid_groups <- sum(group_counts > 5)
    if (valid_groups < 2) {
      result_list[[target_var]] <- c("SampleSize" = sample_size, "Df" = NA,"SumOfSqs" =NA,"F" =NA,"R2" = NA, "Pr" = NA)
      next
    }
  }
  if (sample_size > 1) {
    tryCatch({
      adonis_result <- adonis2(
        d2 ~ ., data = filtered_phe, formula = formula, permutations = 999
      )
        result_list[[target_var]]<-c("SampleSize" = sample_size,"Df" = adonis_result$Df[which(rown
ames(adonis_result)==target_var)],"SumOfSqs" = adonis_result$SumOfSqs[which(rownames(adonis_result
)==target_var)],"F" = adonis_result$F[which(rownames(adonis_result)==target_var)],"R2" = adonis_re
sult$R2[which(rownames(adonis_result)==target_var)],"Pr" = adonis_result$Pr[which(rownames(adonis_
result)==target_var)])
    }, error = function(e) {
      warning(sprintf("Error for variable '%s': %s", target_var, e$message))
      result_list[[target_var]] <- c("SampleSize" = sample_size, "Df" = NA,"SumOfSqs" =NA,"F" =NA,
"R2" = NA, "Pr" = NA)
    })
  } else {
    result_list[[target_var]] <- c("SampleSize" = sample_size, "Df" = NA,"SumOfSqs" =NA,"F" =NA,"R
2" = NA, "Pr" = NA)
  }
}

# Output results
result_df <- as.data.frame(do.call(rbind, result_list))
result_df$Pr <- as.numeric(result_df$Pr)

# Add qvalue column, using BH method for multiple testing correction
result_df$qvalue <- p.adjust(result_df$Pr, method = "BH")
write.table(result_df, file = args[3], sep = "\t", quote = FALSE, col.names = NA)
