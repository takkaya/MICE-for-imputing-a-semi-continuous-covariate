library(data.table)
library(ggplot2)
library(MASS)
library(mice)

expand_imp_where <- function(imp,
                             new_semi_continuous,
                             continuous_part,
                             binary_part) {
  column_names <- c(colnames(imp$where), new_semi_continuous)
  imp$where <- cbind(imp$where ,
                     imp$where[, binary_part] | imp$where[, continuous_part])
  colnames(imp$where) <- column_names
  return(imp)
}

create_imp_iters <- function(imp, new_semi_continuous) {
  n_row <- sum(imp$where[, new_semi_continuous])
  imp$imp[[new_semi_continuous]] <-
    matrix(NA, nrow = n_row, ncol = imp$m)
  return(imp)
}

add_new_semi_continuous_variable_to_imputed_datasets <- function(imp,
                                                                 new_semi_continuous,
                                                                 continuous_part,
                                                                 binary_part) {
  imp <-
    add_new_semi_continuous_variable_to_imputed_datasets_existent(imp, new_semi_continuous, continuous_part, binary_part)
  imp <-
    add_new_semi_continuous_variable_to_imputed_datasets_imputed(imp, new_semi_continuous, continuous_part, binary_part)
  return(imp)
}

add_new_semi_continuous_variable_to_imputed_datasets_existent <-
  function(imp,
           new_semi_continuous,
           continuous_part,
           binary_part) {
    imp$data[, new_semi_continuous] <- (as.integer(levels(imp$data[, binary_part])[imp$data[, binary_part]]) * exp(imp$data[, continuous_part]))
    return(imp)
  }

add_new_semi_continuous_variable_to_imputed_datasets_imputed <-
  function(imp,
           new_semi_continuous,
           continuous_part,
           binary_part) {
    imp <- expand_imp_where(imp, new_semi_continuous, continuous_part, binary_part)
    imp <- create_imp_iters(imp, new_semi_continuous)
    for (i in 1:imp$m) {
      complete_data <- mice::complete(imp, i)
      complete_data[, new_semi_continuous] <-  (as.integer(levels(complete_data[, binary_part])[complete_data[, binary_part]]) * exp(complete_data[, continuous_part]))
      imp$imp[[new_semi_continuous]][, i] <- (complete_data[imp$where[, new_semi_continuous], new_semi_continuous])
    }
    return(imp)
  }



ps_binary <- function(data, exposure, exposures) {
  formula <- as.formula(paste(
    exposures[[exposure]]$dependent,
    " ~ ",
    paste(exposures[[exposure]]$independent, collapse = " + ")
  ))
  return(predict(glm(formula, family="binomial" ,data = data), type="response"))
}

ps <- function(data, exposure, exposures) {
  formula <- as.formula(paste(
    exposures[[exposure]]$dependent,
    " ~ ",
    paste(exposures[[exposure]]$independent, collapse = " + ")
  ))
  return(predict(lm(formula, data = data)))
}


#Function for regression
two_part <- function(data, exposure, exposures) {
  formula <- as.formula(paste(
    exposures[[exposure]]$dependent_one, " ~ ",
    paste(exposures[[exposure]]$independent, collapse = " + ")
  ))
  logistic.fit <-
    glm(formula, family = binomial("logit"), data = data)
  formula <- as.formula(paste(
    exposures[[exposure]]$dependent_two,
    " ~ ",
    paste(exposures[[exposure]]$independent, collapse = " + ")
  ))
  logols.fit <-
    lm(formula,  data = data)
  s_d <- sd(logols.fit$residuals)
  data$"(Intercept)" <- 1
  prediction <- 0
  for (name in names(coef(logols.fit))) {
    prediction <- prediction + coef(logols.fit)[name] * data[name]
  }
  phat <- predict(logistic.fit, data = data, type = "response")
  pred_c <-exp(prediction + 1 / 2 * (s_d ^ 2))
  pred <- as.numeric(pred_c[, 1])*phat
  return(pred)
}




ps_on_mice_imputation <- function(imp,
                                  exposures,
                                  multipliers = NULL,
                                  keep_missings = NULL) {
  #
  # Expand the imputed data structure produced by the MICE package
  # by creating new variables, one for each propensity score, marking
  # them as all observations being missing, and inserting the
  # propensity scores for each of the complete data frames.
  #
  # To do that the imp$imp, imp$data, imp$where, and imp$nmis where
  # modified to accomodate the newly introduced values in the data.
  #
  # Other functions can use these new structure, and will not notice
  # that we changed it before using it. So functions like `with()`
  # or `complete()` will work as expected.
  #
  for (exposure in names(exposures)) {
    imputed_variables <-
      data.frame(row.names = nrow(mice::complete(imp, 1)))
    for (i in 1:imp$m) {
      imputed_variables <- cbind(imputed_variables,
                                 # two_part(mice::complete(imp, i), exposure, exposures))
                                 ps(mice::complete(imp, i), exposure, exposures))
                                 # ps_binary(mice::complete(imp, i), exposure, exposures))
    }
    imp$imp[[exposure]] <- data.frame(imputed_variables)
    colnames(imp$imp[[exposure]]) <- 1:imp$m
    
  }
  if (!is.null(multipliers)) {
    for (mult in multipliers) {
      name <- gsub('mult_', '', mult)
      
      # Change NA values
      m_imputed <- imp$data[[mult]][is.na(imp$data[[name]])]
      if (length(m_imputed) != nrow(imp$imp[[name]])) {
        print('!!!')
        print(mult)
        print(str(m_imputed))
        print(str(imp$imp[[name]]))
        print('!!!')
        #stop("Dimensions are not equal.")
      }
      imp$imp[[name]][,] <- imp$imp[[name]][,] * m_imputed
      
      # Change non-NA values
      imp$data[[name]][!is.na(imp$data[[name]])] <- (imp$data[[name]][!is.na(imp$data[[name]])] *
                                                       imp$data[[mult]][!is.na(imp$data[[name]])])
    }
  }
  if (!is.null(keep_missings)) {
    for (missing in keep_missings) {
      imp$imp[[missing]][,] <- NA
    }
  }
  for (exposure in names(exposures)) {
    imp$data[,exposure] <- NA
    col_names <- colnames(imp$where)
    col_names <- c(col_names, exposure)
    imp$where <- cbind(imp$where, TRUE)
    colnames(imp$where) <- col_names
    imp$nmis[[exposure]] <- nrow(imp$data)
  }
  return(imp)
  
}

ps_on_mice_imputation2 <- function(imp,
                                  exposures,
                                  multipliers = NULL,
                                  keep_missings = NULL) {

  for (exposure in names(exposures)) {
    imputed_variables <-
      data.frame(row.names = nrow(mice::complete(imp, 1)))
    for (i in 1:imp$m) {
      imputed_variables <- cbind(imputed_variables,
                                 two_part(mice::complete(imp, i), exposure, exposures))
    }
    imp$imp[[exposure]] <- data.frame(imputed_variables)
    colnames(imp$imp[[exposure]]) <- 1:imp$m
    
  }
  if (!is.null(multipliers)) {
    for (mult in multipliers) {
      name <- gsub('mult_', '', mult)
      
      # Change NA values
      m_imputed <- imp$data[[mult]][is.na(imp$data[[name]])]
      if (length(m_imputed) != nrow(imp$imp[[name]])) {
        print('!!!')
        print(mult)
        print(str(m_imputed))
        print(str(imp$imp[[name]]))
        print('!!!')
        #stop("Dimensions are not equal.")
      }
      imp$imp[[name]][,] <- imp$imp[[name]][,] * m_imputed
      
      # Change non-NA values
      imp$data[[name]][!is.na(imp$data[[name]])] <- (imp$data[[name]][!is.na(imp$data[[name]])] *
                                                       imp$data[[mult]][!is.na(imp$data[[name]])])
    }
  }
  if (!is.null(keep_missings)) {
    for (missing in keep_missings) {
      imp$imp[[missing]][,] <- NA
    }
  }
  for (exposure in names(exposures)) {
    imp$data[,exposure] <- NA
    col_names <- colnames(imp$where)
    col_names <- c(col_names, exposure)
    imp$where <- cbind(imp$where, TRUE)
    colnames(imp$where) <- col_names
    imp$nmis[[exposure]] <- nrow(imp$data)
  }
  return(imp)
  
}


ps_on_mice_imputation_binary <- function(imp,
                                   exposures,
                                   multipliers = NULL,
                                   keep_missings = NULL) {
  
  for (exposure in names(exposures)) {
    imputed_variables <-
      data.frame(row.names = nrow(mice::complete(imp, 1)))
    for (i in 1:imp$m) {
      imputed_variables <- cbind(imputed_variables,
                                 ps_binary(mice::complete(imp, i), exposure, exposures))
    }
    imp$imp[[exposure]] <- data.frame(imputed_variables)
    colnames(imp$imp[[exposure]]) <- 1:imp$m
    
  }
  if (!is.null(multipliers)) {
    for (mult in multipliers) {
      name <- gsub('mult_', '', mult)
      
      # Change NA values
      m_imputed <- imp$data[[mult]][is.na(imp$data[[name]])]
      if (length(m_imputed) != nrow(imp$imp[[name]])) {
        print('!!!')
        print(mult)
        print(str(m_imputed))
        print(str(imp$imp[[name]]))
        print('!!!')
        #stop("Dimensions are not equal.")
      }
      imp$imp[[name]][,] <- imp$imp[[name]][,] * m_imputed
      
      # Change non-NA values
      imp$data[[name]][!is.na(imp$data[[name]])] <- (imp$data[[name]][!is.na(imp$data[[name]])] *
                                                       imp$data[[mult]][!is.na(imp$data[[name]])])
    }
  }
  if (!is.null(keep_missings)) {
    for (missing in keep_missings) {
      imp$imp[[missing]][,] <- NA
    }
  }
  for (exposure in names(exposures)) {
    imp$data[,exposure] <- NA
    col_names <- colnames(imp$where)
    col_names <- c(col_names, exposure)
    imp$where <- cbind(imp$where, TRUE)
    colnames(imp$where) <- col_names
    imp$nmis[[exposure]] <- nrow(imp$data)
  }
  return(imp)
  
}


