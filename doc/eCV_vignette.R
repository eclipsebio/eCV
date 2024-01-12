## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = FALSE)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("eCV")

## ---- results='hide'----------------------------------------------------------
library("eCV")

## -----------------------------------------------------------------------------
library("tidyverse")

## -----------------------------------------------------------------------------
set.seed(42)
out <- simulate_data(scenario = 1, n_reps = 4, n_features = 1000)
out$sim_data %>% as.data.frame() %>% 
 mutate(`Features group` = as.character(out$sim_params$feature_group)) %>%
 ggplot(aes(x=`Rep 1`,y=`Rep 2`,color=`Features group`)) +
  geom_point(size=1, alpha=0.5) + 
  scale_color_manual(values = c( "#009CA6" , "#F4364C")) + 
  theme_classic()

## -----------------------------------------------------------------------------
# Define parameters for each method.
params <- list(
  eCV = list(max.ite = 100),
  gIDR = list(
    mu = 2,
    sigma = 1.3,
    rho = 0.8,
    p = 0.7,
    eps = 1e-3,
    max.ite = 50
  ),
  mIDR = list(
    mu = 2,
    sigma = 1.3,
    rho = 0.8,
    p = 0.7,
    eps = 1e-3,
    max.ite = 50
  )
)

# Create a list to store results
results <- NULL
# Loop through methods and calculate reproducibility
for (method in c("eCV", "gIDR", "mIDR")) {
  results <- results %>%
    bind_rows(data.frame(
      value =
        mrep_assessment(
          x = out$sim_data,
          method = method,
          param = params[[method]]
        )$rep_index,
      Method = method,
      group = out$sim_params$feature_group
    ))
}

# Plot results
 results %>% 
 mutate(group = ifelse(group == 1,"FALSE","TRUE")) %>%
 ggplot(aes(x=Method, y = value,fill=group)) + 
 scale_fill_manual(values = c( "#009CA6" , "#F4364C")) + 
 geom_boxplot() + 
 theme_classic() + 
 labs(y="Reproducibility assessment", fill="Reproducible\nfeature")

