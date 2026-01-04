library(tidyverse)
library(tidymodels)
library(glmnet)

set.seed(1)

# load, prep, and split data

data <- read_csv("data/processed/gbm_final_dataset.csv") |>
  mutate(MGMT_status = as.factor(if_else(MGMT_status == 1, "methylated", "unmethylated")))

data_split <- initial_split(data, prop = 0.75, strata = MGMT_status)
train_data <- training(data_split)
test_data <- testing(data_split)

# define lasso model and recipe

lasso_spec <- logistic_reg(penalty = tune(), mixture = 1) |>
  set_engine("glmnet")

# recipe and normalization

lasso_recipe <- recipe(MGMT_status ~ ., data = train_data) |>
  update_role(patient_id, new_role = "id") |>
  step_normalize(all_predictors())

lasso_workflow <- workflow() |>
  add_recipe(lasso_recipe) |>
  add_model(lasso_spec)

# cv to find best lambda

val_folds <- vfold_cv(train_data, v = 5, strata = MGMT_status)
lasso_grid <- grid_regular(penalty(), levels = 50)

tune_results <- tune_grid(
  lasso_workflow,
  resamples = val_folds,
  grid = lasso_grid,
  metrics = metric_set(roc_auc)
)

# visualization of tune results

tuning_plot <- autoplot(tune_results) +
  theme_minimal() +
  labs(title = "lasso tuning", caption = "peak auc represents best balance")

# select best penalty

best_penalty <- select_best(tune_results, metric = "roc_auc")

final_lasso <- finalize_workflow(lasso_workflow, best_penalty)
lasso_fit <- fit(final_lasso, data = train_data)

# genes where coefficient != 0

lasso_coefs <- lasso_fit |>
  extract_fit_parsnip() |>
  tidy() |>
  filter(estimate != 0 & term != "(Intercept)") |>
  arrange(desc(abs(estimate)))

coef_plot <- lasso_coefs |>
  slice_max(abs(estimate), n = 20) |>
  ggplot(aes(x = reorder(term, estimate), y = estimate, fill = estimate > 0)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "top 20 genes selected by lasso",
    x = "gene name", y = "impact (coefficient weight)",
    fill = "positive correlation"
  ) +
  theme_minimal()

# save new datasets

selected_genes <- lasso_coefs$term

final_train <- train_data |> select(patient_id, MGMT_status, all_of(selected_genes))
final_test <- test_data |> select(patient_id, MGMT_status, all_of(selected_genes))

write_csv(final_train, "data/processed/lasso_train_data.csv")
write_csv(final_test, "data/processed/lasso_test_data.csv")

print(paste("lasso has selected", length(selected_genes), "genes."))


ggsave("plots/feature-selection/lasso_tuning_plot.png",
  plot = tuning_plot, width = 8, height = 6
)

ggsave("plots/feature-selection/lasso_coefs.png",
  plot = coef_plot, width = 10, height = 8
)
