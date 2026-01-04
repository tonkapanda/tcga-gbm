library(tidyverse)
library(tidymodels)
library(kernlab)

set.seed(1)

train_data <- read_csv("data/processed/lasso_train_data.csv") |>
  mutate(MGMT_status = as.factor(MGMT_status))

test_data <- read_csv("data/processed/lasso_test_data.csv") |>
  mutate(MGMT_status = as.factor(MGMT_status))

# tune for cost (penalty for misclassification) and kernel smoothness

svm_spec <- svm_rbf(
  cost = tune(), 
  rbf_sigma = tune()
) |>
  set_engine("kernlab") |>
  set_mode("classification")

svm_recipe <- recipe(MGMT_status ~ ., data = train_data) |>
  update_role(patient_id, new_role = "id") |>
  step_normalize(all_predictors())

svm_workflow <- workflow() |>
  add_recipe(svm_recipe) |>
  add_model(svm_spec)

# hyperparameter tuning

cv_folds <- vfold_cv(train_data, v = 5, strata = MGMT_status)

svm_grid <- grid_regular(
  cost(),
  rbf_sigma(),
  levels = 5
)

svm_results <- svm_workflow |>
  tune_grid(
    resamples = cv_folds,
    grid = svm_grid,
    metrics = metric_set(accuracy, roc_auc)
  )

# tuning visualization 

svm_tuning_plot <- autoplot(svm_results) +
  theme_minimal() +
  labs(title = "svm parameter tuning: cost vs. sigma",
       subtitle = "finding the optimal margin for mgmt classification")

ggsave("plots/svm-model/svm_tuning_plot.png", plot = svm_tuning_plot, width = 8, height = 6)

# evaluation on test set

best_svm <- select_best(svm_results, metric = "roc_auc")
final_svm_wf <- finalize_workflow(svm_workflow, best_svm)
final_svm_fit <- fit(final_svm_wf, data = train_data)

test_results <- augment(final_svm_fit, test_data)

test_metrics <- test_results |>
  metrics(truth = MGMT_status, estimate = .pred_class)

test_auc <- test_results |>
  roc_auc(truth = MGMT_status, .pred_methylated)

print(bind_rows(test_metrics, test_auc))

# conf matrix + roc curve visualization

cm_plot <- test_results |>
  conf_mat(truth = MGMT_status, estimate = .pred_class) |>
  autoplot(type = "heatmap") +
  labs(title = "svm confusion matrix", subtitle = "predicted vs. actual mgmt status")

ggsave("plots/svm-model/svm_confusion_matrix.png", plot = cm_plot, width = 6, height = 5)

roc_curve_plot <- test_results |>
  roc_curve(truth = MGMT_status, .pred_methylated) |>
  autoplot() +
  labs(title = "svm roc curve", subtitle = "sensitivity vs. specificity for mgmt prediction")

ggsave("plots/svm-model/svm_roc_curve.png", plot = roc_curve_plot, width = 6, height = 6)

# metrics

svm_final_table <- bind_rows(test_metrics, test_auc) |>
  select(Metric = .metric, Estimate = .estimate) |>
  mutate(
    Metric = case_when(
      Metric == "accuracy" ~ "accuracy",
      Metric == "kap" ~ "kappa (agreement)",
      Metric == "roc_auc" ~ "roc-auc (separation power)"
    ),
    Estimate = round(Estimate, 4)
  )
print(svm_final_table)

write_csv(svm_final_table, "plots/svm-model/svm_performance_summary.csv")