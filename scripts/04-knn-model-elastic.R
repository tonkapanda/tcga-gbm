library(tidymodels)
library(tidyverse)

set.seed(1)

# load data

train_data <- read_csv("data/processed/lasso_elastic_train_data.csv") |>
  mutate(MGMT_status = as.factor(MGMT_status))

test_data <- read_csv("data/processed/lasso_elastic_test_data.csv") |>
  mutate(MGMT_status = as.factor(MGMT_status))

# define model & recipe

knn_spec <- nearest_neighbor(neighbors = tune(), weight_func = tune(), dist_power = tune()) |>
  set_engine("kknn") |>
  set_mode("classification")

knn_recipe <- recipe(MGMT_status ~ ., data = train_data) |>
  update_role(patient_id, new_role = "id") |>
  step_normalize(all_predictors())

knn_workflow <- workflow() |>
  add_recipe(knn_recipe) |>
  add_model(knn_spec)

# cv tuning

cv_folds <- vfold_cv(train_data, v = 5, strata = MGMT_status)
k_grid <- grid_regular(
  neighbors(),
  weight_func(),
  dist_power(),
  levels = 10
)

tune_results <- knn_workflow |>
  tune_grid(
    resamples = cv_folds,
    grid = k_grid,
    metrics = metric_set(accuracy, roc_auc)
  )

# select best model & fit

best_knn <- select_best(tune_results, metric = "roc_auc")
final_knn <- finalize_workflow(knn_workflow, best_knn)
knn_fit <- fit(final_knn, data = train_data)

# test set predictions

pred_class <- predict(knn_fit, test_data, type = "class")
pred_prob <- predict(knn_fit, test_data, type = "prob")

results <- test_data |>
  dplyr::select(MGMT_status) |>
  bind_cols(pred_class, pred_prob)

# metrics (accuracy, roc)

metrics_result <- results |>
  metrics(truth = MGMT_status, estimate = .pred_class, .pred_methylated)

roc_result <- results |>
  roc_auc(truth = MGMT_status, .pred_methylated)

print("Final Test Metrics:")
print(metrics_result)
print(roc_result)

# visualizations

# 1. tuning results plot
knn_tuning_plot <- tune_results |>
  collect_metrics() |>
  filter(.metric == "roc_auc") |>
  ggplot(aes(x = neighbors, y = mean, color = weight_func)) +
  geom_line(size = 1) +
  geom_point() +
  facet_wrap(~dist_power, labeller = label_both) +
  labs(
    title = "knn multi-parameter tuning for elastic",
    subtitle = "comparing distance power (1 = manhattan, 2 = euclidean) and weighting",
    x = "number of neighbors (k)",
    y = "mean roc-auc (cross-validation)"
  ) +
  theme_minimal()

# 2. confusion matrix

predictions <- predict(knn_fit, test_data) |>
  bind_cols(test_data)

conf_mat_result <- predictions |>
  conf_mat(truth = MGMT_status, estimate = .pred_class)

print(conf_mat_result)

cm_plot <- autoplot(conf_mat_result, type = "heatmap") +
  labs(
    title = "knn confusion matrix for elastic",
    subtitle = "predicted vs actual MGMT status"
  )

# 3. prediction confidence

pred_dist_plot <- results |>
  ggplot(aes(x = .pred_methylated, fill = MGMT_status)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "prediction distribution (confidence) for elastic",
    subtitle = "density of predicted probability for methylated class",
    x = "predicted probability (methylated)",
    y = "density",
    fill = "true status"
  ) +
  theme_minimal() +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50")

# save plots

ggsave("plots/knn-model/knn_tuning_plot_elastic.png", plot = knn_tuning_plot, width = 10, height = 8)
ggsave("plots/knn-model/confusion_matrix_elastic.png", plot = cm_plot, width = 6, height = 5)
ggsave("plots/knn-model/prediction_distribution_elastic.png", plot = pred_dist_plot, width = 8, height = 6)
