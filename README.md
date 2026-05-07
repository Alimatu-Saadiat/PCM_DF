# Promotion Time Cure Model (PTCM) with Machine Learning Incidence Modeling

This repository contains R code for implementing and evaluating Promotion Time Cure Models (PTCMs) using several Logit, Spline and Decision Tree approaches for the incidence component and proportional hazards regression for the latency component.

The project was developed using the cutaneous melanoma dataset and includes:

- Logistic Regression PTCM
- Spline-based PTCM
- Decision Tree-based PTCM
- EM Algorithm estimation
- Weighted Breslow baseline survival estimation
- Bootstrap inference
- ROC/AUC evaluation

---

# Author

Alimatu Agongo  
PhD Student in Statistics  
University of Texas at Arlington

---

# Overview

The Promotion Time Cure Model assumes that a proportion of subjects are cured and will never experience the event of interest.

The model separates:

1. Incidence Component
   - Models probability of being uncured
   - Implemented using:
     - Logistic regression
     - Natural splines
     - Decision trees

2. Latency Component
   - Models survival among susceptible subjects
   - Implemented using Cox proportional hazards regression

Estimation is performed using an EM algorithm.

---

# Methods Included

## Logistic Regression PTCM
Functions:
- `em.logit.Pois()`
- `smcure.logit.Pois()`

Features:
- EM estimation
- Cox PH latency model
- Bootstrap standard errors
- ROC/AUC evaluation

---

## Spline-based PTCM
Functions:
- `em.spline.Pois()`
- `smcure.spline.Pois()`

Features:
- Natural spline incidence modeling
- Flexible nonlinear covariate effects
- Bootstrap inference

---

## Decision Tree PTCM
Functions:
- `em.dt.Pois()`
- `smcure.dt.Pois()`

Features:
- Decision tree incidence modeling
- Cross-validated hyperparameter tuning
- Variable importance
- Representative tree construction

---

# Main Statistical Components

## Weighted Breslow Estimator
Function:
- `smsurv()`

Used to estimate the baseline survival function within the EM algorithm.

---

## EM Algorithm

The EM algorithm iteratively updates:

### E-step
Posterior probability of being uncured.

### M-step
- Incidence model parameters
- Latency model parameters
- Baseline survival estimates

---

# Data

The analysis uses the cutaneous melanoma dataset.

Variables used:
- `x2` : Age
- `x3` : Nodule Category
- `x4` : Gender
- `x6` : Tumor Thickness

Outcome variables:
- `t` : Survival time
- `d` : Event indicator

---

# Required R Packages

```r
library(survival)
library(survminer)
library(MASS)
library(e1071)
library(pROC)
library(caTools)
library(caret)
library(rpart)
library(neuralnet)
library(randomForest)
library(randomForestSRC)
library(xgboost)
library(splines)
library(readr)
library(timeROC)
library(rpart.plot)
````

---

# Running the Code

## Load Data

```r
Melanoma <- read.table("melanoma-data-1.txt", header = TRUE)
```

## Split Data

```r
set.seed(2026)

idx <- createDataPartition(Melanoma$d, p = 0.7, list = FALSE)

train_df <- Melanoma[idx, ]
test_df  <- Melanoma[-idx, ]
```

## Run Model

```r
methods <- c("Logit","Spline","DT")

res_test <- plot_test_rocs_ptcm(
  methods,
  train_df,
  test_df,
  B = 500,
  emmax = 1000,
  eps = 1e-3
)
```

---

# Outputs

The code produces:

* Estimated uncure probabilities
* Survival estimates
* ROC curves
* AUC summaries
* Bootstrap standard errors
* Variable importance plots
* Representative decision trees


---

# Decision Tree Hyperparameter Tuning

The decision tree model performs cross-validation over:

* Complexity parameter (`cp`)
* Minimum split size (`minsplit`)


The best model is selected using AUC.

---

# Notes

* Age and tumor thickness are standardized before modeling.
* Baseline survival is estimated separately for training and testing sets.
* Bootstrap procedures may require substantial computation time.

---

# Repository Structure

```r
README.md
PTCM_Code.R
melanoma-data-1.txt
```

---

# Future Extensions

Potential future work includes:

* Random Forest PTCM
* XGBoost PTCM
* Neural Network PTCM


---

# Citation

If using this code, please cite the corresponding research work appropriately.
