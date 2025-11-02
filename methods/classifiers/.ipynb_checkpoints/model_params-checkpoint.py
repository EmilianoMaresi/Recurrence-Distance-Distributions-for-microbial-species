# model_params.py

# ========== SMALL dataset: less than 1000 samples and 1500 features ==========
SMALL_RF_PARAMS = {
    "n_estimators": 100,
    "max_depth": 5,
    "max_features": "sqrt",
    "min_samples_split": 5,
    "random_state": 42
}

# NOTE: "tree_metod": "hist" is a faster tree building, uses histogram-based split finding and it's usually the preferred choice for structured/tabular data.
SMALL_XGB_PARAMS = {
    "n_estimators": 150,
    "max_depth": 3,
    "learning_rate": 0.05,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "reg_lambda": 1.0,
    "random_state": 42,
    "eval_metric": "mlogloss",
    "tree_method": "hist"    	
}

SMALL_MLP_PARAMS = {
    "hidden_layer_sizes": (50,),
    "alpha": 1e-3,
    "activation": "relu",
    "max_iter": 300,
    "early_stopping": True,
    "random_state": 42
}


# ========== LARGER dataset: more than 1000 samples and 1500 features ==========
LARGE_RF_PARAMS = {
    "n_estimators": 400,
    "max_depth": None,
    "max_features": "sqrt",
    "random_state": 42
}

# NOTE: "tree_metod": "hist" is a faster tree building, uses histogram-based split finding and it's usually the preferred choice for structured/tabular data.
LARGE_XGB_PARAMS = {
    "n_estimators": 800,
    "max_depth": 6,
    "learning_rate": 0.05,
    "subsample": 0.9,
    "colsample_bytree": 0.8,
    "reg_lambda": 1.0,
    "random_state": 42,
    "eval_metric": "mlogloss",
    "tree_method": "hist"
}

LARGE_MLP_PARAMS = {
    "hidden_layer_sizes": (200, 100),
    "alpha": 1e-4,
    "activation": "relu",
    "max_iter": 400,
    "early_stopping": True,
    "random_state": 42
}

