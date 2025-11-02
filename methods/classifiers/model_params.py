
RF_PARAMS = {
    "n_estimators": 100,
    "max_depth": 5,
    "max_features": "sqrt",
    "min_samples_split": 5,
    "random_state": 42
}

# NOTE: "tree_metod": "hist" is a faster tree building, uses histogram-based split finding and it's usually the preferred choice for structured/tabular data.
XGB_PARAMS = {
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

MLP_PARAMS = {
    "hidden_layer_sizes": (50,),
    "alpha": 1e-3,
    "activation": "relu",
    "max_iter": 300,
    "early_stopping": True,
    "random_state": 42
}

