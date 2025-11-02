#!/usr/bin/env python
# coding: utf-8


import pandas as pd
from time import time
from datetime import timedelta
import numpy as np
from pathlib import Path
from sklearn.model_selection import StratifiedKFold
from methods.classifiers.random_forest import RandomForestModel
from methods.classifiers.mlp import MLPModel
from methods.classifiers.xgboost_model import XGBoostModel
from sklearn.metrics import accuracy_score, f1_score

#classifiers parameters presets
from methods.classifiers.model_params import (
    RF_PARAMS, XGB_PARAMS, MLP_PARAMS
)


#Available classifiers
CLASSIFIERS = {
    "random_forest": RandomForestModel,
    "mlp": MLPModel,
    "xgboost": XGBoostModel
}


#5-fold Cross Validation
def cross_validate(X, y, classifier_name, folds=5, params=None):

    start = time()

    if classifier_name not in CLASSIFIERS:
        raise ValueError(f"Unknown classifier '{classifier_name}'. Choose from {list(CLASSIFIERS.keys())}.")

    clf_class = CLASSIFIERS[classifier_name]
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=42)

    acc_scores = []
    f1_scores = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y), start=1):
        X_train, X_test = X.iloc[train_idx,:], X.iloc[test_idx,:]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        clf = clf_class(**(params or {}))  # create a model instance
        clf.train(X_train, y_train)

        y_pred = clf.predict(X_test)
        acc = accuracy_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred, average="macro")  # change "macro" to "weighted" if unbalanced dataset

        acc_scores.append(acc)
        f1_scores.append(f1)

    elapsed = time() - start

    return {
        "classifier_name": classifier_name,
        "accuraciey_folds": acc_scores,
        "accuracy_mean": np.mean(acc_scores),
        "accuracy_std": np.std(acc_scores),
        "f1_folds": f1_scores,
        "f1_mean": np.mean(f1_scores),
        "f1_std": np.std(f1_scores),
        "time": elapsed
    }


def result_to_csv_row(results_dict, kmer_pair):
    values = [
        kmer_pair,
        results_dict["classifier_name"],
        str(round(results_dict["accuracy_mean"],4)),
        str(round(results_dict["accuracy_std"],4)),
        str(round(results_dict["f1_mean"],4)),
        str(round(results_dict["f1_std"],4)),
        str(round(results_dict["time"],4))
    ]

    #csv_row = "\t".join(values)+"\n"
    csv_row = ",".join(values)+"\n"
    
    print(values[0],values[1],values[2],values[3],values[4],values[5],values[6])
    
    return csv_row


def get_params(model_name, n_samples):
    #Returns the parameters set for the classifiers
    
    if model_name == "random_forest":
        return RF_PARAMS
    elif model_name == "xgboost":
        return XGB_PARAMS
    elif model_name == "mlp":
        return MLP_PARAMS
    else:
        raise ValueError(f"Unknown model '{model_name}'")

def run_classifiers(input_folder, min_fragment_length, max_fragment_length, cpu_processes, binarized_distributions=False, genus_only_labels=False):

    print("[[Step 3 - Classification]]")
    print("binarized distributions:", binarized_distributions)
    print("minimum fragment length:", min_fragment_length)
    print("maximum fragment length:", max_fragment_length)
    print("cpu limit:", cpu_processes)
    print()

    files = [f for f in input_folder.iterdir()]

    out_path = Path(input_folder.parent, "classifiers")
    out_path.mkdir(parents=True, exist_ok=True)
    
    csv_suffix = "_genus" if genus_only_labels else ""
    
    if binarized_distributions == True:
        out_filename = "classifiers_results_binarized"+csv_suffix+".csv"
    else:
        out_filename = "classifiers_results"+csv_suffix+".csv"

    out_file = open(Path(out_path, out_filename), "w")

    columns = ["kmer_pair", "classifier", "accuracy_mean", "accuracy_std", "f1_mean", "f1_std", "running_time" ]
    out_file.write(",".join(columns)+"\n")

    # Running classifiers
    tic = time()
    for i, file in enumerate(input_folder.iterdir(), start=1):

        kmer_pair = file.stem
        print(i,"/",len(files),kmer_pair)

        data = pd.read_csv(file)
        data.set_index(data.columns[0], inplace=True)

        X = data.iloc[:, (4+min_fragment_length):(4+max_fragment_length+1)]

        if binarized_distributions == True:
            X = (X != 0).astype(int)

        
        if genus_only_labels:
            # use genus as labels. Split by underscore "_" or space " "
            y = data.iloc[:, 1].str.split(r'[_\s]').str[0]
            
        else:
            # use species as label
            y = data.iloc[:, 1]
            

        
        n_samples = X.shape[0]

        
        print("=== Running Random Forest Classifier ===")
        rf_params = get_params("random_forest", n_samples)
        rf_params["n_jobs"] = cpu_processes
        result_random_forest = cross_validate(X, y, "random_forest", params=rf_params)    
        out_file.write(result_to_csv_row(result_random_forest, kmer_pair))
        
        
        print("=== Running MLP Classifier ===")
        mlp_params = get_params("mlp", n_samples)
        result_mlp = cross_validate(X, y, "mlp", params=mlp_params)
        out_file.write(result_to_csv_row(result_mlp, kmer_pair))

        
        print("=== Running XGBoost Classifier ===")
        xgboost_params = get_params("xgboost", n_samples)
        xgboost_params["n_jobs"] = cpu_processes
        result_xgboost = cross_validate(X, y, "xgboost", params=xgboost_params)
        out_file.write(result_to_csv_row(result_xgboost, kmer_pair))
        
        
        print(f"Elapsed time: {str(timedelta(seconds=(time()-tic))).split('.')[0]}\n")

    out_file.flush()
    out_file.close()

    elapsed_time = time() - tic
    print("\nStep 3 - Classification: Done!")
    print(f"Results saved in: {out_file.name}")
    print(f"Running time: {str(timedelta(seconds=elapsed_time)).split('.')[0]}")
    print()


