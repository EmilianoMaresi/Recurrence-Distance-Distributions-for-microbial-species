#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np

from time import time
from pathlib import Path
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster

from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score
)

from multiprocessing import Pool, cpu_count


def format_time(t):
    """Convert seconds (float) to a formatted string including days, hours, minutes, 
    and seconds with 4 decimal places."""
    days, remainder = divmod(t, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, remainder = divmod(remainder, 60)
    seconds = remainder  # float with fractional seconds

    # Format seconds to have 2 integer digits and 4 decimal places (e.g., 05.1235)
    sec_str = f"{seconds:07.4f}"

    if days > 0:
        return f"{int(days)}d {int(hours):02d}:{int(minutes):02d}:{sec_str}"
    return f"{int(hours):02d}:{int(minutes):02d}:{sec_str}"



def compute_clustering(X, labels, fragments_min=0, fragments_max=1500,  
                       method="average", compute_silhouette=True):
    """
    Perform hierarchical clustering on a fragment subrange,
    and compute ARI, NMI, and optionally Silhouette scores safely.

    Parameters
    ----------
    X : pandas.DataFrame
        Input data matrix (samples × features).
    fragments_min, fragments_max : int
        Range of fragment indices (inclusive) to select from X.
    labels : array-like
        Ground-truth labels for ARI/NMI computation.
    method : str, optional
        Linkage method (e.g. 'average', 'complete', 'single').
    compute_silhouette : bool, optional
        Whether to compute the silhouette score (default=True).

    Returns
    -------
    dict
        Dictionary with metrics and metadata.
    """
    # Select fragment subrange
    X = np.asarray(X.iloc[:, fragments_min:fragments_max+1])
    labels = np.asarray(labels)

    if X.shape[0] < 2:
        raise ValueError("Not enough valid samples for clustering.")

    # Compute condensed cosine distances
    D_condensed = pdist(X, metric="cosine")
    Z = linkage(D_condensed, method=method)

    n_clusters = len(set(labels))
    clusters = fcluster(Z, t=n_clusters, criterion="maxclust")

    # External validation metrics
    ari = adjusted_rand_score(labels, clusters)
    nmi = normalized_mutual_info_score(labels, clusters)

    # Internal validation metric (optional)
    sil = None
    if compute_silhouette and n_clusters > 1 and len(np.unique(clusters)) > 1:
        try:
            sil = silhouette_score(X, clusters, metric="cosine")
        except Exception:
            sil = np.nan

    # Return all results
    return {
        "fragments_range": (fragments_min, fragments_max),
        "ari": ari,
        "nmi": nmi,
        "silhouette": sil,
        "data_dimensions": X.shape,
    }


def _process_one_file(f, genus_labels):
    """Compute clustering metrics for a single AFLP file."""

    df = pd.read_csv(f)

    # Labels
    if genus_labels:
        aflp_labels = df["species"].str.split(r'[_\s]').str[0]
    else:
        aflp_labels = df["species"]

    # Feature matrix
    aflp_values = df.iloc[:, 5:]

    # Run hierarchical clustering
    res = compute_clustering(X=aflp_values, labels=aflp_labels)
    res["kmer_pair"] = f.stem

    return res


def _worker_wrapper(args):
    """Wrapper so multiprocessing can pickle arguments cleanly."""
    #f, binarized_aflp, genus_labels = args
    #return _process_one_file(f, binarized_aflp, genus_labels)
    f, genus_labels = args
    return _process_one_file(f, genus_labels)


def aflp_clustering(aflp_folder, nproc=2, genus_labels=False):

    aflp_files_lst = [f for f in aflp_folder.rglob('*') if f.is_file()]
    total = len(aflp_files_lst)

    flag_genus = '_genus' if genus_labels else ''

    print(f"\nClustering AFLPs{flag_genus} using {nproc} processes...")
    print(f"Total AFLP files: {total}\n")

    # args passed to each worker
    args_list = [(f, genus_labels) for f in aflp_files_lst]

    results_list = []
    processed = 0

    # Run parallel pool
    with Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(_worker_wrapper, args_list):
            results_list.append(res)

            processed += 1

            if processed % 50 == 0:
                print(f"Processed: {processed} of {total}")

    # Build results
    results_df = pd.DataFrame(results_list)
    results_df = results_df[
        ["kmer_pair", "ari", "nmi", "silhouette"]
    ]
    results_df = results_df.sort_values(
        by=["ari", "nmi"], ascending=[False, False]
    ).reset_index(drop=True)

    return results_df


def categorize_metrics_with_silhouette(row):
    ari = row['ari']
    nmi = row['nmi']
    sil = row['silhouette']

    if ari >= 0.7 and nmi >= 0.8:
        if sil > 0.2:
            return 'Top Performer (Tight & Accurate)'
        else:
            return 'Top Performer (Geometrically Messy)'
    elif sil < 0:
        return 'Geometrically Broken (Negative Silhouette)'
    elif ari < 0.4 and nmi >= 0.6:
        return 'Over-Splitter (Pure Sub-clusters)'
    elif ari < 0.2 and nmi < 0.3:
        return 'Poor / Random'
    else:
        return 'Moderate / Mixed'


def run_clustering(input_folder, cpu_processes=2):

    print("[[Step 4 - Clustering]]\n")
    
    tic = time()

    #run clustering
    results_clustering = aflp_clustering(aflp_folder=input_folder, nproc = cpu_processes)

    #make kmer pairs clustering categories
    results_clustering['category'] = results_clustering.apply(categorize_metrics_with_silhouette, axis=1)

    # save as csv file
    out_folder = Path(input_folder.parent, "clusters")
    out_folder.mkdir(parents=True, exist_ok=True)

    out_file = Path(out_folder, f"clustering_results.csv")
    results_clustering.to_csv(out_file, index=False)

    print("Saved in:", str(out_file), "\n")

    ### make clusterings performance summary
    # Explicitly define all possible categories (because some may not be present in the results and need to be counted as 0)
    all_categories = [
        'Top Performer (Tight & Accurate)',
        'Top Performer (Geometrically Messy)',
        'Over-Splitter (Pure Sub-clusters)',
        'Moderate / Mixed',
        'Geometrically Broken (Negative Silhouette)',
        'Poor / Random'    
    ]

    # Count occurrences, reindex to include all categories with 0 if missing, and reset index
    summary_df = results_clustering['category'].value_counts().reindex(all_categories, fill_value=0).reset_index()

    # Rename columns
    summary_df.columns = ['category', 'counts']

    # Calculate percentage over the total number of items
    total_items = summary_df['counts'].sum()
    summary_df['percentage'] = (summary_df['counts'] / total_items) * 100
    summary_df['percentage'] = summary_df['percentage'].round(2)

    # Save to CSV
    summary_df.to_csv(Path(out_folder, "clustering_summary.csv"), index=False)

    toc = time()
    print("Clustering elapsed time:", format_time(toc-tic))
    print("Clustering - Done!")



