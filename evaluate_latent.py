import scanpy as sc
import pandas as pd
import numpy as np
import os
import re
import argparse
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from collections import Counter


def label_imbalance(cluster_labels, batch_labels, mode="mean", epsilon=1e-8):
    df = pd.concat([cluster_labels, batch_labels], axis=1)
    df.columns = ['cluster', 'batch']

    global_batch_distribution = df['batch'].value_counts(normalize=True).sort_index()
    kl_scores = []

    for cluster, cluster_df in df.groupby('cluster'):
        cluster_batch_distribution = cluster_df['batch'].value_counts(normalize=True).sort_index()

        p = cluster_batch_distribution + epsilon
        q = global_batch_distribution + epsilon
        kl = np.sum(p * np.log(p / q))
        kl_scores.append(kl)

    if mode == "mean":
        return round(np.mean(kl_scores), 3)
    elif mode == "max":
        return round(np.max(kl_scores), 3)
    else:
        raise ValueError("Invalid mode. Use 'mean' or 'max'.")


def evaluate_latents_from_csv(
    adata,
    latent_dir: str,
    resolutions: list,
    donor_key: str,
    pool_key: str,
    min_cells: int,
    output_csv: str
):
    results = []
    pattern = re.compile(r"latent_embedding_nlayers(\d+)_nlatent(\d+)\.csv")

    for fname in os.listdir(latent_dir):
        match = pattern.match(fname)
        if not match:
            continue

        n_layers = int(match.group(1))
        n_latent = int(match.group(2))
        print(f"Evaluating latent embedding: n_layers={n_layers}, n_latent={n_latent}")

        latent_df = pd.read_csv(os.path.join(latent_dir, fname), index_col=0)
        adata.obsm["scvi"] = latent_df.loc[adata.obs_names].values  # Ensure order consistency

        sc.pp.neighbors(adata, use_rep="scvi", key_added="scvi_neighbors")

        for res in resolutions:
            leiden_key = f"scvi_leiden_{res:.2f}"
            sc.tl.leiden(adata, neighbors_key="scvi_neighbors", resolution=res, key_added=leiden_key)

            labels = adata.obs[leiden_key].astype(int)
            X = adata.obsm["scvi"]
            cluster_sizes = Counter(labels)
            n_small_clusters = sum(1 for size in cluster_sizes.values() if size < min_cells)

            result = {
                'n_layers': n_layers,
                'n_latent': n_latent,
                'resolution': res,
                'n_clusters': len(cluster_sizes),
                'n_small_clusters': n_small_clusters,
                'silhouette': silhouette_score(X, labels),
                'calinski_harabasz': calinski_harabasz_score(X, labels),
                'davies_bouldin': davies_bouldin_score(X, labels),
                'dn_im_mean': label_imbalance(labels, adata.obs[donor_key], 'mean'),
                'dn_im_max': label_imbalance(labels, adata.obs[donor_key], 'max'),
                'pool_im_mean': label_imbalance(labels, adata.obs[pool_key], 'mean'),
                'pool_im_max': label_imbalance(labels, adata.obs[pool_key], 'max'),
            }
            results.append(result)

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_csv, index=False)
    print(f"Evaluation complete. Results saved to {output_csv}")
    return results_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate precomputed SCVI latent embeddings.")
    parser.add_argument("--latent", type=str, required=True, help="Path to directory with latent CSV files.")
    parser.add_argument("--output", type=str, default="grid_search_results.csv", help="Path to output CSV file.")
    args = parser.parse_args()

    adata = sc.read("hvg.h5ad")

    resolutions = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6]

    evaluate_latents_from_csv(
        adata=adata,
        latent_dir=args.latent,
        resolutions=resolutions,
        donor_key="donor",
        pool_key="pool",
        min_cells=500,
        output_csv=args.output
    )
