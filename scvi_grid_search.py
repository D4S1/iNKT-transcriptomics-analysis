import scanpy as sc
import scvi
import pandas as pd
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
import numpy as np
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
        return round(np.mean(kl_scores),3)
    elif mode == "max":
        return round(np.max(kl_scores),3)
    else:
        raise ValueError("Invalid mode. Use 'mean' or 'max'.")
    

def run_scvi_grid_search(
    adata,
    batch_key: str,
    labels_key: str,
    n_layers_list: list,
    latent_dim_list: list,
    resolutions: list,
    min_cells: int,
    latent_dir: str = "latent_embeddings",
    output_csv: str = "grid_search_scvi_metrics.csv"
):

    # Setup anndata
    scvi.model.SCVI.setup_anndata(adata, batch_key="donor", labels_key="treatment", categorical_covariate_keys=["pool"])

    results = []

    for n_layers in n_layers_list:
        for n_latent in latent_dim_list:
            print(f"Training model with n_layers={n_layers}, n_latent={n_latent} ...")
            
            # Initialize and train SCVI
            model = scvi.model.SCVI(
                adata,
                n_layers=n_layers,
                n_latent=n_latent,
                gene_likelihood="zinb",
                dispersion="gene-label"
            )
            model.train(max_epochs=200, early_stopping=True)

            # Get latent representation
            latent = model.get_latent_representation()
            adata.obsm["scvi"] = latent

            # Save latent embedding
            latent_df = pd.DataFrame(latent, index=adata.obs_names)
            latent_df.to_csv(f"{latent_dir}/latent_embedding_nlayers{n_layers}_nlatent{n_latent}.csv")

            # Compute neighbors
            sc.pp.neighbors(adata, use_rep="scvi", key_added="scvi_neighbors")

            for res in resolutions:
                leiden_key = f"scvi_leiden_{res:.2f}"
                sc.tl.leiden(adata, neighbors_key="scvi_neighbors", resolution=res, key_added=leiden_key)

                labels = adata.obs[leiden_key].astype(int)
                X = latent

                # Count small clusters
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
                    'dn_im_mean': label_imbalance(labels, adata.obs['donor'], 'mean'),
                    'dn_im_max': label_imbalance(labels, adata.obs['donor'], 'max'),
                    'pool_im_mean': label_imbalance(labels, adata.obs['pool'], 'mean'),
                    'pool_im_max': label_imbalance(labels, adata.obs['pool'], 'max'),
                }
                results.append(result)

    # Save and return results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_csv, index=False)
    print(f"Grid search finished. Results saved to {output_csv}")
    return results_df



if __name__ == "__main__":

    adata = sc.read("hvg.h5ad")

    # Parameters
    n_layers_list = [2, 5, 8, 10, 12]
    latent_dim_list = [10, 20, 30, 40, 50]
    resolutions = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6]

    results_df = run_scvi_grid_search(
        adata=adata,
        batch_key="donor_treatment",
        labels_key="treatment",
        n_layers_list=n_layers_list,
        latent_dim_list=latent_dim_list,
        resolutions=resolutions,
        min_cells=500,
        latent_dir="run_zinb_cov/latent_embeddings",
        output_csv="run_zinb_cov/grid_search_results.csv"
    )
