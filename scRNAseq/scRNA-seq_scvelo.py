
import scvelo as scv
import scanpy as sc
import sys

if __name__ == '__main__':
    path = sys.argv[1]
    adata = scv.read(f"{path}/Hocine_prevelo_diet.h5ad")
    s = scv.read(f"{path}/Hocine_prevelo_diet_spliced.csv")
    u = scv.read(f"{path}/Hocine_prevelo_diet_unspliced.csv")
    # It would be better to check cell names and gene names
    adata.layers['spliced'] = s.X
    adata.layers['unspliced'] = u.X

    adata_BADC = adata[adata.obs.Genotype == "BADC"].copy()

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=4000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    scv.tl.recover_dynamics(adata)
    scv.tl.latent_time(adata)
    adata.write(f"{path}/Hocine_velo_diet.h5ad", compression='gzip')

    # Do the same with BADC only:
    scv.pp.filter_and_normalize(adata_BADC, min_shared_counts=20, n_top_genes=4000)
    sc.pp.pca(adata_BADC)
    sc.pp.neighbors(adata_BADC, n_pcs=30, n_neighbors=30)
    scv.pp.moments(adata_BADC, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata_BADC)
    scv.tl.velocity_graph(adata_BADC)
    scv.tl.recover_dynamics(adata_BADC)
    scv.tl.latent_time(adata_BADC)
    adata_BADC.write(f"{path}/Hocine_velo_diet_BADC.h5ad", compression='gzip')
