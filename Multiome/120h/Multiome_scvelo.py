
import scvelo as scv
import scanpy as sc
import sys

if __name__ == '__main__':
   path = sys.argv[1]
   adata = scv.read(f"{path}/Multiome_prevelo.h5ad")
   scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=4000)
   sc.pp.pca(adata)
   sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
   scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
   scv.tl.velocity(adata)
   scv.tl.velocity_graph(adata)
   scv.tl.recover_dynamics(adata)
   scv.tl.latent_time(adata)
   # see https://stackoverflow.com/questions/70234014/valueerror-index-is-a-reserved-name-for-dataframe-columns
   del adata.raw
   adata.write(f"{path}/Multiome_velo.h5ad", compression='gzip')
   adata = scv.read(f"{path}/Multiome_caudal_prevelo.h5ad")
   scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=4000)
   sc.pp.pca(adata)
   sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
   scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
   scv.tl.velocity(adata)
   scv.tl.velocity_graph(adata)
   scv.tl.recover_dynamics(adata)
   scv.tl.latent_time(adata)
   # see https://stackoverflow.com/questions/70234014/valueerror-index-is-a-reserved-name-for-dataframe-columns
   del adata.raw
   adata.write(f"{path}/Multiome_caudal_velo.h5ad", compression='gzip')
