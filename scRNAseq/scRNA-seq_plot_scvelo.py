
import scvelo as scv
import scanpy as sc
import sys

if __name__ == '__main__':
    path = sys.argv[1]
    output_path = sys.argv[2]
    adata = scv.read(f"{path}/Hocine_velo_diet_BADC.h5ad")
    scv.pl.velocity_embedding_stream(
        adata, basis='umap', color='Fate',
        palette=["#1f77b4", "#5254a3", "#e377c2", "#ffbb78", "#e6550d",
                 "#98df8a", "#f7b6d2", "#ff9896", "#9467bd", "#c5b0d5",
                 "#8c564b", "#c49c94", "#bcbd22", "#dbdb8d", "#1bca8d"][::-1],
        size=40,
        save=f"{output_path}/velo_on_umap_BADC.pdf")
