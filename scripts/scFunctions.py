import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def get_cluster_proportions(adata,
                            cluster_key="cluster_final",
                            sample_key="replicate",
                            drop_values=None):
    """
    Input
    =====
    adata : AnnData object
    cluster_key : key of `adata.obs` storing cluster info
    sample_key : key of `adata.obs` storing sample/replicate info
    drop_values : list/iterable of possible values of `sample_key` that you don't want
    
    Returns
    =======
    pd.DataFrame with samples as the index and clusters as the columns and 0-100 floats
    as values
    """
    
    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([cluster_key, sample_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index() 
    props = props.pivot(columns=sample_key, index=cluster_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)
    
    if drop_values is not None:
        for drop_value in drop_values:
            props.drop(drop_value, axis=0, inplace=True)
    return props


def plot_cluster_proportions(cluster_props, 
                             cluster_palette=None,
                             xlabel_rotation=0): 
    fig, ax = plt.subplots(dpi=300)
    fig.patch.set_facecolor("white")
    
    cmap = None
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette, 
            n_colors=len(cluster_palette), 
            as_cmap=True)
   
    cluster_props.plot(
        kind="bar", 
        stacked=True, 
        ax=ax, 
        legend=None, 
        colormap=cmap
    )
    
    ax.legend(bbox_to_anchor=(1.01, 1), frameon=False, title="")
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=xlabel_rotation)
    ax.set_xlabel(cluster_props.index.name.capitalize())
    ax.set_ylabel("Proportion")
    fig.tight_layout()
    
    return fig



def pp_PCA(a):


    """
    Input
    =====
    a : AnnData object
    
    Returns
    =======
    AnnData object with counts log-normalized and scaled, HVG found and pca computed
    """

    adata = a.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    adata.raw = adata

    adata = adata[:, adata.var.highly_variable]

    # sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    return(adata)



def embedding_density_gene(adata, gene):

    """
    Input
    =====
    adata : AnnData object
    gene : key of `adata.var` that correspond to the wanted gene 


    Show
    =======
    a density embedding graph splitting the gene expression 
    in half across all cells and showing on the umap the density of gene expression

    """

    row = adata.raw[:,gene].X.todense()

    row_percent = []
    med = np.percentile(list(row), 50)
    for i in row:
        if i > med:
            row_percent.append('High')
        else:
            row_percent.append('Low')


    adata.obs[f"{gene}_percent"] = row_percent
    print("____ Computing embedding density ____")
    sc.tl.embedding_density(adata, groupby=f"{gene}_percent")
    sc.pl.embedding_density(adata, groupby=f"{gene}_percent")


    return 0





