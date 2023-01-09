import seaborn as sns
from mira.plots.base import map_colors
import numpy as np

def disentanglement_plot(adata, ax, gene, color = None, palette = None, size = 0.5, alpha = 1.,
                        add_legend = False, hue_order = None, vmin = None, vmax = None):
    
    if color is None:
        color = gene
    
    color_vec = adata.obs_vector(color).reshape(-1)
    #order = np.argso#np.argsort(color_vec)
    x = adata.obs_vector(gene, layer = 'batch_effect')#[order]
    y = adata.obs_vector(gene, layer = 'imputed')#[order]

    sns.scatterplot(
        x = x,
        y = y,
        alpha = alpha,
        c = map_colors(
            ax, color_vec, 
            palette,
            add_legend=add_legend, hue_order = hue_order, vmin = vmin, vmax = vmax,
        ),
        linewidth=0,
        s = size,
        ax = ax,
    )
    
    ax.set(yscale = 'log', xlabel = 'Batch Effect', ylabel = 'Imputed Expression',
        yticks = [], xticks = [0])

    sns.despine()

    return ax