import seaborn as sns
from mira.plots.base import map_colors


def disentanglement_plot(adata, ax, gene, color = None, palette = None, size = 0.5,
                        add_legend = False, hue_order = None):
    
    if color is None:
        color = gene
    
    sns.scatterplot(
        x = adata.obs_vector(gene, layer = 'batch_effect'),
        y = adata.obs_vector(gene, layer = 'imputed'),
        c = map_colors(
            ax, adata.obs_vector(color).reshape(-1), palette,
            add_legend=add_legend, hue_order = hue_order,
        ),
        s = size,
        ax = ax,
    )
    
    ax.set(yscale = 'log', xlabel = 'Batch Effect', ylabel = 'Imputed Expression',
        yticks = [], xticks = [0])

    sns.despine()

    return ax