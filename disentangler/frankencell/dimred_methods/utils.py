import frankencell as fc
import scanpy as sc
import matplotlib.pyplot as plt

def read_and_process(dynframe):
    
    testdata = fc.read_dynframe(dynframe)
    
    sc.pp.filter_genes(testdata, min_cells = 15)
    testdata.raw = testdata

    sc.pp.normalize_total(testdata, target_sum=1e4)
    sc.pp.log1p(testdata)
    sc.pp.highly_variable_genes(testdata)

    testdata.layers['counts'] = testdata.raw.to_adata().X.copy()

    sc.pp.highly_variable_genes(testdata, min_disp=0.5)

    return testdata


def plot_umaps(data, outfile):

    sc.pl.umap(data, 
          color = ['mix_weight_' + str(i) for i in range(5)] \
                + ['pseudotime', 'batch'],
          frameon=False, show = False)
    
    plt.savefig(outfile)