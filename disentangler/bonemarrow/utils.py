from . import config
import scanpy as sc

def load_GEX_data():

    data = sc.read_h5ad(config.RAW_DATA)

    rna_data = data[:, data.var.feature_types == 'GEX'].copy()
    del data
    sc.pp.calculate_qc_metrics(rna_data, inplace=True, log1p=False)

    sc.pp.filter_cells(rna_data, min_genes = 400)
    sc.pp.filter_genes(rna_data, min_cells=30)
    sc.pp.normalize_total(rna_data, target_sum=1e4)
    sc.pp.log1p(rna_data)

    sc.pp.highly_variable_genes(rna_data, min_disp= 0.3)
    rna_data.var['exog'] = rna_data.var.highly_variable.copy()
    
    rna_data = rna_data[:, rna_data.var.exog].copy()

    return rna_data


def load_ATAC_data():

    data = sc.read_h5ad('data/bonemarrow/bonemarrow_data.h5ad')

    atac_data = data[:, data.var.feature_types != 'GEX'].copy()
    del data
    sc.pp.filter_cells(atac_data, min_genes = 400)
    sc.pp.filter_genes(atac_data, min_cells=30)

    return atac_data