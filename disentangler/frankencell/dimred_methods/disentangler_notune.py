
import mira
assert mira.__file__ == '/liulab/alynch/projects/multiomics/BatchEffect/MIRA/mira/__init__.py'

from scipy import sparse
import shutil
import frankencell as fc
import scanpy as sc
from .utils import read_and_process, plot_umaps
import os
import optuna

def run_mira(dynframe, out_h5, plot_file, threads = 1):

    shutil.copy(dynframe, out_h5)

    data = read_and_process(out_h5)
    data.layers['sparse_counts'] = sparse.csr_matrix(data.layers['counts'])

    model = mira.topics.TopicModel(
        *data.shape,
        feature_type='expression',
        exogenous_key='highly_variable',
        counts_layer='sparse_counts',
        categorical_covariates='batch',
        cost_beta = 2.
    )

    model.set_learning_rates(3e-3, 0.25)

    model.set_params(num_topics = 6, decoder_dropout = 0.05, 
                 seed = 2556 + 22).fit(data)

    model.predict(data)
    model.get_umap_features(data, box_cox=0.33)

    sc.pp.neighbors(data, use_rep='X_umap_features', metric='manhattan')
    sc.tl.umap(data, min_dist=0.1)

    plot_umaps(data, plot_file)

    fc.add_dimred_prior(out_h5, data.obsm['X_umap_features'])

def main(args):
    
    run_mira(
        args.dynframe,
        args.outh5,
        args.plotfile,
        threads = args.threads,
    )

def add_arguments(parser):
    parser.add_argument('dynframe', type = str)
    parser.add_argument('outh5', type = str)
    parser.add_argument('plotfile', type = str)
    parser.add_argument('--threads','-t', type = int, default = 1)