
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

    def faux_print(*x):
        return 'Trial completed.'

    mira.topic_model.trainer._print_study = faux_print

    train_data = dynframe+'_train'
    test_data = dynframe+'_test'

    if os.path.isdir(train_data):
        shutil.rmtree(train_data)
        shutil.rmtree(test_data)

    train, test = mira.topics.SpeedyTuner.train_test_split(data, 
            train_size=0.8, 
            stratify=data.obs_vector('batch'), seed = 0
        )

    model.write_ondisk_dataset(train, dirname= train_data)
    model.write_ondisk_dataset(test, dirname= test_data)

    del train, test

    try:
        optuna.delete_study(
            storage = 'sqlite:///mira-BENCHMARKING.db',
            study_name = dynframe
        )
    except KeyError:
        pass

    tuner = mira.topics.SpeedyTuner(
        model = model,
        save_name = dynframe,
        min_topics = 3,
        max_topics = 10,
        seed = 2556,
        min_trials = 32,
        max_trials = 64,
        n_jobs = threads,
        stop_condition = 8,
        storage = 'sqlite:///mira-BENCHMARKING.db',
    )

    tuner.fit(train_data, test_data)

    model = tuner.fetch_best_weights()

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