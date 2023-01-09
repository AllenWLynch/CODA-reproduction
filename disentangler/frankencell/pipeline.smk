
test_params = [] # NK_ratio, B_ratio, k, seed, method
#methods = ['mira', 'mira-notune','scvi-notune-3','scvi-notune-4',
#            'scvi-notune-5','scvi-notune-6', 'harmony','scanorama']
methods = ['scvi-notune-3','scvi-notune-4','scvi-notune-5','scvi-notune-6','scanorama']
# single-arm depletion, deplete NK arm from control
seed = 0
for r in [1/2, 1/6, 1/18, 0.]:
    for k in [0., 0.05, 0.1, 0.15]:
        for method in methods:
            test_params.append((1/2, r, k, seed, method))
        seed += 1

# double depletion, deplete B arm from KO
for r in [1/6, 1/18, 0.]:
    for k in [0., 0.05, 0.1, 0.15]:
        for method in methods:
            test_params.append((r, 0., k, seed, method))
        seed += 1

DATA_DIR = 'data/frankencell'
RESULTS_FILE = DATA_DIR + '/results_summary/{biobatch}_{method}/{NK}-{B}-{K}-{seed}.tsv'

outputs = [
    RESULTS_FILE.format(method = method, NK = NK, B = B, K = k, seed = seed, biobatch = biobatch)
    for NK, B, k, seed, method in test_params
    for biobatch in ['bio','batch']
]

rule all:
    input : outputs

rule gen_test:
    output : DATA_DIR + '/test/{NK}-{B}-{K}-{seed}.h5'
    params:
        K=lambda w : w.K,
        B = lambda w : w.B,
        NK = lambda w : w.NK,
        seed = lambda w : w.seed,
        dset1 = '../BatcheffectFrankencell/data/datasets/s3d10_gex.h5ad',
        dset2 = '../BatcheffectFrankencell/data/datasets/s4d1_gex.h5ad',
    threads : 1
    shell :
        './dt-cmd frankencell-gen-test {params.dset1} {params.dset2} {output} '
        '-k {params.K} -B {params.B} -NK {params.NK} --seed {params.seed} -t {threads}'


rule get_dimred:
    input : 
        rules.gen_test.output[0]
    output: 
        h5=DATA_DIR + '/process/{method}/{NK}-{B}-{K}-{seed}.h5',
        plots=DATA_DIR + '/plots/{method}/{NK}-{B}-{K}-{seed}.png',
        dimred = DATA_DIR + '/dimred/{method}/{NK}-{B}-{K}-{seed}.tsv'
    params:
        method = lambda w : w.method
    threads : 1
    shell:
        './dt-cmd frankencell-{params.method} {input} {output.h5} {output.plots} '
        '{output.dimred} -t {threads}'


rule evaluate_bio_metrics:
    input:
        rules.get_dimred.output.h5
    output:
        h5=DATA_DIR + '/results/{method}/{NK}-{B}-{K}-{seed}.h5',
        tsv=DATA_DIR + '/results_summary/bio_{method}/{NK}-{B}-{K}-{seed}.tsv',
    params:
        method = lambda w : w.method,
        metric = lambda w : 'manhattan' if 'mira' in w.method else 'euclidean'
    shell:
        './lib/frankencell-python/cmd-frankencell eval mira '
        '-t {input} -o {output.h5} -r {output.tsv} -p "threshold=1., metric=\'{params.metric}\'"'


rule evaluate_batch_metrics:
    input:
        rules.get_dimred.output.dimred
    output:
        tsv=DATA_DIR + '/results_summary/batch_{method}/{NK}-{B}-{K}-{seed}.tsv',
    params:
        method = lambda w : w.method,
        metric = lambda w : 'manhattan' if 'mira' in w.method else 'euclidean'
    shell:
        './dt-cmd frankencell-eval-batch {input} {output} --metric {params.metric}'