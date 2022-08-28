
test_params = [] # NK_ratio, B_ratio, k, seed, method

# single-arm depletion, deplete NK arm from control
seed = 2556
for method in ['mira', 'mira-notune']:
    for r in [1/2, 1/6, 1/18, 0.]:
        for k in [0., 0.05, 0.1, 0.15]:
            test_params.append((1/2, r, k, seed, method))

    # double depletion, deplete B arm from KO
    for r in [1/6, 1/18, 0.]:
        for k in [0., 0.05, 0.1, 0.15]:
            test_params.append((r, 0., k, seed, method))

DATA_DIR = 'data/frankencell'
RESULTS_FILE = DATA_DIR + '/results/{method}/{NK}-{B}-{K}-{seed}.h5'

outputs = [
    RESULTS_FILE.format(method = method, NK = NK, B = B, K = k, seed = seed)
    for NK, B, k, seed, method in test_params
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
    threads : 5
    shell :
        './dt-cmd frankencell-gen-test {params.dset1} {params.dset2} {output} '
        '-k {params.K} -B {params.B} -NK {params.NK} --seed {params.seed} -t {threads}'


rule get_dimred:
    input : 
        rules.gen_test.output[0]
    output: 
        h5=DATA_DIR + '/process/{method}/{NK}-{B}-{K}-{seed}.h5',
        plots=DATA_DIR + '/plots/{method}/{NK}-{B}-{K}-{seed}.png'
    params:
        method = lambda w : w.method
    threads : 1
    shell:
        './dt-cmd frankencell-{params.method} {input} {output.h5} {output.plots} -t {threads}'


rule evaluate:
    input:
        rules.get_dimred.output.h5
    output:
        h5=RESULTS_FILE,
        tsv=DATA_DIR + '/results_summary/{method}/{NK}-{B}-{K}-{seed}.tsv',
    params:
        method = lambda w : w.method
    shell:
        './lib/frankencell-python/cmd-frankencell eval mira '
        '-t {input} -o {output.h5} -r {output.tsv} -p "threshold=1."'