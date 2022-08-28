import frankencell as fc

from scipy import sparse
import networkx as nx

def make_graph(k):
    G = nx.DiGraph()
    G.add_node('HSC', mixing_weights = [1,0,0,0,0,0,0])
    G.add_node('1', mixing_weights = [0.5, 0.125, 0.125, 0.125, 0, 0.125,0])
    
    G.add_node('3', mixing_weights = [0.2, 0.4, 0, 0, 0, 0.4,0])
    G.add_node('Mono', mixing_weights = [0,1-k,0,0,0,k/2,k/2])
    G.add_node('pDC', mixing_weights = [0,k,0,0,0,(1-k)/2,(1-k)/2])
    
    
    G.add_node('2', mixing_weights = [0.2,0,0.4,0.4,0,0,0])
    G.add_node('B cell', mixing_weights = [0,0,1-k,k/2,k/2,0,0])
    G.add_node('NK', mixing_weights = [0,0,k,(1-k)/2,(1-k)/2,0,0])
    
    G.add_edge('HSC','1', weight = 1)
    G.add_edge('1','2', weight = 3)
    G.add_edge('1','3', weight = 2)
    G.add_edge('3','Mono', weight = 1)
    G.add_edge('3','pDC', weight = 1)
    G.add_edge('2','B cell', weight = 1)
    G.add_edge('2','NK', weight = 1)
    
    return G


def get_cell_info(G, NK_ratio, B_ratio, **generation_kwargs):
    
    wt_graph = G.copy()
    wt_graph.edges['2','NK']['weight'] = NK_ratio
    
    ko_graph = G.copy()
    ko_graph.edges['2','B cell']['weight'] = B_ratio
    
    wt_cell_info = fc.fill_scaffold(wt_graph, 'HSC', **generation_kwargs)
    ko_cell_info = fc.fill_scaffold(ko_graph, 'HSC', **generation_kwargs)
    
    return wt_cell_info, ko_cell_info

def generate_batched_perturbation(
    wt_refdata, ko_refdata, outfile,
    seed = 2556,
    k = 0.,
    NK_ratio = 0.,
    B_ratio = 0.,
    generation_kwargs = dict(),
    test_kwargs = dict(),
):
    G = make_graph(k)
    fc.check_definition(G)

    wt_info, ko_info = get_cell_info(G, NK_ratio, B_ratio, **generation_kwargs)

    wt_expression, wt_feature_info = fc.mix_expression(
        cell_info= wt_info,
        datasets = [wt_refdata],
        seed = seed,
        **test_kwargs
    )

    ko_expression, ko_feature_info = fc.mix_expression(
        cell_info=ko_info,
        datasets = [ko_refdata],
        seed = seed + 1,
        **test_kwargs
    )

    assert wt_feature_info == ko_feature_info

    cell_info = fc.append_cell_info(G, wt_info, ko_info)
    expression = sparse.vstack([wt_expression, ko_expression])

    fc.format_dynverse_dataset(
        G = G, 
        cell_info = cell_info,
        feature_info = wt_feature_info, 
        expression = expression, 
        root_node = 'HSC', 
        output_path = outfile,
    )
    
    return G, cell_info, wt_feature_info, expression


def gen_test(*,
    k, B_ratio, NK_ratio, 
    dset1, dset2, 
    outfile,
    seed = 2556,
    threads = 1):

    generation_kwargs = dict(n_cells=2000, ptime_beta= 0.5, seed=2556, sigmoid_aggression=2)

    test_kwargs = dict(
        rd_means = [1500],
        rd_stds = [0.35],
        pure_states = ['HSC', 'CD16+ Mono','B1 B', 'NK', 'CD8+ T','pDC','cDC2','G/M prog'],
        feature_types=['RNA'],
        cell_state_col='cell_type',
        n_jobs=threads,
    )

    generate_batched_perturbation(
        dset1, 
        dset2, 
        outfile,
        seed = seed,
        k = k,
        NK_ratio = NK_ratio,
        B_ratio = B_ratio,
        generation_kwargs = generation_kwargs,
        test_kwargs = test_kwargs,
    )

def main(args):
    
    gen_test(
        k = args.k,
        B_ratio=args.B_ratio,
        NK_ratio = args.NK_ratio,
        dset1 = args.dataset1,
        dset2 = args.dataset2,
        outfile = args.outfile,
        seed= args.seed,
        threads = args.threads,
    )


def add_arguments(parser):

    parser.add_argument('dataset1', type = str)
    parser.add_argument('dataset2', type = str)
    parser.add_argument('outfile', type = str)
    parser.add_argument('--seed', '-s', type = int, default = 0)
    parser.add_argument('--k', '-k', type = float, default = 0.)
    parser.add_argument('--NK_ratio', '-NK', type = float, default = 0.)
    parser.add_argument('--B_ratio', '-B', type = float, default = 0.)
    parser.add_argument('--threads', '-t', type = int, default = 1)