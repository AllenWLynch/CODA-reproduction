
import argparse
import sys
from functools import partial
import logging

from disentangler.frankencell import generate_KO_test
from disentangler.frankencell.dimred_methods import disentangler as mira_dt
from disentangler.frankencell.dimred_methods import disentangler_notune as notune
from disentangler.frankencell.dimred_methods import scvi
from disentangler.frankencell.dimred_methods import scvi_notune
from disentangler.frankencell.dimred_methods import harmony
from disentangler.frankencell.dimred_methods import scanorama

parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
subparsers = parser.add_subparsers(help = 'commands')

def add_subcommand(definition_file, cmd_name, **func_kwargs):

    subparser = subparsers.add_parser(cmd_name)
    definition_file.add_arguments(subparser)
    subparser.set_defaults(func = partial(definition_file.main, **func_kwargs))

add_subcommand(generate_KO_test, 'frankencell-gen-test')
add_subcommand(mira_dt, 'frankencell-mira', cost_beta = 2.)
add_subcommand(mira_dt, 'frankencell-mira-beta1', cost_beta = 1.)
add_subcommand(notune, 'frankencell-mira-notune')
add_subcommand(scvi, 'frankencell-scvi')
add_subcommand(harmony, 'frankencell-harmony')
add_subcommand(scanorama, 'frankencell-scanorama')

try:
    from disentangler.frankencell import eval_batch
    add_subcommand(eval_batch,'frankencell-eval-batch')
except ImportError:
    logging.warn('Could not load scib package, command frankencell-eval-batch will not be available')

for i in range(2, 7):
    add_subcommand(scvi_notune, 'frankencell-scvi-notune-' + str(i), n_latent = i)

def main():
    #____ Execute commands ___

    args = parser.parse_args()

    try:
        args.func #first try accessing the .func attribute, which is empty if user tries ">>>lisa". In this case, don't throw error, display help!
    except AttributeError:
        print(parser.print_help(), file = sys.stderr)
    else:
        args.func(args)
