
import argparse
import sys
from disentangler.frankencell import generate_KO_test
from disentangler.frankencell.dimred_methods import disentangler as mira_dt
from disentangler.frankencell.dimred_methods import disentangler_notune as notune

parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
subparsers = parser.add_subparsers(help = 'commands')

def add_subcommand(definition_file, cmd_name):

    subparser = subparsers.add_parser(cmd_name)
    definition_file.add_arguments(subparser)
    subparser.set_defaults(func = definition_file.main)

add_subcommand(generate_KO_test, 'frankencell-gen-test')
add_subcommand(mira_dt, 'frankencell-mira')
add_subcommand(notune, 'frankencell-mira-notune')

def main():
    #____ Execute commands ___

    args = parser.parse_args()

    try:
        args.func #first try accessing the .func attribute, which is empty if user tries ">>>lisa". In this case, don't throw error, display help!
    except AttributeError:
        print(parser.print_help(), file = sys.stderr)
    else:
        args.func(args)
