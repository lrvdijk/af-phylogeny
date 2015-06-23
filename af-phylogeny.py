import argparse
import logging

from skbio import SequenceCollection
from skbio.tree import nj

from afphylogeny.statistics.distance import create_distance_matrix
from afphylogeny.statistics import d2


def get_argument_parser():
    parser = argparse.ArgumentParser(
        description="Construct phylogenetic trees using an alignment free "
                    "method."
    )

    parser.add_argument(
        'infile', type=argparse.FileType('r'),
        help="The file containing the evolutionary related sequences."
    )

    parser.add_argument(
        'outfile', type=argparse.FileType('w'),
        help="The tree output filename."
    )

    parser.add_argument(
        '--format', '-f', default="fasta",
        help="The sequence collection input file format. See the scikit-bio "
             "documentation for more information about the supported types. "
             "Defaults to fasta."
    )

    parser.add_argument(
        '--target', '-t', default="newick",
        help="The tree output file format. See the scikit-bio "
             "documentation for more information about the supported types. "
             "Defaults to newick."
    )

    log_choices = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
    parser.add_argument(
        '--log-level', '-l', default="INFO", choices=log_choices,
        help="Set logging level. Default is info."
    )

    return parser

if __name__ == '__main__':
    parser = get_argument_parser()
    args = parser.parse_args()

    level = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.basicConfig(level=level)

    sequences = SequenceCollection.read(args.infile, format=args.format)
    dmatrix = create_distance_matrix(sequences, d2.distance,
                                     statistic=d2.d2_neighbourhood_dna)

    print(dmatrix)
    phylo_tree = nj(dmatrix)
    print(phylo_tree.ascii_art())
    phylo_tree.write(args.outfile, format=args.target)
