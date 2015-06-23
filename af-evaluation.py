import os
import argparse
import logging
from pathlib import Path

from skbio import TreeNode
from skbio.io import FileFormatError


logger = logging.getLogger(__name__)


def get_argument_parser():
    parser = argparse.ArgumentParser(
        description="Compare sets of generated trees to eachother using the "
                    "Robinson-Foulds metric. This program compares all tree "
                    "files in one directory with matching tree files in "
                    "another directory."
    )

    parser.add_argument('dir1', help="First directory containing tree files")
    parser.add_argument('dir2', help="Second directory containing tree files")

    parser.add_argument(
        '--format', '-f', default="newick",
        help="The format of the tree files. Defaults to newick. See the "
             "scikit-bio toolkit for more information about the supported "
             "formats."
    )

    parser.add_argument(
        '--output-dir', '-o', default="",
        help="Output directory for evaluation files. Defaults to the first "
             "directory."
    )

    log_choices = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
    parser.add_argument(
        '--log-level', '-l', default="INFO", choices=log_choices,
        help="Set logging level. Default is info."
    )

    return parser


def get_tree_files(directory):
    current_path = Path(directory)

    return current_path.glob('*.tree')


if __name__ == '__main__':
    parser = get_argument_parser()
    args = parser.parse_args()
    output_dir = args.output_dir or args.dir1

    level = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.basicConfig(level=level)

    rf_scores = []

    for tree in get_tree_files(args.dir1):
        file1 = os.path.join(args.dir1, tree.name)
        file2 = os.path.join(args.dir2, tree.name)

        logger.debug('%s %s', file1, file2)

        if not os.path.exists(file2):
            logger.info('Skipping %s as the corresponding file %s does not '
                        'exists.', file1, file2)
            continue

        with open(file1) as fh1, open(file2) as fh2:
            try:
                tree1 = TreeNode.read(fh1, format=args.format)
                tree2 = TreeNode.read(fh2, format=args.format)
            except FileFormatError as e:
                logger.warning('Skipping evalutation of %s and %s due to '
                               'exception', file1, file2)
                logger.exception(e)
                continue

            rf = tree1.compare_rfd(tree2)

            outfile = os.path.join(
                output_dir, file1[:-4] + "evaluation")

            with open(outfile, 'w') as ofile:
                ofile.write(
                    "Tree 1\n"
                    "------\n\n"
                )
                ofile.write(file1 + "\n")
                ofile.write(tree1.ascii_art())
                ofile.write("\n\n")

                ofile.write(
                    "Tree 2\n"
                    "------\n\n"
                )
                ofile.write(file2 + "\n")
                ofile.write(tree2.ascii_art())
                ofile.write("\n\n")
                ofile.write("Robinson-Foulds: {}".format(rf))

            rf_scores.append(rf)
            logging.info('Robinson-Foulds score for tree %s and %s: %.2f',
                         file1, file2, rf)

    logging.info('Average Robinson-Foulds score: %.2f',
                 sum(rf_scores) / len(rf_scores))
