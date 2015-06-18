import os
import argparse
from io import StringIO

from skbio import SequenceCollection
from skbio.io import phylip


def get_argument_parser():
    parser = argparse.ArgumentParser(
        description="Split a concatenated PHYLIP file. Currently only PHYLIP "
                    "files in sequential format are supported."
    )

    parser.add_argument(
        'infile', type=argparse.FileType('r'),
        help="The PHYLIP formatted file to read from"
    )

    parser.add_argument('outdir', help="The output directory")

    parser.add_argument(
        '--format', '-f', default="fasta",
        help="The output format. See the scikit-bio documentation for more "
             "information (the skbio.io module)."
    )

    parser.add_argument(
        '--prefix', '-p', default="",
        help="Filename prefix. Each output filename will be prefixed with "
             "this string"
    )

    return parser


def concatenated_phylip_generator(fh):
    chunks = []

    for line in fh:
        if phylip._parse_phylip_header(line):
            # Valid PHYLIP header
            # It's either the start of the file or the start of a new
            # alignment block
            if chunks:
                yield "".join(chunks)

            chunks = [line]
        else:
            if line.strip():
                chunks.append(line)

    # Yield the last alignment block
    if chunks:
        yield "".join(chunks)


def convert_phylip(infile, outfile, format):
    seqs = SequenceCollection.read(
        infile, format='phylip',
        data_parser=phylip.relaxed_ids
    )

    seqs.write(outfile, format=format)


if __name__ == '__main__':
    parser = get_argument_parser()
    args = parser.parse_args()

    prefix = args.prefix or os.path.basename(
        args.infile.name[:args.infile.name.rfind('.')])

    counter = 0
    for phylip_chunk in concatenated_phylip_generator(args.infile):
        filename = os.path.join(args.outdir, "{}{}.{}".format(
            prefix, counter, args.format))

        with open(filename, 'w') as outfile:
            infile = StringIO(phylip_chunk)

            convert_phylip(infile, outfile, args.format)

        counter += 1
