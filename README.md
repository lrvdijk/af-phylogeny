Alignment-Free Phylogeny
========================

This is a Python package which provides several scripts for reconstructing 
phylogenetic trees from a collection of DNA sequences without using a sequence 
alignment method.

Currently, this module creates the phylogenetic tree based on the D2 
neighbourhood statistic [(1)][d2]. 

Dependencies
------------

* Python 3.4+
* NumPy
* SciPy
* Pandas
* Custom fork of [scikit-bio][].

This program has only been tested under Linux.

Installation
------------

Because we currently use some custom additions to [scikit-bio][] which are not 
incorporated in the master branch yet, we recommend to create a custom Python 
virtual environment.

    $ pyvenv af-phylogeny
    $ cd af-phylogeny
    $ source bin/activate

Now you're inside your activated Python virtual environment. Clone the 
repository and  install the dependencies for this package using Pip, this also 
automatically downloads the required [scikit-bio][] branch.

    $ git clone https://github.com/sh4wn/af-phylogeny
    $ cd af-phylogeny
    $ pip install -r requirements.txt

Usage
-----

### `phylip-splitter.py` - Split concatenated PHYLIP files

This script can split a PHYLIP file with multiple sets of sequences into 
separate files again. The script allows you to specify the output directory and 
the output format. The default output format is FASTA. 

Note that this script currently only supports a PHYLIP file in sequential 
format, which means that the whole sequence needs to be on the same line.

    usage: phylip-splitter.py [-h] [--format FORMAT] [--prefix PREFIX]
                              infile outdir

    Split a concatenated PHYLIP file. Currently only PHYLIP files in sequential
    format are supported.

    positional arguments:
      infile                The PHYLIP formatted file to read from
      outdir                The output directory

    optional arguments:
      -h, --help            show this help message and exit
      --format FORMAT, -f FORMAT
                            The output format. See the scikit-bio documentation
                            for more information (the skbio.io module).
      --prefix PREFIX, -p PREFIX
                            Filename prefix. Each output filename will be prefixed
                            with this string

### `af-phylogeny.py` - Construct phylogenetic tree

This is the actual work horse of this Python package. This script creates a 
phylogenetic tree from a collection of sequences using an alignment free 
method. 

It supports several input and output file formats, please refer to the 
[skbio.io module documentation][skbio-io] for more information.

This script is also able to use multiple cores. It tries to be smart about when 
to enable this: when there are more than 16 sequences in the file, it enables 
multicore support and will automatically spawn a number of child processes 
corresponding to the number of cores. This behaviour can be overwritten by 
setting the number of child processes manually.

    usage: af-phylogeny.py [-h] [--format FORMAT] [--target FORMAT] [--parallel 
    N]
                           [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                           infile outfile

    Construct phylogenetic trees using an alignment free method.

    positional arguments:
      infile                The file containing the evolutionary related
                            sequences.
      outfile               The tree output filename.

    optional arguments:
      -h, --help            show this help message and exit
      --format FORMAT, -f FORMAT
                            The sequence collection input file format. See the
                            scikit-bio documentation for more information about
                            the supported types. Defaults to fasta.
      --target FORMAT, -t FORMAT
                            The tree output file format. See the scikit-bio
                            documentation for more information about the supported
                            types. Defaults to newick.
      --parallel N, -p N    Enable multicore support. This option specifies the
                            number of child processes to spawn. By default it
                            automatically detects the number of cpu cores and
                            enables multicore support when there are more than 16
                            sequences in the file to analyze.
      --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}, -l {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                            Set logging level. Default is info.

### `af-evaluation.py` - Evaluate the generated trees

This script is useful when you want to compare two sets of tree files to each 
other.

It expects the tree files generated with *method1* and *method2* in separate 
directories, with corresponding filenames. It compares tree with the same 
filename using the Robinson-Foulds metric. For each tree it generates an 
*evaluation file*, containing ASCII art of both trees and the Robinson-Foulds 
distance. By default it stores this file in the first directory specified, but 
you can override the output directory with the corresponding command line 
option.

    usage: af-evaluation.py [-h] [--format FORMAT] [--output-dir OUTPUT_DIR]
                            [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                            dir1 dir2

    Compare sets of generated trees to eachother using the Robinson-Foulds metric.
    This program compares all tree files in one directory with matching tree files
    in another directory.

    positional arguments:
      dir1                  First directory containing tree files
      dir2                  Second directory containing tree files

    optional arguments:
      -h, --help            show this help message and exit
      --format FORMAT, -f FORMAT
                            The format of the tree files. Defaults to newick. See
                            the scikit-bio toolkit for more information about the
                            supported formats.
      --output-dir OUTPUT_DIR, -o OUTPUT_DIR
                            Output directory for evaluation files. Defaults to the
                            first directory.
      --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}, -l {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                            Set logging level. Default is info.

### `mafft-tree-fixer.py` - Fix tree files generated by MAFFT

For our internal evaluation we compared the trees generated by our alignment 
free method to tree files generated with [MAFFT-LINSI][mafft]. MAFFT is able to 
output a tree file, but scikit-bio is quite strict in its Newick formatting, 
and is unable to read it because the generated tree by MAFFT misses a closing 
semicolon. Furthermore, MAFFT also modifies the sequence names, and to be able 
to compare two trees the node names must be the same.

This script fixes all of that: it restores the original sequence names, and 
outputs a properly formatted Newick file.

    usage: mafft-tree-fixer.py [-h] infile [outfile]

    Fix the tree files generated by MAFFT. Tree files generated by MAFFT do not
    contain a closing semi-colon, and therefore scikit-bio is unable to read them.
    They also modify the name of the sequences, so we revert that change to be
    able to compare the trees using Robinson-Foulds.

    positional arguments:
      infile      The tree file to fix
      outfile     The output file. If not specified file will be modified inplace.

    optional arguments:
      -h, --help  show this help message and exit


Workflow
--------

This repository already contains a set of sequences: several synthetic trees 
and one real tree. To automate most tasks a Makefile is included. Before we 
explain how to use this Makefile let us explain the directory structure:

- af-phylogeny/ - Container folder
    - afphylogeny/ - Python package source code
    - data/ - Directory holding all data files (sequences, trees)
        - af/ - Output directory for tree files created with out alignment free 
          method.
        - mafft/ - Output directory for tree files created with MAFFT-LINSI.
        - sequences/ - The folder containing all sequence collections.

Using the Makefile you can easily run MAFFT or our own alignment free method on 
all the sequence collections. 

### `make mafft` - Create tree files using MAFFT

For each FASTA file with multiple sequences in the data/sequences/ folder it 
runs MAFFT-LINSI on that file, and outputs a Newick formatted tree file in 
data/mafft/. It also automatically runs the `mafft-tree-fixer.py` script, so 
you don't need to worry about that.

### `make af` - Create tree files using our Alignment Free method

For each FASTA file with multiple sequences in the data/sequences/ folder it 
runs the `af-phylogeny.py` script, and the tree is stored in data/af/.

### `make all` - Create all trees

Runs both `make mafft` and `make af`.

### Multicore support

Make is also able to run multiple builds in parallel. This is useful when you 
have a lot of small trees which are quickly generated (and for which 
`af-phylogeny.py` does **not** enable multicore support). 

To enable multicore support in Make, use the `-j` command line option:

    make -j4 all

This runs four builds in parallel. 

[scikit-bio]: http://scikit-bio.org
[skbio-io]: http://scikit-bio.org/docs/latest/io.html#supported-file-formats
[mafft]: http://mafft.cbrc.jp/alignment/software/
