# Variable to be used for substiturion
empty:=

# Collect all the sequence files
SEQUENCES = $(subst data/sequences/,$(empty),$(wildcard data/sequences/*.fasta))

# Tree files generated with MAFFT
MAFFT_TREES = $(addprefix data/mafft/, $(addsuffix .tree, $(SEQUENCES)))

# Tree files generated with out own alignment-free method
AF_TREES = $(addprefix data/af/, $(addsuffix .tree, $(SEQUENCES)))

all: trees

trees: $(MAFFT_TREES) $(AF_TREES)
	
data/mafft/%.fasta.tree: data/sequences/%.fasta
	mafft-linsi --treeout --thread -1 $< > $@

data/af/%.fasta.tree: data/sequences/%.fasta
	@echo 'python af-tree-construction.py $< $@'

print-%:
	@echo '$*=$($*)'
	@echo '  origin = $(origin $*)'
	@echo '  flavor = $(flavor $*)'
	@echo '   value = $(value  $*)'

clean:
	rm -f $(MAFFT_TREES)
	rm -f $(AF_TREES)

.PHONY: all distances trees clean
