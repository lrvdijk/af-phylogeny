# Variable to be used for substitution
empty:=

# Collect all the sequence files
SEQUENCES = $(subst data/sequences/,$(empty),$(wildcard data/sequences/*.fasta))

# Tree files generated with MAFFT
MAFFT_TREES = $(addprefix data/mafft/, $(addsuffix .tree, $(SEQUENCES)))

# Tree files generated with out own alignment-free method
AF_TREES = $(addprefix data/af/, $(addsuffix .tree, $(SEQUENCES)))

all: trees

mafft: $(MAFFT_TREES)
af: $(AF_TREES)
trees: mafft af
	
data/mafft/%.fasta.tree: data/sequences/%.fasta
	mafft-linsi --treeout --thread -1 $< > $@
	mv $(addsuffix .tree, $<) $@

data/af/%.fasta.tree: data/sequences/%.fasta
	python af-phylogeny.py -l WARNING $< $@

print-%:
	@echo '$*=$($*)'
	@echo '  origin = $(origin $*)'
	@echo '  flavor = $(flavor $*)'
	@echo '   value = $(value  $*)'

clean:
	rm -f $(MAFFT_TREES)
	rm -f $(AF_TREES)

.PHONY: all mafft af trees clean
