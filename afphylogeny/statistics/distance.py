import itertools

import numpy as np
from skbio.stats.distance import DistanceMatrix

from afphylogeny.utils import sliding_window
from afphylogeny.statistics import d2


def create_distance_matrix(sequences):
    index_map = {seq: i for i, seq in enumerate(sequences)}

    data = np.zeros((len(sequences), len(sequences)))

    for seq1, seq2 in itertools.combinations(sequences, 2):
        if seq1 == seq2:
            continue

        i = index_map[seq1]
        j = index_map[seq2]

        data[i, j] = d2.d2_neighbourhood_dna(sliding_window(seq1, 8),
                                             sliding_window(seq2, 8))

    return DistanceMatrix(data, sequences.ids())
