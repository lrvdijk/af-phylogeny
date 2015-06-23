import itertools
import logging

import numpy as np
from skbio.stats.distance import DistanceMatrix


logger = logging.getLogger(__name__)


def create_distance_matrix(sequences, distance_func, *args, **kwargs):
    index_map = {seq.metadata['id']: i for i, seq in enumerate(sequences)}

    data = np.zeros((len(sequences), len(sequences)))

    for seq1, seq2 in itertools.combinations(sequences, 2):
        logger.info('Calculating distances between %s and %s',
                    seq1.metadata['id'], seq2.metadata['id'])

        if seq1 == seq2:
            continue

        i = index_map[seq1.metadata['id']]
        j = index_map[seq2.metadata['id']]

        dist = distance_func(seq1, seq2, *args, **kwargs)

        logger.debug('Distance between %s and %s is %.2f', seq1.metadata['id'],
                     seq2.metadata['id'], dist)
        data[i, j] = dist
        data[j, i] = dist

    return DistanceMatrix(data, sequences.ids())
