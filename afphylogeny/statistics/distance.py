import itertools
import logging
import multiprocessing as mp

import numpy as np
from skbio.stats.distance import DistanceMatrix


logger = logging.getLogger(__name__)


def calculate_distance(seq1, seq2, distance_func, args, kwargs):
    logger.info('Calculating distances between %s and %s',
                seq1.metadata['id'], seq2.metadata['id'])

    dist = distance_func(seq1, seq2, *args, **kwargs)

    logger.debug('Distance between %s and %s is %.2f', seq1.metadata['id'],
                 seq2.metadata['id'], dist)

    return (seq1.metadata['id'], seq2.metadata['id'], dist)


def calculate_distance_kwargs(kwargs):
    return calculate_distance(**kwargs)


def create_distance_matrix(sequences, distance_func, pool=1, *args, **kwargs):
    index_map = {seq.metadata['id']: i for i, seq in enumerate(sequences)}

    data = np.zeros((len(sequences), len(sequences)))

    if pool > 1:
        with mp.Pool(pool) as p:
            async_iter = p.imap_unordered(
                calculate_distance_kwargs, (
                    {
                        'seq1': seq1,
                        'seq2': seq2,
                        'distance_func': distance_func,
                        'args': args,
                        'kwargs': kwargs
                    } for seq1, seq2 in itertools.combinations(sequences, 2)
                ),
                chunksize=32,
            )

            for seq1_id, seq2_id, dist in async_iter:
                i = index_map[seq1_id]
                j = index_map[seq2_id]

                data[i, j] = dist
                data[j, i] = dist
    else:
        for seq1, seq2 in itertools.combinations(sequences, 2):
            seq1_id, seq2_id, dist = calculate_distance(
                seq1, seq2, distance_func, args, kwargs)

            i = index_map[seq1_id]
            j = index_map[seq2_id]

            data[i, j] = dist
            data[j, i] = dist

    return DistanceMatrix(data, sequences.ids())
