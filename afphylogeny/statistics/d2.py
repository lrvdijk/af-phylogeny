import math
import logging
from collections import defaultdict, Counter


logger = logging.getLogger(__name__)


def d2(seq1_iter, seq2_iter):
    """Compute D2 score for given sequences `seq1` and `seq2`.

    Parameters
    ----------
    seq1_iter : iterator yielding each word
        This parameter should be an iterator returning each k-length word in
        the first sequence.
    seq2_iter : string
        This parameter should be an iterator returning each k-length word in
        the second sequence.

    Notes
    -----
    D2 is defined as:

    .. math:: D_2 = \sum_{w \in A^k}X_w Y_w

    Where:

    :math:`X_w`
        The number of occurences of word `w` in sequence 1

    :math:`Y_w`
        The number of occurences of word `w` in sequence 2

    :math:`A^k`
        The set of all possible words of length k with alphabet A
    """

    word_count1 = Counter(seq1_iter)
    word_count2 = Counter(seq2_iter)

    self_matching1 = sum(x**2 for x in word_count1.values())
    self_matching2 = sum(x**2 for x in word_count2.values())

    logger.debug('%s', str(word_count1))
    logger.debug('%s', str(word_count2))

    return abs(math.log(
        sum(
            word_count1[word] * word_count2[word] for word in word_count1
        ) / math.sqrt(self_matching1 * self_matching2)
    ))


class D2Normalised:
    def __init__(self, probs1, probs2):
        self.probs1 = probs1
        self.probs2 = probs2

    def __call__(self, seq1_iter, seq2_iter):
        word_count1 = defaultdict(int)
        word_count2 = defaultdict(int)

        for word in seq1_iter:
            word_count1[word] += 1

        for word in seq2_iter:
            word_count2[word] += 1

        n = len(word_count1)
        m = len(word_count2)

        score = 0.0
        self_matching1 = 0
        self_matching2 = 0

        for word in word_count1:
            normalised1 = word_count1[word] - (n * self.probs1[word])
            normalised2 = word_count2[word] - (m * self.probs2[word])

            score += (normalised1 * normalised2) / math.sqrt(
                normalised1**2 + normalised2**2
            )

            self_matching1 += (normalised1 * normalised1) / math.sqrt(
                normalised1**2 + normalised1**2
            )

        for word in word_count2:
            normalised2 = word_count2[word] - (m * self.probs2[word])

            self_matching2 += (normalised2 * normalised2) / math.sqrt(
                normalised2**2 + normalised2**2
            )

        return abs(math.log(
            score / math.sqrt(self_matching1 * self_matching2)
        ))


class D2Neighbourhood:
    def __init__(self, alphabet="ACTG"):
        self.alphabet = alphabet

    def _get_neighbourhood(self, word):
        logger.debug('Generating neighbourhood for %s', word)
        word = str(word)
        for i in range(len(word)):
            for letter in self.alphabet:
                yield "".join([word[:i], letter, word[i+1:]])

    def __call__(self, seq1_iter, seq2_iter):
        word_count1 = Counter(
            neighbour for word in seq1_iter for neighbour
            in self._get_neighbourhood(word)
        )

        word_count2 = Counter(
            neighbour for word in seq2_iter for neighbour
            in self._get_neighbourhood(word)
        )

        logger.debug('Word count 1: %s', str(word_count1))
        logger.debug('Word count 2: %s', str(word_count2))

        self_matching1 = sum(x**2 for x in word_count1.values())
        self_matching2 = sum(x**2 for x in word_count2.values())

        return abs(math.log(
            sum(
                word_count1[word] * word_count2[word] for word in word_count1
            ) / math.sqrt(self_matching1 * self_matching2)
        ))


d2_neighbourhood_dna = D2Neighbourhood()
d2_neighbourhood_rna = D2Neighbourhood("ACUG")


def distance(seq1, seq2, statistic=d2, k=8):
    return statistic(seq1.sliding_window(k), seq2.sliding_window(k))
