import math
from collections import defaultdict

from afphylogeny.utils import sliding_window


def d2(seq1, seq2):
    """Compute D2 score for given sequences `seq1` and `seq2`.

    Parameters
    ----------
    seq1_iter : iterator yielding each word
        This parameter should be an iterator returning each k-length word in
        the first sequence.
    seq2 : string
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

    word_count1 = defaultdict(int)
    word_count2 = defaultdict(int)

    for word in sliding_window(seq1):
        word_count1[word] += 1

    for word in sliding_window(seq2):
        word_count2[word] += 1

    return sum([word_count1[word] * word_count2[word] for word in word_count1])


def d2_normalised(seq1_iter, seq2_iter, probs1, probs2):
    word_count1 = defaultdict(int)
    word_count2 = defaultdict(int)

    for word in seq1_iter:
        word_count1[word] += 1

    for word in seq2_iter:
        word_count2[word] += 1

    n = len(word_count1)
    m = len(word_count2)

    score = 0.0

    for word in word_count1:
        normalised1 = word_count1[word] - (n * probs1[word])
        normalised2 = word_count2[word] - (m * probs2[word])

        score += (normalised1 * normalised2) / math.sqrt(
            normalised1**2 + normalised2**2
        )

    return score


def _get_neighbourhood(word, alphabet):
    for i in range(len(word)):
        for letter in alphabet:
            yield "".join(word[:i], letter, word[i+1:])


def d2_neighbourhood(seq1_iter, seq2_iter, alphabet="ACTG"):
    word_count1 = defaultdict(int)
    word_count2 = defaultdict(int)

    for word in seq1_iter:
        for neighbour in _get_neighbourhood(word, alphabet):
            word_count1[neighbour] += 1

    for word in seq2_iter:
        for neighbour in _get_neighbourhood(word, alphabet):
            word_count2[neighbour] += 1

    return sum([word_count1[word] * word_count2[word] for word in word_count1])
