import numpy
from scipy.stats.stats import pearsonr

def distance_matrix_from_allele_counts(
    allele_dicts,
    min_reads=10,
    metric="pearson",
    alleles=None):
    """
    Return a distance matrix.
    """
    num_samples = len(allele_dicts[0])

    distance_matrix = numpy.zeros((num_samples, num_samples))
    for i in range(num_samples):
        for j in range(i, num_samples):
            vector_i = []
            vector_j = []
            for (row_num, row) in enumerate(allele_dicts):
                row_i = row[i]
                sum_i = float(sum(row_i.values()))
                row_j = row[j]
                sum_j = float(sum(row_j.values()))
                
                if min(sum_i, sum_j) < min_reads:
                    continue

                if alleles is None:
                    alleles_to_use = list(set(row_i).union(row_j))
                else:
                    alleles_to_use = alleles[row_num]
                vector_i.extend([row_i[a] / max(sum_i, 1) for a in alleles_to_use])
                vector_j.extend([row_j[a] / max(sum_j, 1) for a in alleles_to_use])
            assert len(vector_i) == len(vector_j)
            vector_i = numpy.array(vector_i)
            vector_j = numpy.array(vector_j)
            
            if metric == "pearson":
                (corr, p_value) = pearsonr(vector_i, vector_j)
                distance = numpy.sqrt(1 - corr)
            elif metric == "euclidean":
                distance = numpy.sqrt(((vector_i - vector_j)**2).sum())
            else:
                raise ValueError("Invalid metric: %s" % metric)

            distance_matrix[i, j] = distance_matrix[j, i] = distance

    return distance_matrix

