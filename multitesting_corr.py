""" Implement methods that adjust p-values to account for multiple comparisons.
"""

import numpy as np

def p_multitest(pvals, method='holm', alpha=0.05):
    """Remedy multiple comparison problem by controlling
    familywise error rate to be LEQ <alpha>.
    
    Parameters
    ----------
    pvals: np.array of floats, p-values to adjust
    method: str, method to use
    alpha: float, desired probability of Type-I error
    
    Returns
    ----------
    p_map: dict of float, holds original and adjusted p-values
    """
    supported_methods = ('holm', 'bonferroni')
    
    if method == 'holm':
        return holm_correct(pvals, alpha)
    elif method == 'bonferroni':
        return bonferroni_correct(pvals)
    else:
        raise ValueError(f"Unknown method {method}, must be one of {','.join(supported_methods)}")

def holm_correct(pvals, alpha=0.05):
    """Implementation of Holm-Bonferroni correction.
    References:
    [1] http://www.gs.washington.edu/academics/courses/akey/56008/lecture/lecture10.pdf
    """
    pvals_sorted = np.sort(pvals)
    pvals_adjust = np.empty(shape=len(pvals))
    m = len(pvals)

    for k in range(m):
        alpha_adjusted = alpha / (m + 1 - (k + 1))
        pvals_adjust[k] = min(pvals_sorted[k] * (m + 1 - (k + 1)), 1)
    
    # ensure that highest adjusted p-value is chosen for duplcate p-values
    if len(pvals_adjust) > 1: 
        i = 1
        for a, b in zip(pvals_adjust[1:], pvals_adjust):
            pvals_adjust[i] = max(a, b)
            i += 1    
    
    p_map = dict()
    for p in range(m):
        p_map[pvals_sorted[p]] = pvals_adjust[p]
                
    return p_map

def bonferroni_correct(pvals):
    """Assumes equal alocation of test size alpha.
    [1] http://www.gs.washington.edu/academics/courses/akey/56008/lecture/lecture10.pdf
    """
    n = len(pvals)
    pvals_adjust = {pvalue_old: min(pvalue_old * n, 1) for pvalue_old in pvals}
    return pvals_adjust
