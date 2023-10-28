## Function: <code>solve\_r\_for\_factors(r, N, ..)</code>
Attempts to factor N completely given the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N, where g need not be explicitly specified.

This by using the algorithm in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).

## Import directive
```python
from quaspy.factoring.general.postprocessing.ekera import solve_r_for_factors
```

## Parent module
- [<code>ekera</code>](README.md)

## Prototype
```python
def solve_r_for_factors(r,
                        N,
                        c = 1,
                        k = None,
                        timeout = None,
                        verbose = False,
                        opt_split_factors_with_multiplicity = True,
                        opt_report_accidental_factors = True,
                        opt_abort_early = True,
                        opt_square = True,
                        opt_exclude_one = True,
                        opt_process_composite_factors = OptProcessCompositeFactors.SEPARATELY_MOD_Np)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| r | The order r of g, or a positive integer multiple of r. |
| N | The integer N. |
| c | A parameter c >= 1 that specifies the maximum size of the missing cm-smooth component in lambda'(N) when solving r for the complete factorization of N. In this context, m is the bit length of N, and cm-smoothness is defined as in [E21b, E22p].<br><br>As is explained in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1), increasing c increases the success probability, at the expense of increasing the runtime. |
| k | The maximum number of iterations to perform. Defaults to None. If k is set to None, as many iterations as are necessary to completely factor N will be performed. If k is explicitly specified, and the complete factorization of N has not been found after k iterations, an exception of type IncompleteFactorizationException will be raised. |
| timeout | A timeout in seconds. Defaults to None. If the timeout is set to None, as much time as is necessary to completely factor N will be used. If a timeout is explicitly specified, and the complete factorization of N has not been found when the timeout elapses, an exception of type IncompleteFactorizationException will be raised. |
| verbose | A flag that may be set to True to print intermediary results. Defaults to False.<br><br>All other parameters control various optimizations. They are documented below, and in the code. It is recommended to use the defaults.<br><br> |
| opt_split_factors_with_multiplicity | A flag that may either be set to True (default option) or to False.<br><br>When set to True, as is the default, the solver initially tests if gcd(r, N) yields a non-trivial factor of N. If so, the non-trivial factor is reported.<br><br>When set to False, the aforementioned test is not performed.<br><br>To understand the test, note that if p^e divides N, for e > 1 an integer and p a large prime, then p^(e-1) is likely to also divide r. Note furthermore that it is relatively inexpensive to split N by computing gcd(r, N), compared to exponentiating modulo N.<br><br>For random N void of small factors, it is very unlikely for prime factors to occur with multiplicity. For such N, this optimization is hence not of any practical use. It is only if N is likely to have factors that occur with multiplicity, for some special reason, that this optimization is useful in practice. This optimization may also yield non-trivial factors if N has small factors. However, such factors would typically first be removed, e.g. by trial division, before calling upon these more elaborate factoring techniques.<br><br>Note that for N as in the problem instances setup in Appendix A.3 of [[E21b]](https://doi.org/10.1007/s11128-021-03069-1), prime factors are intentionally forced to occur with multiplicity with high probability (when e_max > 1). This so as to test that such special cases are handled correctly by the solver. For such N, this optimization is likely to report non-trivial factors. This served as our rationale for including it as an option.<br><br>This optimization is described in earlier works in the literature, see for instance [GLMS15] for one such description.<br><br>[GLMS15] Grosshans, F., T. Lawson, F. Morain and B. Smith: "Factoring Safe Semiprimes with a Single Quantum Query". ArXiv 1511.04385 (2015).<br><br> |
| opt_report_accidental_factors | A flag that may either be set to True (default option) or to False.<br><br>When set to True, as is the default, the solver reports non-trivial factors of N found "by accident" when sampling x (denoted x_j in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1)) uniformly at random from Z_{N'}^*, for N' equal to N, or to the product of all pairwise coprime composite factors of N in the factor collection, depending on which option is selected for opt_process_composite_factors. When set to False, such non-trivial factors are not reported.<br><br>Note that it is very unlikely for non-trivial factors to be found by accident if N is void of small prime factors. It is only if small factors of N are not first removed, e.g. by trial division, that factors are likely to be found by accident. On the other hand, if we do find factors by accident, it is only logical to report them. This served as our rationale for including this optimization as an option.<br><br>Note furthermore that if N has small factors when this optimization is enabled, and if opt_process_composite_factors is furthermore set to OptProcessCompositeFactors.JOINTLY_MOD_N, then the same non-trivial factor of N may be repeatedly found by accident when sampling. This may generate long printouts.<br><br>Also note that factors that are found "by accident" when sampling g uniformly at random from Z_N^* are not reported even if opt_report_accidental_factors is set to True. This is because such factors, if found in practice, would typically affect for which N order finding is performed in the first place.<br><br> |
| opt_abort_early | A flag that may either be set to True (default option) or to False.<br><br>When set to True, as is the default, the solver computes x^(2^i o) for 0, 1, .., min(s, t), for t as in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) and s the least non-negative integer such that x^(2^s o) = 1. When set to False, the solver instead computes x^(2^i o) for i = 0, 1, .., t.<br><br>Note that when opt_square (see below) is set to True, the solver first computes x^o. It then takes consecutive squares to form x^(2^i o) for each i. It follows that the saving incurred by this optimization is fairly minor when the opt_square flag is set to True, as it is trivial to square one. It is only if the opt_square flag is set to False that the saving is substantial.<br><br> |
| opt_square | A flag that may either be set to True (default option) or to False.<br><br>When set to True, as is the default, the solver first computes x^o. It then takes consecutive squares to form x^(2^i o) for each i. When set to False, the solver naïvely computes x^(2^i o) from scratch for each i.<br><br>This optimization is described in Sect. 3.2.1 of [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).<br><br> |
| opt_exclude_one | A flag that may be set either to True (default option) or to False.<br><br>When set to False, the solver selects g and x (denoted x_j in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1)) uniformly at random from Z_N^* and Z_{N'}^*, respectively, for N' equal to N, or to the product of all pairwise coprime composite factors of N in the factor collection, depending on which option is selected for the opt_process_composite_factors flag.<br><br>When set to True, as is the default, the solver excludes one by instead selecting g and x uniformly at random from Z_N^* \ {1} and Z_{N'}^* \ {1}, respectively.<br><br>This optimization is described in Sect. 3.2.1 of [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).<br><br> |
| opt_process_composite_factors | An enumeration entry from the OptProcessCompositeFactors class that specifies how x should be selected, and composite factors processed, by the solver.<br><br>As is described in OptProcessCompositeFactors, there are three options:<br><br>1. OptProcessCompositeFactors.JOINTLY_MOD_N<br><br>The solver selects x (denoted x_j in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1)) uniformly at random from Z_N^*. It then exponentiates x modulo N. This is how the unoptimized algorithm is described in Sect. 3.2 of [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).<br><br>2. OptProcessCompositeFactors.JOINTLY_MOD_Np<br><br>The solver selects x uniformly at random from Z_{N'}^*. It then exponentiates x modulo N'.<br><br>Above N' is the product of all pairwise coprime composite factors of N currently stored in the factor collection.<br><br>3. OptProcessCompositeFactors.SEPARATELY_MOD_Np (default option)<br><br>The solver selects x uniformly at random from Z_{N'}^*, for N' the product of all pairwise coprime composite factors of N currently stored in the factor collection.<br><br>The solver then exponentiates x modulo N' where N' now runs over all pairwise coprime composite factors of N currently stored in the factor collection. This is the default option.<br><br>Note that all three options are equivalent with respect to their ability to find non-trivial factors of N. The options differ only in terms of their arithmetic complexity. Although several exponentiations may be required when the default option is used, this fact is, in general, more than compensated for by the fact that the moduli are smaller, leading the default option to outperform the other two options.<br><br>This optimization is described in Sect. 3.2.1 of [[E21b]](https://doi.org/10.1007/s11128-021-03069-1). |

## Return value
A set of all distinct prime factors that divide N.
