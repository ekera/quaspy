## Class: <code>OptProcessCompositeFactors</code>
An enumeration of optimization options for how the solver is to select x, and to process composite factors.

There are three optimization options:

1. JOINTLY_MOD_N

Select x uniformly at random from Z_N^*, for N the number to be factored, and exponentiate x modulo N to 2^t o.

This is as described in the algorithm in Sect. 3.2 of [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).

2. JOINTLY_MOD_Np

Select x uniformly at random from Z_N'^*, for N' the product of all pairwise coprime composite factors of N currently stored in the collection, and exponentiate x modulo N' to 2^t o.

This is as above, but with optimizations from Sect. 3.2.1 of [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).

3. SEPARATELY_MOD_Np

Select x uniformly at random from Z_N'^*, for N' the product of all pairwise coprime composite factors of N currently stored in the collection. Exponentiate x modulo N' to 2^t o, as N' runs over the pairwise coprime composite factors of N currently stored in the collection.

This is as above, with more optimizations from Sect. 3.2.1 of [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).

Note that all three options are equivalent with respect to their ability to find non-trivial factors of N. The options differ only in terms of their arithmetic complexity. Although several exponentiations may be required when the default option is used, this fact is, in general, more than compensated for by the fact that the moduli are smaller, leading the default option to outperform the other two options.

For further details, see Sect. 3.2.1 of [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).

## Import directive
```python
from quaspy.factoring.general.postprocessing.ekera import OptProcessCompositeFactors
```

## Parent module
- [<code>ekera</code>](README.md)

## Members
- <code>JOINTLY_MOD_N</code>
- <code>JOINTLY_MOD_Np</code>
- <code>SEPARATELY_MOD_Np</code>

