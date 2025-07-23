## Module: <code>sampling</code>
A module for heuristically sampling frequency pairs from the distribution induced by the quantum part of Shor's algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for finding discrete logarithms in groups of known order, modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

This module may also be used to sample frequency pairs from the distribution induced by the quantum parts of Ekerå–Håstad's [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) and Ekerå's [[E21]](https://doi.org/10.1515/jmc-2020-0006) variations of said algorithm, depending on how parameters are selected.

## Import directive
```python
import quaspy.logarithmfinding.sampling
```

## Parent module
- [<code>logarithmfinding</code>](../README.md)

## Functions
- [<code>sample_j_k_given_d_r_heuristic(d, r, m, sigma, l, ..)</code>](sample_j_k_given_d_r_heuristic.md)

  Samples a frequency pair (j, k) heuristically from the distribution induced by the quantum part of Shor's algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for finding a discrete logarithm d in a group of known order r, modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

