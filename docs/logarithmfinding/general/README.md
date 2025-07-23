## Module: <code>general</code>
A module for finding general discrete logarithms.

## Import directive
```python
import quaspy.logarithmfinding.general
```

## Parent module
- [<code>logarithmfinding</code>](../README.md)

## Submodules
- [<code>postprocessing</code>](postprocessing/README.md)

  A module for solving frequency pairs yielded by the quantum part of Shor's algorithm for the general discrete logarithm problem [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084) and [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see Sect. 5.5), for the logarithm given the order.

- [<code>sampling</code>](sampling/README.md)

  A module for sampling a frequency pair (j, k) heuristically from the distribution induced by the quantum part of Shor's algorithm for finding a discrete logarithm d in a group of known order r, or from the distribution induced by the quantum part of Ekerå's variation thereof, depending on how parameters are selected.

