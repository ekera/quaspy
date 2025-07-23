## Module: <code>sampling</code>
A module for sampling a frequency j from the distribution induced by the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for a given order r. This by using the sampling procedure described in [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see in particular Sects. 5.3.5 and 5.4.3).

This module can also be used to sample from the distribution induced by the quantum part of Seifert's variation [[Seifert01]](https://doi.org/10.1007/3-540-45353-9_24) of Shor's order-finding algorithm for a given order r, depending on how parameters are selected.

Throughout this module, the algorithms are as described in [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf), [[E24]](https://doi.org/10.1145/3655026) and [[E21]](https://doi.org/10.1515/jmc-2020-0006). The notation is also inherited from said works.

## Import directive
```python
import quaspy.orderfinding.general.sampling
```

## Parent module
- [<code>general</code>](../README.md)

## Functions
- [<code>optimal_j_for_z_r(z, r, m, l)</code>](optimal_j_for_z_r.md)

  Computes and returns the optimal frequency j0(z) for z in [0, r).

- [<code>sample_j_given_r(r, m, l, ..)</code>](sample_j_given_r.md)

  Samples a frequency j from the distribution induced by the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for a given order r. This by using the sampling procedure described in [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see in particular Sects. 5.3.5 and 5.4.3).

