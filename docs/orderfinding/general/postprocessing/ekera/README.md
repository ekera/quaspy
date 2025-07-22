## Module: <code>ekera</code>
A module for solving a frequency j yielded by the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for the order r, or optionally for a positive integer multiple of r, by using the classical post-processing algorithms from [[E24]](https://doi.org/10.1145/3655026).

This module furthermore contains functions for solving a list of n frequencies [j_1, ..., j_n] yielded by n runs of the quantum part of Seifert's variation [[Seifert01]](https://doi.org/10.1007/3-540-45353-9_24) of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for r, or for a positive integer multiple of r. This by using the classical post-processing algorithms from [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see Sect. 5.4) and [[E21]](https://doi.org/10.1515/jmc-2020-0006), with supporting functions from [[E24]](https://doi.org/10.1145/3655026).

Throughout this module, the algorithms are as described in [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf), [[E24]](https://doi.org/10.1145/3655026) and [[E21]](https://doi.org/10.1515/jmc-2020-0006). The notation is also inherited from said works.

## Import directive
```python
import quaspy.orderfinding.general.postprocessing.ekera
```

## Parent module
- [<code>postprocessing</code>](../README.md)

## Submodules
- [<code>internal</code>](internal/README.md)

  A module for internal classes and functions used by the functions in the parent module for solving a frequency j yielded by the quantum part of Shor's order-finding algorithm for the order r.

## Functions
- [<code>solve_j_for_r(j, m, l, g, ..)</code>](solve_j_for_r.md)

  Attempts to compute the order r of g, or a positive integer multiple thereof, given a frequency j yielded by the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), by using the post-processing algorithms described in detail in [[E24]](https://doi.org/10.1145/3655026).

- [<code>solve_j_for_r_mod_N(j, m, l, g, N, ..)</code>](solve_j_for_r_mod_N.md)

  Attempts to compute the order r of g mod N, or a positive integer multiple thereof, given a frequency j yielded by the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), by using the post-processing algorithms described in detail in [[E24]](https://doi.org/10.1145/3655026).

- [<code>solve_multiple_j_for_r(j_list, m, l, g, ..)</code>](solve_multiple_j_for_r.md)

  Attempts to compute the order r of g, or a positive integer multiple thereof, given a list of n frequencies [j_1, ..., j_n] yielded by n independent runs of the quantum part of Seifert's variation [[Seifert01]](https://doi.org/10.1007/3-540-45353-9_24) of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) as described in [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see Sect. 5.4) and [[E21]](https://doi.org/10.1515/jmc-2020-0006) (see App. A). This by using Ekerå's lattice-based classical post-processing from [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) and [[E21]](https://doi.org/10.1515/jmc-2020-0006), with supporting functions from [[E24]](https://doi.org/10.1145/3655026).

- [<code>solve_multiple_j_for_r_mod_N(j_list, m, l, g, N, ..)</code>](solve_multiple_j_for_r_mod_N.md)

  Attempts to compute the order r of g mod N, or a positive integer multiple thereof, given a list of n frequencies [j_1, ..., j_n] yielded by n independent runs of the quantum part of Seifert's variation [[Seifert01]](https://doi.org/10.1007/3-540-45353-9_24) of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) as described in [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see Sect. 5.4) and [[E21]](https://doi.org/10.1515/jmc-2020-0006) (see App. A). This by using Ekerå's lattice-based classical post-processing from [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) and [[E21]](https://doi.org/10.1515/jmc-2020-0006), with supporting functions from [[E24]](https://doi.org/10.1145/3655026).

## Classes
- [<code>EnumerationOptions</code>](EnumerationOptions.md)

  An enumeration of options for solving a list of n frequencies [j_1, ..., j_n] for r.
- [<code>SolutionMethods</code>](SolutionMethods.md)

  An enumeration of methods for solving j for r.
