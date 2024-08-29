## Module: <code>ekera</code>
A module for solving a frequency j yielded by the quantum part of Shor's order-finding algorithm for the order r, or optionally for a positive integer multiple of r. This by using the post-processing algorithms in [[E24]](https://doi.org/10.1145/3655026).

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

  Attempts to compute the order r of g, or a positive integer multiple thereof, given a frequency j yielded by the quantum part of Shor's order-finding algorithm, by using the post-processing algorithms described in detail in [[E24]](https://doi.org/10.1145/3655026).

- [<code>solve_j_for_r_mod_N(j, m, l, g, N, ..)</code>](solve_j_for_r_mod_N.md)

  Attempts to compute the order r of g mod N, or a positive integer multiple thereof, given a frequency j yielded by the quantum part of Shor's order-finding algorithm, by using the post-processing algorithms described in detail in [[E24]](https://doi.org/10.1145/3655026).

## Classes
- [<code>SolutionMethods</code>](SolutionMethods.md)

  An enumeration of methods for solving j for r.
