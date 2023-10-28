## Module: <code>algorithms</code>
A module for implementing Algorithms 1–4 from [[E22p]](https://doi.org/10.48550/arXiv.2201.07791).

## Import directive
```python
import quaspy.orderfinding.general.postprocessing.ekera.internal.algorithms
```

## Parent module
- [<code>internal</code>](../README.md)

## Functions
- [<code>algorithm1(g, r_tilde, m, ..)</code>](algorithm1.md)

  Recovers a multiple rp of r, assuming r_tilde is such that r = d * r_tilde where d is cm-smooth.

- [<code>algorithm2(g, r_tilde, m, ..)</code>](algorithm2.md)

  Recovers r, assuming r_tilde is such that r = d * r_tilde where d is cm-smooth.

- [<code>algorithm3(g, r_tilde, m, ..)</code>](algorithm3.md)

  Recovers r, assuming r_tilde is such that r = d * r_tilde where d is cm-smooth.

- [<code>algorithm4(g, S, m, ..)</code>](algorithm4.md)

  Returns a subset Sp of S consisting of all r_tilde in S that are such that d * r_tilde is a positive integer multiple of r, where d is cm-smooth.

- [<code>is_valid_r_tilde(r_tilde, m)</code>](is_valid_r_tilde.md)

  Checks if r_tilde is an integer on [1, 2^m).

