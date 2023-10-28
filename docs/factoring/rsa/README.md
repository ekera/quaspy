## Module: <code>rsa</code>
A module for factoring large random RSA integers.

This module uses Ekerå–Håstad's algorithm to factor RSA integers, as described in [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20), with improvements from [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754).

## Import directive
```python
import quaspy.factoring.rsa
```

## Parent module
- [<code>factoring</code>](../README.md)

## Submodules
- [<code>postprocessing</code>](postprocessing/README.md)

  A module for splitting N = pq into the large l-bit prime factors p and q given d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1).

## Functions
- [<code>setup_d_given_p_q(p, q)</code>](setup_d_given_p_q.md)

  Sets up d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1) given p and q.

- [<code>setup_x_given_g_N(g, N)</code>](setup_x_given_g_N.md)

  Sets up x = g^d' for d' = (N - 1) / 2 - 2^(l - 1) given g and N, for N the product of two large random distinct l-bit primes.

