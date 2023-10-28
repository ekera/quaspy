## Module: <code>sampling</code>
A module for exactly or heuristically sampling an element g uniformly at random from the multiplicative group of the ring of integers modulo N and returning the order r of g. This when the complete factorization of N is known.

The procedures in this module are described in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1), in the [factoritall repository](https://www.github.com/ekera/factoritall) (available at https://www.github.com/ekera/factoritall), and in the short note M. Ekerå: "A note on sampling random elements of known order from $Z^*_N$ when the factors of $N$ are known" which is to appear.

## Import directive
```python
import quaspy.factoring.sampling
```

## Parent module
- [<code>factoring</code>](../README.md)

## Functions
- [<code>sample_g_r_given_N(N, N_factors, ..)</code>](sample_g_r_given_N.md)

  Returns [g, r], for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N, and r either a heuristic estimate of the order of g, or the exact order of g, depending on if optional parameters are specified.

- [<code>sample_r_given_N(N, factors)</code>](sample_r_given_N.md)

  Returns the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N, without explicitly computing and returning the element g.

