## Module: <code>postprocessing</code>
A module for factoring an integer N given the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N.

## Import directive
```python
import quaspy.factoring.general.postprocessing
```

## Parent module
- [<code>general</code>](../README.md)

## Submodules
- [<code>ekera</code>](ekera/README.md)

  A module for factoring N completely given the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N, where g need not be explicitly specified. This by using the algorithm in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).

- [<code>shor</code>](shor/README.md)

  A module for splitting N, given the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N, where g must be explicitly specified. This by using the original algorithm in [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700).

