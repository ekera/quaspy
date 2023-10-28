## Function: <code>setup\_x\_given\_g\_N(g, N)</code>
Sets up x = g^d' for d' = (N - 1) / 2 - 2^(l - 1) given g and N, for N the product of two large random distinct l-bit primes.

This is a convenience function for Ekerå–Håstad's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) that factors a large random RSA integer N = pq into p and q.

More specifically, as is explained in App. A.2 of [[E20]](https://doi.org/10.1007/s10623-020-00783-2), it holds that

x = g^d' = g^((N - 1) / 2 - 2^(l - 1)) = g^((p - 1) / 2 + (q - 1) / 2 - 2^(l - 1)) = g^d.

Provided that the order r of g is sufficiently large, it holds that

d = d' mod r = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1),

allowing d to be computed as the discrete logarithm of x to the base g.

This may be done efficiently by using Ekerå–Håstad's quantum algorithm that computes short discrete logarithms in groups of unknown order.

Given d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1) and N = pq, it is then trivial to compute p and q by solving a quadratic equation.

This convenience function computes x in the above procedure given g and N.

## Import directive
```python
from quaspy.factoring.rsa import setup_x_given_g_N
```

## Parent module
- [<code>rsa</code>](README.md)

## Prototype
```python
def setup_x_given_g_N(g,
                      N)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| g | The generator g selected uniformly at random from the multiplicative group of the ring of integers modulo N. |
| N | The integer N = pq, for p and q two large distinct random l-bit prime numbers. |

## Return value
The element x = g^((N - 1) / 2 - 2^(l - 1)).

