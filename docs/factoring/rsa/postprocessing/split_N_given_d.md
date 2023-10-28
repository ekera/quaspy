## Function: <code>split\_N\_given\_d(d, N)</code>
Splits N = pq into the large l-bit prime factors p and q given d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1).

This is a convenience function for Ekerå–Håstad's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) that factors a large random RSA integer N = pq into p and q.

More specifically, as is explained in App. A.2 of [[E20]](https://doi.org/10.1007/s10623-020-00783-2), it holds that

x = g^d' = g^((N - 1) / 2 - 2^(l - 1)) = g^((p - 1) / 2 + (q - 1) / 2 - 2^(l - 1)) = g^d

Provided that the order r of g is sufficiently large, it holds that

d = d' mod r = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1),

allowing d to be computed as the discrete logarithm of x to the base g.

This may be done efficiently by using Ekerå–Håstad's quantum algorithm that computes short discrete logarithms in groups of unknown order.

Given d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1) and N = pq, it is then trivial to compute p and q by solving a quadratic equation.

This convenience function performs the last step in the above procedure.

## Import directive
```python
from quaspy.factoring.rsa.postprocessing import split_N_given_d
```

## Parent module
- [<code>postprocessing</code>](README.md)

## Prototype
```python
def split_N_given_d(d,
                    N)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| d | The logarithm d. |
| N | The integer N. It is required that N is odd, and the product of two large random distinct l-bit prime numbers. |

## Return value
A set {p, q} of two non-trivial factors of N such that N = pq, or None if the algorithm failed to split N.

