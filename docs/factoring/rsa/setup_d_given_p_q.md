## Function: <code>setup\_d\_given\_p\_q(p, q)</code>
Sets up d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1) given p and q.

This is a convenience function for Ekerå–Håstad's algorithm [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) that factors a large random RSA integer N = pq into p and q.

More specifically, as is explained in App. A.2 of [[E20]](https://doi.org/10.1007/s10623-020-00783-2), it holds that

x = g^d' = g^((N - 1) / 2 - 2^(l - 1)) = g^((p - 1) / 2 + (q - 1) / 2 - 2^(l - 1)) = g^d

Provided that the order r of g is sufficiently large, it holds that

d = d' mod r = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1),

allowing d to be computed as the discrete logarithm of x to the base g.

This may be done efficiently by using Ekerå–Håstad's quantum algorithm that computes short discrete logarithms in groups of unknown order.

Given d = (p - 1) / 2 + (q - 1) / 2 - 2^(l - 1) and N = pq, it is then trivial to compute p and q by solving a quadratic equation.

This convenience function computes d in the above procedure given p and q.

## Import directive
```python
from quaspy.factoring.rsa import setup_d_given_p_q
```

## Parent module
- [<code>rsa</code>](README.md)

## Prototype
```python
def setup_d_given_p_q(p,
                      q)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| p | A large random l-bit prime. |
| q | A large random l-bit prime not equal to p. |

## Return value
The logarithm d = ((p - 1) / 2) + ((q - 1) / 2) - 2^(l - 1).

