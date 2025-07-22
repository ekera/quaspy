## Function: <code>sample\_r\_given\_N(N, factors)</code>
Returns the order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N, without explicitly computing and returning the element g.

Suppose that N = p1^e1 * ... * pn^en, for p1, ..., pn pairwise distinct odd prime factors, and e1, ..., en positive integer exponents.

For i in [1, n], suppose that gi is selected uniformly at random from the multiplicative group of the ring of integers modulo pi^ei, and that the order ri of gi is computed. Then r = lcm(r1, ..., rn) is the order of g, where g may be computed via the Chinese remainder theorem, by using that it must hold that gi = g mod pi^ei for i in [1, n].

(Note that the above map is an isomorphism: Hence g will be selected uniformly at random from the multiplicative group of the ring of integers modulo N, as the gi are selected uniformly at random from the multiplicative group of the ring of integers modulo pi^ei.)

A problem with the above approach is that it is hard to compute the order ri, and hence to compute r, unless the factorization of pi - 1 for i in [1, n] is known. To circumvent this problem, instead of directly selecting the gi uniformly at random from the multiplicative group of the ring of integers modulo pi^ei, suppose exponents di are instead selected uniformly at random on the interval [0, lambda(pi^ei)), where

lambda(pi^ei) = (pi - 1) pi^(ei - 1)

as all pi are odd, and gi = Gi^di computed, for Gi some fixed generator of the multiplicative group the ring of integers modulo pi^ei. Then, gi is of order ri = lambda(pi^ei) / gcd(lambda(pi^ei), di) for i in [1, n], so the ri are easy to compute, as is r = lcm(r1, ..., rn).

Again, unless the factorization of pi - 1 for i in [1, n] is known, it is hard to prove that an element Gi is a generator, and hence to explicitly compute g. This explains why this function only returns r.

The above procedure is not described in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1), but in [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see in particular Sect. 5.2.3) and in the [Factoritall repository](https://www.github.com/ekera/factoritall) (available at https://www.github.com/ekera/factoritall).

## Import directive
```python
from quaspy.factoring.sampling import sample_r_given_N
```

## Parent module
- [<code>sampling</code>](README.md)

## Prototype
```python
def sample_r_given_N(N : int | gmpy2.mpz,
                     factors : list)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| N | The integer N. |
| factors | The factors of N = p1^e1 * ... * pn^en, represented on the form [[p1, e1], ..., [pn, en]], for p1, ..., pn pairwise distinct prime factors, and for e1, ..., en positive integer exponents. |

## Return value
The order r of an element g selected uniformly at random from the multiplicative group of the ring of integers modulo N.

