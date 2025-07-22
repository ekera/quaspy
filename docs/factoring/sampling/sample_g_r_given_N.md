## Function: <code>sample\_g\_r\_given\_N(N, N_factors, ..)</code>
Returns [g, r], for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N, and r either a heuristic estimate of the order of g, or the exact order of g, depending on if optional parameters are specified.

Suppose that N = p1^e1 * ... * pn^en, for p1, ..., pn pairwise distinct odd prime factors, and e1, ..., en positive integer exponents.

For i in [1, n], this function then selects gi from the multiplicative group of the ring of integers modulo pi^ei.

1. If the factorization of pi - 1 for i in [1, n] is *not* specified, this function then heuristically estimates the order ri of gi by using the method in App. A of [[E21b]](https://doi.org/10.1007/s11128-021-03069-1): Specificially, by using that

lambda(pi^ei) = (pi - 1) pi^(ei - 1),

as pi is odd, and by using a factor base of primes <= B to find all small factors of pi - 1 via trial division. It then computes g via the Chinese remainder theorem, by requiring that gi = g mod pi^ei, along with a heuristic estimate r = lcm(r1, ..., rn) of the order of g, that is correct with high probability, as is explained in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1).

2. If the factorization of pi - 1 for i in [1, n] is specified, this function then exactly computes the order ri of gi. Specifically, by using

lambda(pi^ei) = (pi - 1) pi^(ei - 1)

as an initial guess ri' for the order ri of gi. Then, for each prime factor f that divide ri', for as long as f divides ri' and gi^(ri' / f) = 1 (mod N), let ri' <- ri' / f. It follows that ri = ri' at the end of the procedure. The order of g is then r = lcm(r1, ..., rn).

The above procedure is described in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1), and in the [factoritall repository](https://www.github.com/ekera/factoritall) (available at https://www.github.com/ekera/factoritall).

## Import directive
```python
from quaspy.factoring.sampling import sample_g_r_given_N
```

## Parent module
- [<code>sampling</code>](README.md)

## Prototype
```python
def sample_g_r_given_N(N : int | gmpy2.mpz,
                       N_factors : list,
                       pi_minus_one_factors : list[list[list[int | gmpy2.mpz]]] | None = None,
                       B : int = 1000000)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| N | The integer N. |
| N_factors | The factors of N = p1^e1 * ... * pn^en, represented on the form [[p1, e1], ..., [pn, en]], for p1, ..., pn pairwise distinct prime factors, and for e1, ..., en positive integer exponents. |
| pi_minus_one_factors | The factors of pi-1 = qi1^di1 * ... * qim^dim, for i in [1, n], represented on the form [F1, ..., Fn], where each Fi is on the form [[qi1, qi1], ..., [qim, qim]], for qi1, ..., qim pairwise distinct prime factors, and for di1, ..., dim positive integer exponents. May be set to None, in which case r will be computed deterministically as described above. If explicitly specified, the order r will be computed exactly. |
| B | The upper bound on the prime factors to consider when performing trial division. Has no effect if pi_minus_one_factors is explicitly specified, as trial division is then not performed. |

## Return value
The pair [g, r], for g an element selected uniformly at random from the multiplicative group of the ring of integers modulo N, and r a heuristic estimate of the order of g, or the exact order of g, depending on if optional parameters are specified.

