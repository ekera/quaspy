## Module: <code>primes</code>
A module for computing prime ranges and prime power products, for sampling random primes and for testing if integers are smooth.

## Import directive
```python
import quaspy.math.primes
```

## Parent module
- [<code>math</code>](../README.md)

## Functions
- [<code>is_B_smooth(d, B)</code>](is_B_smooth.md)

  Tests if the integer d is B-smooth.

- [<code>prime_power_product(B)</code>](prime_power_product.md)

  Returns the product of q^e, as q runs over all primes <= B, for e the largest non-negative integer exponent such that q^e <= B.

- [<code>prime_range(B)</code>](prime_range.md)

  Returns an ordered list of all primes less than B.

- [<code>sample_l_bit_prime(l)</code>](sample_l_bit_prime.md)

  Returns an l-bit prime selected uniformly at random from the set of all such primes.

