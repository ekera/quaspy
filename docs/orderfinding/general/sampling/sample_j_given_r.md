## Function: <code>sample\_j\_given\_r(r, m, l, ..)</code>
Samples a frequency j from the distribution induced by Shor's order-finding algorithm for a given order r.

To sample the distribution, one of the r peaks is first picked uniformly at random by picking an index z uniformly at random from [0, r). This function furthermore selects a pivot uniformly at random from [0, 1).

For j0(z) the optimal frequency for this peak, it then subtracts the probabilities r * P((j0(z) + offset) mod 2^(m + l)) of frequencies that are offset by 0, ±1, ±2, .., ±(B - 1) from the optimal offset from the pivot, returning j = (j0(z) + offset) mod 2^(m + l) as soon as pivot <= 0.

## Import directive
```python
from quaspy.orderfinding.general.sampling import sample_j_given_r
```

## Parent module
- [<code>sampling</code>](README.md)

## Prototype
```python
def sample_j_given_r(r,
                     m,
                     l,
                     B = 1000,
                     verbose = False,
                     extended_result = False)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| r | The order r. |
| m | A positive integer m such that r < 2^m. |
| l | An integer l on [0, m), such that m + l is the length of the control register in the quantum order-finding algorithm. |
| B | A parameter B that upper-bounds the offset from the optimal frequency j0(z) as explained above. |
| verbose | A flag that may be set to True to print intermediary results when sampling. |
| extended_result | A flag that may be set to True to not only return the frequency j, but [j, [z, j0(z), offset]]. |

## Return value
The frequency j sampled if the extended_result flag is set to False, or [j, [z, j0(z), offset]] if the extended_result flag is set to True, or None if sampling failed because the upper bound B on the offset from the optimal frequency j0(z) was reached.

