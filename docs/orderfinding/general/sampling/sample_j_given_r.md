## Function: <code>sample\_j\_given\_r(r, m, l, ..)</code>
Samples a frequency j from the distribution induced by the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) for a given order r. This by using the sampling procedure described in [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf) (see in particular Sects. 5.3.5 and 5.4.3).

This fuction can also be used to sample from the distribution induced by Seifert's variation [[Seifert01]](https://doi.org/10.1007/3-540-45353-9_24) of Shor's order-finding algorithm for a given order r, depending on how parameters are selected.

Throughout this function, the algorithms are as described in [[E21]](https://doi.org/10.1515/jmc-2020-0006), [[E24]](https://doi.org/10.1145/3655026) and [[E24t]](https://diva-portal.org/smash/get/diva2:1902626/FULLTEXT01.pdf). The notation is also inherited from said works.

To sample the distribution, one of the r peaks is first picked uniformly at random by picking an index z uniformly at random from [0, r). This function furthermore selects a pivot uniformly at random from [0, 1).

For j0(z) the optimal frequency for this peak, it then subtracts the probabilities r * P((j0(z) + offset) mod 2^(m + l)) of frequencies that are offset by 0, ±1, ±2, ..., ±(B - 1) from the optimal offset from the pivot, returning j = (j0(z) + offset) mod 2^(m + l) as soon as pivot <= 0.

## Import directive
```python
from quaspy.orderfinding.general.sampling import sample_j_given_r
```

## Parent module
- [<code>sampling</code>](README.md)

## Prototype
```python
def sample_j_given_r(r : int | gmpy2.mpz,
                     m : int,
                     l : int,
                     B : int = 100000,
                     timeout : int | None | quaspy.utils.timeout.Timeout = None,
                     verbose : bool = False,
                     extended_result : bool = False)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| r | The order r. |
| m | A positive integer m such that r < 2^m. |
| l | An integer l on [0, m), such that m + l is the length of the control register in the quantum part of the algorithm. |
| B | A parameter B that upper-bounds the offset from the optimal frequency j0(z) as explained above. |
| timeout | A timeout after which a TimeoutError will be raised and the sampling procedure aborted.<br><br>The timeout may be represented as an integer specifying the timeout in seconds, or as an instance of the Timeout class. May be set to None, as is the default, in which case no timeout is enforced. |
| verbose | A flag that may be set to True to print intermediary results when sampling. |
| extended_result | A flag that may be set to True to not only return the frequency j, but [j, [z, j0(z), offset]]. |

## Return value
The frequency j sampled if the extended_result flag is set to False, or [j, [z, j0(z), offset]] if the extended_result flag is set to True, or None if sampling failed because the upper bound B on the offset from the optimal frequency j0(z) was reached.

