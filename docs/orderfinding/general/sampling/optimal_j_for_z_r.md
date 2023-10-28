## Function: <code>optimal\_j\_for\_z\_r(z, r, m, l)</code>
Computes and returns the optimal frequency j0(z) for z in [0, r).

## Import directive
```python
from quaspy.orderfinding.general.sampling import optimal_j_for_z_r
```

## Parent module
- [<code>sampling</code>](README.md)

## Prototype
```python
def optimal_j_for_z_r(z,
                      r,
                      m,
                      l)
```

## Parameters
| <b>Name</b> | <b>Description</b> |
| ----------- | ------------------ |
| z | The peak index z in [0, r). |
| r | The order r. |
| m | A positive integer m such that r < 2^m. |
| l | An integer l on [0, m), such that m + l is the length of the control register in the quantum order-finding algorithm. |

## Return value
The optimal frequency j0(z) = round(2^(m + l) / r * z).

