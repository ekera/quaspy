## Class: <code>CandidateCollection</code>
A class representing a collection of candidates for r_tilde, that is kept reduced in the sense that if r_tilde is in the collection, then so is c * r_tilde, for c any positive integer.

The idea is to only add candidates for r_tilde that meet the requirement that g^(e * r_tilde) = 1 to the collection, for e a product of cm-smooth prime powers. If r_tilde meets this requirement, then of course so does c * r_tilde, so it makes sense to keep the collection reduced.

In the internal set of candidates that represent the collection, only the minimum set of candidates for r_tilde required to represent the collection are stored. The collection furthermore has features for testing if a candidate is already in the collection. This is useful when searching for r_tilde, as it reduces the number of exponentiations required, and the amount of storage needed to represent candidates, and so forth.

## Import directive
```python
from quaspy.orderfinding.general.postprocessing.ekera.internal.collection import CandidateCollection
```

## Parent module
- [<code>collection</code>](README.md)

## Methods
- [<code>\_\_contains\_\_(self, candidate)</code>](CandidateCollection/__contains__.md)

  Checks if a candidate for r_tilde is in this collection.

- [<code>\_\_init\_\_(self)</code>](CandidateCollection/__init__.md)

  Initializes the collection of candidates for r_tilde.

- [<code>\_\_iter\_\_(self)</code>](CandidateCollection/__iter__.md)

  Returns an iterator for the candidates for r_tilde that represent this collection.

- [<code>\_\_len\_\_(self)</code>](CandidateCollection/__len__.md)

  Returns the number of candidates for r_tilde that represent this collection.

- [<code>\_\_repr\_\_(self)</code>](CandidateCollection/__repr__.md)

  Returns a string representation of this collection.

- [<code>\_\_str\_\_(self)</code>](CandidateCollection/__str__.md)

  Returns a string representation of this collection.

- [<code>add(self, candidate)</code>](CandidateCollection/add.md)

  Adds a candidate for r_tilde to this collection.

