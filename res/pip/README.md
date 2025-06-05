![Quaspy](https://raw.githubusercontent.com/ekera/quaspy/main/docs/quaspy.png)

# The Quaspy library for Python
The [Quaspy](https://github.com/ekera/quaspy) (<i>Quantum algorithm simulations in Python</i>) library for [Python](https://www.python.org) contains modules that implement:

- Simulators for the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E24]](https://doi.org/10.1145/3655026), and the classical post-processing in [[E24]](https://doi.org/10.1145/3655026) that recovers the order in a single run with very high success probability.

- Simulators for factoring general integers via order-finding, and the post-processing in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) and [[E24]](https://doi.org/10.1145/3655026) that factors any integer completely in a single order-finding run with very high success probability.

- Simulators for the quantum part of Shor's algorithm for computing general discrete logarithms [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084), and the post-processing in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084) that recovers the logarithm given the order in a single run with very high probability of success.

- Simulators for the quantum part of Ekerå–Håstad's algorithm for computing short discrete logarithms [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20), modified as in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754), and the post-processing in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that recovers the logarithm in a single run with very high probability of success. This algorithm does not require the order to be known.

- Simulators for factoring RSA integers via short discrete logarithms, by using the reduction in [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20), modified as in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754), and the post-processing in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that factors random RSA integers in a single run with very high probability of success.

All modules, classes, methods and functions in [Quaspy](https://github.com/ekera/quaspy) are [documented](https://github.com/ekera/quaspy/blob/main/docs/README.md) using [Python](https://www.python.org) docstrings.

Note that [Quaspy](https://github.com/ekera/quaspy) does not implement support for tradeoffs in the aforementioned algorithms. Support for tradeoffs may potentially be added in the future. For the time being, see instead the [Qunundrum](https://www.github.com/ekera/qunundrum) repository with its suite of MPI programs. Note furthermore that portions of [Quaspy](https://github.com/ekera/quaspy) are inherited from the [Factoritall](https://www.github.com/ekera/factoritall) repository.

[Quaspy](https://github.com/ekera/quaspy) is a work in progress, and may be subject to major changes without prior notice. [Quaspy](https://github.com/ekera/quaspy) was developed for academic research purposes. It grew out of our research project in an organic manner as research questions were posed and answered. It is distributed "as is" without warranty of any kind, either expressed or implied. For further details, see the [license](https://github.com/ekera/quaspy/blob/main/LICENSE.md).

## Examples
For examples that illustrate how to use [Quaspy](https://github.com/ekera/quaspy), please see the [<code>examples</code>](https://github.com/ekera/quaspy/blob/main/examples) directory in the [Quaspy](https://github.com/ekera/quaspy) repository.

See also the [documentation](https://github.com/ekera/quaspy/blob/main/docs/README.md) for [Quaspy](https://github.com/ekera/quaspy) for help on how to use the library.

## About and acknowledgments
The [Quaspy](https://github.com/ekera/quaspy) library was developed by [Martin Ekerå](mailto:ekera@kth.se), in part at [KTH, the Royal Institute of Technology](https://www.kth.se/en), in Stockholm, [Sweden](https://www.sweden.se). Valuable comments and advice were provided by Johan Håstad throughout the development process.

Funding and support was provided by the Swedish NCSA that is a part of the [Swedish Armed Forces](https://www.mil.se).