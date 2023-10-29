![Quaspy](docs/quaspy.png)

# Quaspy
The [Quaspy](https://github.com/ekera/quaspy) (<i>Quantum algorithm simulations in Python</i>) library for [Python3](https://www.python.org) contains modules that implement:

- Simulators for the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791), and the classical post-processing in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791) that recovers the order in a single run with very high success probability.

   [Examples](examples/orderfinding.ipynb) | [Documentation](docs/orderfinding/general/README.md)

- Simulators for factoring general integers via order-finding, and the post-processing in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) and [[E22p]](https://doi.org/10.48550/arXiv.2201.07791) that factors any integer completely in a single order-finding run with very high success probability.

   [Examples](examples/factoring.ipynb) | [Documentation](docs/factoring/general/README.md)

- Simulators for the quantum part of Shor's algorithm for computing general discrete logarithms [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084), and the post-processing in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084) that recovers the logarithm given the order in a single run with very high probability of success.

   [Examples](examples/logarithmfinding.ipynb) | [Documentation](docs/logarithmfinding/general/README.md)

- Simulators for the quantum part of Ekerå–Håstad's algorithm for computing short discrete logarithms [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20), modified as in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754), and the post-processing in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that recovers the logarithm in a single run with very high probability of success. This algorithm does not require the order to be known.

   [Examples](examples/logarithmfinding-short.ipynb) | [Documentation](docs/logarithmfinding/short/README.md)

- Simulators for factoring RSA integers via short discrete logarithms, by using the reduction in [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20), modified as in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754), and the post-processing in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that factors random RSA integers in a single run with very high probability of success.

   [Examples](examples/factoring-rsa.ipynb) | [Documentation](docs/factoring/rsa/README.md)

All modules, classes, methods and functions in [Quaspy](https://github.com/ekera/quaspy) are [documented](docs/README.md) using [Python3](https://www.python.org) docstrings.

Note that [Quaspy](https://github.com/ekera/quaspy) does not implement support for tradeoffs in the aforementioned algorithms. Support for tradeoffs may potentially be added in the future. For the time being, see instead the [Qunundrum](https://www.github.com/ekera/qunundrum) repository with its suite of MPI programs. Note furthermore that portions of [Quaspy](https://github.com/ekera/quaspy) are inherited from the [Factoritall](https://www.github.com/ekera/factoritall) repository.

[Quaspy](https://github.com/ekera/quaspy) is a work in progress, and may be subject to major changes without prior notice.
[Quaspy](https://github.com/ekera/quaspy) was developed for academic research purposes. It grew out of our research project in an organic manner as research questions were posed and answered. It is distributed "as is" without warranty of any kind, either expressed or implied. For further details, see the [license](LICENSE.md).

## Prerequisites
To install [Python3](https://www.python.org) under [Ubuntu 22.04 LTS](https://releases.ubuntu.com/22.04), along with required dependencies, execute:

```console
$ sudo apt install python3 python3-pip python3-venv
$ sudo apt install libgmp-dev libmpfr-dev libmpc-dev
```

For other Linux and Unix distributions, or operating systems, you may need to [download Python3](https://www.python.org/download) and install it manually along with the required dependencies.

### Installing the library
To install the [Quaspy](https://github.com/ekera/quaspy) library, execute:

```console
$ pip3 install dist/quaspy-0.9.2-py3-none-any.whl
```

You may also install [Quaspy](https://github.com/ekera/quaspy) via [Pip3](https://pypi.org) by executing:

```console
$ pip3 install quaspy
```

## Examples
For examples that illustrate how to use [Quaspy](https://github.com/ekera/quaspy), please see the [<code>examples</code>](examples) directory.

See also the [documentation](docs/README.md) for [Quaspy](https://github.com/ekera/quaspy) for help on how to use the library.

## About and acknowledgments
The [Quaspy](https://github.com/ekera/quaspy) library was developed by [Martin Ekerå](mailto:ekera@kth.se), in part at [KTH, the Royal Institute of Technology](https://www.kth.se/en), in Stockholm, [Sweden](https://www.sweden.se). Valuable comments and advice were provided by Johan Håstad throughout the development process.

Funding and support was provided by the Swedish NCSA that is a part of the [Swedish Armed Forces](https://www.mil.se).