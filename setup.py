import setuptools;

long_description  = "![Quaspy](https://raw.githubusercontent.com/ekera/quaspy/main/docs/quaspy.png)\n\n";

long_description += "# Quaspy\n";
long_description += "The [Quaspy](https://github.com/ekera/quaspy) (<i>Quantum algorithm simulations in Python</i>) library for [Python3](https://www.python.org) contains modules that implement:\n\n";

long_description += "- Simulators for the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791), and the classical post-processing in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791) that recovers the order in a single run with very high success probability.\n\n";
long_description += "- Simulators for factoring general integers via order-finding, and the post-processing in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) and [[E22p]](https://doi.org/10.48550/arXiv.2201.07791) that factors any integer completely in a single order-finding run with very high success probability.\n\n";
long_description += "- Simulators for the quantum part of Shor's algorithm for computing general discrete logarithms [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084), and the post-processing in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084) that recovers the logarithm given the order in a single run with very high probability of success.\n\n";
long_description += "- Simulators for the quantum part of Ekerå–Håstad's algorithm for computing short discrete logarithms [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20), modified as in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754), and the post-processing in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that recovers the logarithm in a single run with very high probability of success. This algorithm does not require the order to be known.\n\n";
long_description += "- Simulators for factoring RSA integers via short discrete logarithms, by using the reduction in [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20), modified as in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754), and the post-processing in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that factors random RSA integers in a single run with very high probability of success.\n\n";

long_description += "All modules, classes, methods and functions in [Quaspy](https://github.com/ekera/quaspy) are [documented](https://github.com/ekera/quaspy/blob/main/docs/README.md) using [Python3](https://www.python.org) docstrings.\n\n";
long_description += "Note that [Quaspy](https://github.com/ekera/quaspy) does not implement support for tradeoffs in the aforementioned algorithms. Support for tradeoffs may potentially be added in the future. For the time being, see instead the [Qunundrum](https://www.github.com/ekera/qunundrum) repository with its suite of MPI programs. Note furthermore that portions of [Quaspy](https://github.com/ekera/quaspy) are inherited from the [Factoritall](https://www.github.com/ekera/factoritall) repository.\n\n";
long_description += "[Quaspy](https://github.com/ekera/quaspy) is a work in progress, and may be subject to major changes without prior notice. [Quaspy](https://github.com/ekera/quaspy) was developed for academic research purposes. It grew out of our research project in an organic manner as research questions were posed and answered. It is distributed \"as is\" without warranty of any kind, either expressed or implied. For further details, see the [license](https://github.com/ekera/quaspy/blob/main/LICENSE.md).\n\n";

long_description += "## Examples\n";
long_description += "For examples that illustrate how to use [Quaspy](https://github.com/ekera/quaspy), please see the [<code>examples</code>](https://github.com/ekera/quaspy/blob/main/examples) directory in the [Quaspy](https://github.com/ekera/quaspy) repository.\n\n";
long_description += "See also the [documentation](https://github.com/ekera/quaspy/blob/main/docs/README.md) for [Quaspy](https://github.com/ekera/quaspy) for help on how to use the library.\n\n";

long_description += "## About and acknowledgments\n";
long_description += "The [Quaspy](https://github.com/ekera/quaspy) library was developed by [Martin Ekerå](mailto:ekera@kth.se), in part at [KTH, the Royal Institute of Technology](https://www.kth.se/en), in Stockholm, [Sweden](https://www.sweden.se). Valuable comments and advice were provided by Johan Håstad throughout the development process.\n\n";
long_description += "Funding and support was provided by the Swedish NCSA that is a part of the [Swedish Armed Forces](https://www.mil.se).";

setuptools.setup(
  name="quaspy",
  version="0.9.2",
  author="Martin Ekerå",
  author_email="ekera@kth.se",
  description="A package for post-processing Shor's order-finding and factoring algorithms.",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/ekera/quaspy",
  project_urls={
    "Bug Tracker": "https://github.com/ekera/quaspy/issues",
  },
  classifiers=[
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
  ],
  package_dir={"": "src"},
  packages=setuptools.find_packages(where="src"),
  python_requires=">=3.6",
  install_requires=['gmpy2', 'secret']
);
