import setuptools;

long_description  = "![Quaspy](https://github.com/ekera/quaspy/blob/main/docs/quaspy.png)\n\n";
long_description += "# Quaspy\n";
long_description += "[Quaspy](https://github.com/ekera/quaspy) is a [Python3](https://www.python.org) library. [Quaspy](https://github.com/ekera/quaspy) is an abbreviation of <i>Quantum algorithm simulations in Python</i>.\n\n";
long_description += "[Quaspy](https://github.com/ekera/quaspy) contains modules that implement:\n\n";
long_description += "- Simulators for the quantum part of Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791), and the classical post-processing in [[E22p]](https://doi.org/10.48550/arXiv.2201.07791) that recovers the order in a single run with very high success probability.\n\n";
long_description += "- Simulators for factoring general integers via order-finding, and the post-processing in [[E21b]](https://doi.org/10.1007/s11128-021-03069-1) and [[E22p]](https://doi.org/10.48550/arXiv.2201.07791) that factors any integer completely in a single order-finding run with very high success probability.\n\n";
long_description += "- Simulators for the quantum part of Shor's algorithm for computing general discrete logarithms [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700), modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084), and the post-processing in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084) that recovers the logarithm given the order in a single run with very high probability of success.\n\n";
long_description += "- Simulators for the quantum part of Ekerå–Håstad's algorithm for computing short discrete logarithms [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20), modified as in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754), and the post-processing in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that recovers the logarithm in a single run with very high probability of success. This algorithm does not require the order to be known.\n\n";
long_description += "- Simulators for factoring RSA integers via short discrete logarithms, by using the reduction in [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20), modified as in [[E20]](https://doi.org/10.1007/s10623-020-00783-2) and [[E23p]](https://doi.org/10.48550/arXiv.2309.01754), and the post-processing in [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) that factors random RSA integers in a single run with very high probability of success.\n\n";
long_description += "For further details, please see the [Quaspy repository](https://github.com/ekera/quaspy) on [GitHub](https://github.com).";

setuptools.setup(
  name="quaspy",
  version="0.9.0",
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
