""" @brief  The root module for the Quaspy (Quantum algorithm simulations in
            Python) library for Python.

    Quaspy contains modules that implement:

    - Simulators for the quantum part of Shor's order-finding algorithm
      [Shor94], modified as in [E24], and the classical post-processing from
      [E24] that recovers the order in a single run with very high success
      probability. Furthermore, these modules implement simulators for the
      quantum part of Seifert's variation [Seifert01] of Shor's order-finding
      algorithm [Shor94], as described as in [E24t] and [E21], and with the
      classical post-processing from [E24t] and [E21] with supporting functions
      from [E24] that efficiently recovers the order, or a positive integer
      multiple of the order, in multiple runs when making tradeoffs.

    - Simulators for factoring general integers via order-finding, and the
      classical post-processing from [E21b] and [E24] that factors any integer
      completely in a single order-finding run with very high success
      probability.

    - Simulators for the quantum part of Shor's algorithm for computing general
      discrete logarithms [Shor94], modified as in [E19p], and the classical
      post-processing from [E19p] that recovers the logarithm given the order in
      a single run with very high probability of success, and that also
      efficiently recovers the logarithm given the order when making tradeoffs.

    - Simulators for the quantum part of Ekerå–Håstad's algorithm for computing
      short discrete logarithms [EH17], modified as in [E20] and [E23p], and the
      classical post-processing from [E23p] that recovers the logarithm in a
      single run with very high probability of success. This algorithm does not
      require the order to be known. Furthermore, these modules implement the
      classical post-processing from [E20] that efficiently recovers the
      logarithm in multiple runs of the quantum part of Ekerå–Håstad's algorithm
      when making tradeoffs.

    - Simulators for factoring RSA integers via short discrete logarithms, by
      using the reduction in [EH17], modified as in [E20] and [E23p], and the
      classical post-processing from [E23p] that factors random RSA integers in
      a single run of the quantum part of Ekerå–Håstad's algorithm with very
      high probability of success. Furthermore, these modules implement the
      classical post-processing from [E20] that efficiently factors random RSA
      integers in multiple runs of the quantum part of Ekerå–Håstad's algorithm
      when making tradeoffs.

    All modules, classes, methods and functions in Quaspy are documented using
    Python docstrings (see https://github.com/ekera/quaspy/blob/main/docs).

    Note that Quaspy implements basic support for tradeoffs via a native Python
    implementation of LLL that is stable and resasonable performant. See also
    the Qunundrum repository (see https://github.com/ekera/qunundrum) with its
    suite of MPI programs that implements support for tradeoffs via LLL and BKZ
    as implemented by fpLLL (see https://github.com/fplll/fplll).

    Note furthermore that portions of Quaspy are inherited from the Factoritall
    repository (see https://github.com/ekera/factoritall).

    Quaspy is a work in progress, and may be subject to major changes without
    prior notice. Quaspy was developed for academic research purposes. It grew
    out of our research project in an organic manner as research questions were
    posed and answered. It is distributed "as is" without warranty of any kind,
    either expressed or implied. For further details, see the license (see
    https://github.com/ekera/quaspy/blob/main/LICENSE.md).

    Quaspy was developed by Martin Ekerå, in part at KTH, the Royal Institute of
    Technology, in Stockholm, Sweden. Valuable comments and advice were provided
    by Johan Håstad throughout the development process. Funding and support was
    provided by the Swedish NCSA that is a part of the Swedish Armed Forces.

    For further details on Quaspy, see the Quaspy repository on GitHub
    (available at https://github.com/ekera/quaspy).

    [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                          Logarithms and Factoring".
                         In: Proceedings from FOCS '94, pp. 124–134 (1994).

    [Seifert01] Seifert, J.-P.: "Using fewer qubits in Shor's factorization
                                 algorithm via simultaneous Diophantine
                                 approximation". In: CT-RSA 2001.
                                Springer LNCS 2020, pp. 319–227 (2001).

    [EH17] Ekerå, M. and Håstad, J.: "Quantum Algorithms for Computing Short
                                      Discrete Logarithms and Factoring RSA
                                      Integers.". In: PQCrypto 2017.
                                     Springer LNCS 10346, pp. 347–363 (2017).

    [E19p] Ekerå, M.: "Revisiting Shor's quantum algorithm for computing
                       general discrete logarithms".
                      ArXiv 1905.09084v4 (2024).

    [E20] Ekerå, M.: "On post-processing in the quantum algorithm for
                      computing short discrete logarithms".
                     Des. Codes Cryptogr. 88, pp. 2313–2335 (2020).

    [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                       single run of an order-finding algorithm".
                      Quantum Inf. Process. 20(6):205 (2021).

    [E23p] Ekerå, M.: "On the success probability of the quantum algorithm for
                       the short DLP". ArXiv 2309.01754v2 (2025).

    [E24] Ekerå, M.: "On the success probability of quantum order finding".
                     ACM Trans. Quantum Comput. 5(2):11 (2024).

    [E24t] Ekerå, M.: "On factoring integers, and computing discrete logarithms
                       and orders, quantumly".
                      PhD thesis, KTH Royal Institute of Technology (2024). """

# Shor's algorithm for splitting N if r is even and g^(r/2) != -1 (mod N).
from .factoring.general.postprocessing.shor import split_N_given_g_r;

# The algorithm from [E21b] for completely factoring N given r.
from .factoring.general.postprocessing.ekera import solve_r_for_factors;

# Convenience functions implementing algorithms from [E21b] and [E24] for first
# solving j for r, and then solving r and N for the complete factorization of N.
from .factoring.general.postprocessing.ekera import solve_j_for_factors;
from .factoring.general.postprocessing.ekera import solve_j_for_factors_mod_N;

# Convenience function implementing the algorithm from [EH17], with improvements
# from [E20], for splitting N given d by solving a quadratic equation.
from .factoring.rsa.postprocessing import split_N_given_d;

# Convenience function implementing the algorithm from [EH17], with improvements
# from [E20], for setting up x = g^d given g and N = pq, and for computing d
# given the primes p and q.
from .factoring.rsa import setup_d_given_p_q;
from .factoring.rsa import setup_x_given_g_N;

# The algorithms in [E21b], [E24t] (see in particular Sect. 5.2.3) and the
# Factoritall repository (available at https://github.com/ekera/factoritall) for
# sampling g uniformly at random from the multiplicative group of the ring of
# integers modulo N, and for returning [g, r] or r, for r the (in some cases
# heuristically computed) order of g.
from .factoring.sampling import sample_g_r_given_N;
from .factoring.sampling import sample_r_given_N;

# The post-processing algorithms in [E24] for solving j for r.
from .orderfinding.general.postprocessing.ekera import solve_j_for_r;
from .orderfinding.general.postprocessing.ekera import solve_j_for_r_mod_N;

# Algorithms for sampling j given r to test the above order-finding algorithms.
from .orderfinding.general.sampling import sample_j_given_r;

# The post-processing algorithm from [E19p] for solving (j, k) for d given r.
from .logarithmfinding.general.postprocessing import solve_j_k_for_d_given_r;

# Algorithms for heuristically sampling (j, k) given d and r to test the above
# post-processing algorithm.
from .logarithmfinding.general.sampling import sample_j_k_given_d_r_heuristic;

# The post-processing algorithm from [E23p] for solving (j, k) for d.
from .logarithmfinding.short.postprocessing import solve_j_k_for_d;

# Algorithms for heuristically sampling (j, k) given d, r and tau to test the
# above post-processing algorithm.
from .logarithmfinding.short.sampling import sample_j_k_given_d_r_tau;
