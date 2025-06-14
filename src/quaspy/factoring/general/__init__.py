""" @brief  A module for factoring general integers.

    This module uses Shor's order-finding algorithm to factor general integers,
    as described in [Shor94], with improvements from [E21b] and [E24].

    [Shor94] Shor, P.W.: "Algorithms for Quantum Computation: Discrete
                          Logarithms and Factoring".
                         In: Proceedings from FOCS '94, pp. 124–134 (1994).

    [E21b] Ekerå, M.: "On completely factoring any integer efficiently in a
                       single run of an order-finding algorithm".
                      Quantum Inf. Process. 20(6):205 (2021).

    [E24] Ekerå, M.: "On the success probability of quantum order finding".
                     ACM Trans. Quantum Comput. 5(2):11 (2024). """
