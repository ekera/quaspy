""" @brief  A module for collecting the non-trivial factors of N. """

from gmpy2 import iroot;
from gmpy2 import is_prime;

from gmpy2 import mpz;
from gmpy2 import gcd;

from ......utils.timer import Timer;

class FactorCollection:

  """ @brief  A class for collecting the non-trivial factors of N. """

  def __init__(self, N):

    """ @brief  Initializes the collection of non-trivial factors of N.

        @param N  The integer N. """

    # The integer N to be factored.
    self.N = mpz(N);

    # The set of factors found thus far, reduced so that all factors in the set
    # are pairwise coprime to each other. This property is enforced by add().
    self.found_factors = set();

    # The set of prime factors found thus far; a subset of found_factors.
    self.found_primes = set();

    # A timer for measuring the time spent performing primality tests.
    self.timer_test_primality = Timer();

    # A timer for measuring the time spent detecting perfect powers.
    self.timer_test_perfect_power = Timer();

    # The residual; the product of the composite pairwise coprime factors in the
    # collection, or one if there are no composite factors in the collection.
    self.residual = 1;

    # Add N as a factor.
    self.add(N);

  def is_complete(self):

    """ @brief  Checks if all prime factors of N have been found.

        @return True if all prime factors have been found, False otherwise. """

    return self.residual == 1;

  def add(self, d):

    """ @brief  Adds a factor d to this collection.

        @param d  The factor to add to this collection. """

    d = mpz(d);

    # Check that the factor is non-trivial and has not already been found.
    if (d == 1) or (d in self.found_factors):
      return;

    # Test if d shares a factor with any of the factors found.
    D = mpz(1);

    for f in self.found_factors:
      D = gcd(f, d);

      if D != 1:
        break;

    if D != 1:
      # If so, remove f, split f and d, and add the resulting factors.
      self.found_factors.remove(f);
      if f not in self.found_primes:
        # Also remove f from the residual when removing f from the collection.
        self.residual //= f;

      f //= D;
      d //= D;

      self.add(D);

      if f != 1:
        self.add(f);

      if d != 1:
        self.add(d);
    else:
      # Check if d is a perfect power, and if so reduce d.
      self.timer_test_perfect_power.start();

      if d.is_power():
        for e in range(2, d.bit_length()):
          (d_reduced, result) = iroot(d, e);
          if result:
            d = d_reduced;
            break;

      self.timer_test_perfect_power.stop();

      # Add d to the factors found.
      self.found_factors.add(d);

      # Check if d is prime, and if so register it.
      self.timer_test_primality.start();
      result = is_prime(d);
      self.timer_test_primality.stop();

      if result:
        self.found_primes.add(d);
      else:
        # If d is not prime, multiply d onto the residual.
        self.residual *= d;

  def print_status(self):

    """ @brief  Prints status information for this collection. """

    print("Found factors:", len(self.found_factors));
    print("Found primes:", len(self.found_primes));

    found_factors = list(self.found_factors);
    found_factors.sort();

    for i in range(len(found_factors)):
      print(" Factor " + str(i) + ":", found_factors[i]);
    print("");

  def __eq__(self, other):

    """ @brief  Compares this collection to another collection.

        @param other  The other collection to which to compare this collection.

        @return   True if the collections are equal, False otherwise. """

    if self.N != other.N:
      return False;
    if self.residual != other.residual:
      return False;
    if self.found_factors != other.found_factors:
      return False;
    if self.found_primes != other.found_primes:
      return False;

    return True; # Note: Disregards timer differences.

  def __str__(self):

    """ @brief  Returns a string representation of this collection.

        @return   A string representation of this collection. """

    return str(self.found_factors);

  def __repr__(self):

    """ @brief  Returns a string representation of this collection.

        @return   A string representation of this collection. """

    return str(self);