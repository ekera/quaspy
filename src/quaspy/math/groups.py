""" @brief  A module for performing or simulating arithmetic operations in
            cyclic groups. """

from abc import abstractmethod;

from gmpy2 import mpz;

from gmpy2 import gcd;
from gmpy2 import powmod;

from gmpy2 import invert as mpz_inv;

class CyclicGroupElement:

  """ @brief  An abstract class that represents a group element in a finite
              cyclic group, that is written multiplicatively. """

  @abstractmethod
  def __pow__(self, e : int):

    """ @brief  Returns the group element g^e, for g this group element.

        @param e  The exponent e.

        @return   The group element g^e, for g this group element. """

    pass;

  @abstractmethod
  def __mul__(self, x):

    """ @brief  Returns the group element g * x, for g this group element.

        @param x  The group element x.

        @return   The group element g * x, for g this group element. """

    pass;

  @abstractmethod
  def __eq__(self, x):

    """ @brief  Compares this group element to another group element.

        The two elements must be in the same group.

        @param x  The group element to which to compare this group element, or
                  one (1) to compare this group element to the identity element.

        @return   True if x is equal to this group element, False otherwise. """

    pass;

  @abstractmethod
  def __hash__(self):

    """ @brief  Returns the hash digest of this group element.

        @return   The hash digest of this group element. """

    pass;


class IntegerModRingMulSubgroupElement(CyclicGroupElement):

  """ @brief  A class that represents the group element g in Z_N^*. """

  def __init__(self, g, N):

    """ @brief  Constructs the group element g in Z_N^*.

        @param g  The element g. An integer in [0, N), coprime to N.

        @param N  The modulus N. A positive integer. """

    self.g = mpz(g);
    self.N = mpz(N);

    if not ((1 <= self.g < self.N) and (1 == gcd(self.g, self.N))):
      raise Exception("Error: Incorrect parameters.");

  def __pow__(self, e : int):

    """ @brief  Returns the group element g^e, for g this group element.

        @param e  The exponent e.

        @return   The group element g^e, for g this group element. """

    if e < 0:
      x = mpz_inv(self.g, self.N);
      x = powmod(x, mpz(-e), self.N);
    else:
      x = powmod(self.g, mpz(e), self.N);

    return IntegerModRingMulSubgroupElement(x, self.N);

  def __mul__(self, x):

    """ @brief  Returns the group element g * x, for g this group element.

        @param x  The group element x.

        @return   The group element g * x, for g this group element. """

    return IntegerModRingMulSubgroupElement((self.g * x.g) % self.N, self.N);

  def __eq__(self, x):

    """ @brief  Compares this group element to another group element.

        @param x  The group element to which to compare this group element, or
                  one (1) to compare this group element to the identity element.

        @return   True if x is equal to this group element, False otherwise. """

    if (x == 1) and (self.g == 1):
      return True;

    if type(x) != type(self):
      return False;

    if x.N != self.N:
      return False;

    return self.g == x.g;

  def __hash__(self):

    """ @brief  Returns the hash digest of this group element.

        @return   The hash digest of this group element. """

    return hash((self.g, self.N));

  def __str__(self):

    """ @brief  Returns a string representation of the group element.

        @return   A string representation of the group element. """

    return str(self.g) + " (mod " + str(self.N) + ")";

  def __repr__(self):

    """ @brief  Returns a string representation of the group element.

        @return   A string representation of the group element. """

    return str(self);


class ShortWeierstrassCurveOverPrimeField:

  """ @brief  A class that represents an elliptic curve on short Weierstrass
              form y^2 = x^3 + ax + b (mod p) for p a prime number. """

  def __init__(self, a, b, p):

    """ @brief  Constructs an elliptic curve on short Weierstrass form
                y^2 = x^3 + ax + b (mod p) for p a prime number.

        @param a  The Weierstrass coefficient a.

        @param b  The Weierstrass coefficient b.

        @param p  The prime p. """

    self.p = mpz(p);

    self.a = mpz(a) % self.p;
    self.b = mpz(b) % self.p;

  def identity_y(self):

    """ @brief  Returns the y coordinate of the point used to represent the
                point at infinity on this curve.

        The point at infinity is represented by (0, self.identity_y()), where
        the y coordinate is 0 if (0, 0) is not in E, and (0, 1) otherwise.

        @return   The y coordinate of the point used to represent the
                  point at infinity on this curve."""

    if 0 == self.b:
      return 1; # Point at infinity is represented as (0, 1).
    else:
      return 0; # Point at infinity is represented as (0, 0).

  def __eq__(self, E):

    """ @brief  Compares this curve to another curve.

        @param E  The curve to which to compare this curve.

        @return   True if E is equal to this curve, False otherwise. """

    return (self.p == E.p) and (self.a == E.a) and (self.b == E.b);

  def __str__(self):

    """ @brief  Returns a string representation of the curve.

        @return   A string representation of the curve. """

    return "y^2 = x^3 + " + str(self.a) + " x + " + str(self.b) + \
      " (mod " + str(self.p) + ")";

  def __repr__(self):

    """ @brief  Returns a string representation of the curve.

        @return   A string representation of the curve. """

    return str(self);


class PointOnShortWeierstrassCurveOverPrimeField(CyclicGroupElement):

  """ @brief  A class that represents a point (x, y) on a short Weierstrass
              curve over a prime field.

      Note that we write elliptic curve groups multiplicatively: This is
      non-standard, but beneficial since it allows elliptic curve groups to
      be substituted for multiplicative groups, and vice versa, in function
      calls. More specifically, P + Q, [e] P and [-1] P becomes P * Q, P^e
      and P^-1, respectively, when writing the group multiplicatively. """

  def __init__(self, x, y, E):

    """ @brief  Constructs the point (x, y) on the elliptic curve E.

        @param x  The x coordinate.

        @param y  The y coordinate.

        @param E  The elliptic curve E. """

    self.x = mpz(x) % E.p;
    self.y = mpz(y) % E.p;
    self.E = E;

    if not self.is_identity():
      if 0 != (((self.x ** 3) + E.a * self.x + E.b - (self.y ** 2)) % E.p):
        raise Exception("Error: The point is not on the curve.");

  def is_identity(self):

    """ @brief  Returns True if this point is the identity, False otherwise.

        @return   True if this point is the identity, False otherwise. """

    return (0 == self.x) and (self.E.identity_y() == self.y);

  def invert(self):

    """ @brief  Returns the point P^-1, for P this point.

        @return   The point P^-1, for P this point. """

    return PointOnShortWeierstrassCurveOverPrimeField(self.x, -self.y, self.E);

  def __pow__(self, e : int):

    """ @brief  Returns the point P^e, for P this point.

        @param e  The exponent e.

        @return   The point P^e, for P this point. """

    # Set P to this point.
    P = PointOnShortWeierstrassCurveOverPrimeField(self.x, self.y, self.E);

    if e < 0:
      P = P.invert();
      e = -e;

    # Set R to the point at infinity on E.
    R = PointOnShortWeierstrassCurveOverPrimeField(
          0, self.E.identity_y(), self.E);

    while 0 != e:
      if 1 == (e % 2):
        R = R * P;
      P = P * P;
      e = e // 2;

    return R;

  def __mul__(self, Q):

    """ @brief  Returns the point P * Q, for P this point.

        @param Q  The point Q.

        @return   The point P * Q, for P this point. """

    if Q == 1:
      return self;

    if self == 1:
      return Q;

    if self == Q.invert():
      if 0 == self.E.b:
        return PointOnShortWeierstrassCurveOverPrimeField(0, 1, self.E);
      else:
        return PointOnShortWeierstrassCurveOverPrimeField(0, 0, self.E);

    if self.E != Q.E:
      raise Exception("Error: The two points are not on the same curve.");

    if self == Q:
      s = (3 * (self.x ** 2) + self.E.a) * mpz_inv(2 * self.y, self.E.p);
      s = s % self.E.p;
      x = ((s ** 2) - 2 * self.x) % self.E.p;
      y = (-s * (x - self.x) - self.y) % self.E.p;
    else:
      s = ((self.y - Q.y) * mpz_inv(self.x - Q.x, self.E.p)) % self.E.p;
      x = ((s ** 2) - self.x - Q.x) % self.E.p;
      y = ( s * (self.x - x) - self.y) % self.E.p;

    return PointOnShortWeierstrassCurveOverPrimeField(x, y, self.E);

  def __eq__(self, Q):

    """ @brief  Compares this point to another point Q.

        @param Q  The point to which to compare this point, or one (1) to
                  compare this point to the identity element.

        @return   True if Q is equal to this point, False otherwise. """

    if (Q == 1) and (self.x == 0) and (self.y != self.E.b):
      return True;

    if type(Q) != type(self):
      return False;

    if Q.E != self.E:
      return False;

    return (Q.x == self.x) and (Q.y == self.y);

  def __hash__(self):

    """ @brief  Returns the hash digest of this pointgroup element.

        @return   The hash digest of this group element. """

    return hash((self.x, self.y, self.E.p, self.E.a, self.E.b));

  def __str__(self):

    """ @brief  Returns a string representation of the group element.

        @return   A string representation of the group element. """

    if self.is_identity():
      return "Point at infinity on " + str(self.E);

    return "(" + str(self.x) + ", " + str(self.y) + ") on " + str(self.E);

  def __repr__(self):

    """ @brief  Returns a string representation of the group element.

        @return   A string representation of the group element. """

    return str(self);


class SimulatedCyclicGroupElement(CyclicGroupElement):

  """ @brief  A class that represents the simulated group element G^d, for G a
              generator of a cylic group of order r. """

  def __init__(self, r, d = 1):

    """ @brief  Constructs the simulated group element G^d, for G a generator of
                a cylic group of order r.

        @param r  The order r. A positive integer.

        @param d  The exponent d on [0, r). """

    self.r = mpz(r);
    self.d = mpz(d);

    if not (0 <= self.d < self.r):
      raise Exception("Error: Incorrect parameters.");

  def __pow__(self, e : int):

    """ @brief  Returns the group element g^e, for g this group element.

        @param e  The exponent e.

        @return   The group element g^e, for g this group element. """

    return SimulatedCyclicGroupElement(self.r, (self.d * e) % self.r);

  def __mul__(self, x):

    """ @brief  Returns the group element g * x, for g this group element.

        @param x  The group element x.

        @return   The group element g * x, for g this group element. """

    if x == 1:
      return SimulatedCyclicGroupElement(self.r, self.d);

    if type(self) != type(x):
      raise Exception("Error: The two elements must be in the same group.");

    if self.r != x.r:
      raise Exception("Error: The two elements must be in the same group.");

    return SimulatedCyclicGroupElement(self.r, (self.d + x.d) % self.r);

  def __eq__(self, x):

    """ @brief  Compares this group element to another group element.

        @param x  The group element to which to compare this group element, or
                  one to compare this group element to the identity element.

        @return   True if x is equal to this group element, False otherwise. """

    if (x == 1) and (self.d == 0):
      return True;

    if type(x) != type(self):
      return False;

    if x.r != self.r:
      return False;

    return self.d == x.d;

  def __hash__(self):

    """ @brief  Returns the hash digest of this group element.

        @return   The hash digest of this group element. """

    return hash((self.d, self.r));

  def __str__(self):

    """ @brief  Returns a string representation of the group element.

        @return   A string representation of the group element. """

    return "C(" + str(self.d) + " / " + str(self.r) + ")";

  def __repr__(self):

    """ @brief  Returns a string representation of the group element.

        @return   A string representation of the group element. """

    return str(self);