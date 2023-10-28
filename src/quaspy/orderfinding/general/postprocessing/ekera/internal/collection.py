""" @brief  A module for collecting candidates for r_tilde. """

class CandidateCollection:

  """ @brief  A class representing a collection of candidates for r_tilde, that
              is kept reduced in the sense that if r_tilde is in the collection,
              then so is c * r_tilde, for c any positive integer.

      The idea is to only add candidates for r_tilde that meet the requirement
      that g^(e * r_tilde) = 1 to the collection, for e a product of cm-smooth
      prime powers. If r_tilde meets this requirement, then of course so does
      c * r_tilde, so it makes sense to keep the collection reduced.

      In the internal set of candidates that represent the collection, only the
      minimum set of candidates for r_tilde required to represent the collection
      are stored. The collection furthermore has features for testing if a
      candidate is already in the collection. This is useful when searching for
      r_tilde, as it reduces the number of exponentiations required, and the
      amount of storage needed to represent candidates, and so forth. """

  def __init__(self):

    """ @brief  Initializes the collection of candidates for r_tilde. """

    self.candidates = set();

  def add(self, candidate):

    """ @brief  Adds a candidate for r_tilde to this collection.

        @param candidate  The candidate for r_tilde.

        @return True if the collection was modified, False otherwise. """

    if candidate in self:
      # The candidate is already in the collection.
      return False;

    for x in self.candidates.copy():
      if (x % candidate) == 0:
        self.candidates.remove(x);

    self.candidates.add(candidate);

    return True;

  def __contains__(self, candidate):

    """ @brief  Checks if a candidate for r_tilde is in this collection.

        @param candidate  The candidate for r_tilde.

        @return True if the candidate is in the collection, False otherwise. """

    if candidate in self.candidates:
      return True;

    for x in self.candidates:
      if (candidate % x) == 0:
        return True;

    return False;

  def __len__(self):

    """ @brief  Returns the number of candidates for r_tilde that represent this
                collection.

        @remark   Note that the collection is kept reduced, so the number of
                  candidates for r_tilde that represent the collection may be
                  fewer than the number of candidates added to the collection.

        @return   The number of candidates for r_tilde that represent the
                  collection. """

    return len(self.candidates);

  def __iter__(self):

    """ @brief  Returns an iterator for the candidates for r_tilde that
                represent this collection.

        @remark   Note that the collection is kept reduced, so the number of
                  candidates for r_tilde that represent the collection may be
                  fewer than the number of candidates added to the collection.

        @return   An iterator for the candidates for r_tilde that represent the
                  collection. """

    return self.candidates.__iter__();

  def __str__(self):

    """ @brief  Returns a string representation of this collection.

        @return   A string representation of this collection. """

    return str(self.candidates);

  def __repr__(self):

    """ @brief  Returns a string representation of this collection.

        @return   A string representation of this collection. """

    return str(self);