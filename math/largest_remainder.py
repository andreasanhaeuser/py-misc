#!/usr/bin/python3

import numpy as np

def largest_remainder(votes, total_seats):
    """Return largest remainder destribution of seats.

        Parameters
        ----------
        votes : Iterable of float, length N
        total_seats : int > 0

        Returns
        -------
        seats : array of int, length N
    """
    votes = np.array(votes)
    total_votes = np.sum(votes)
    total_seats = int(total_seats)
    assert total_seats > 0

    # allocate automatic seats
    seats_per_vote = total_seats / total_votes
    fractional_seats = votes * seats_per_vote
    automatic_seats = np.array(np.floor(fractional_seats), dtype=int)

    # allocate highest remainder seats
    remaining_seats = int(total_seats - np.sum(automatic_seats))
    remainders = fractional_seats - automatic_seats
    idx_add = np.argsort(remainders)[-remaining_seats:]
    additional_seats = np.zeros_like(automatic_seats, dtype=int)
    additional_seats[idx_add] = 1

    return automatic_seats + additional_seats
