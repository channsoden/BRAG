#!/usr/bin/env python

# Nonstandard Modules
import numpy as np
from scipy.misc import factorial

def composite_poisson_likelihood(counts, likelihoods, times, rates):
    # compute the likelhood landscape for the rate parameter of a poisson process
    # when there is uncertainty about the number of counts observed
    # and uncertainty about the amount of time observed.
    if (len(counts) != len(likelihoods) or
        len(counts) != len(times)):
        raise ValueError("counts, likelihoods, and times must be of equal length")
    counts = np.array(counts)
    likelihoods = np.array(likelihoods)
    times = np.array(times)
    
    # Calculate the expectation = maximum likelihood estimator
    E = np.sum(counts * likelihoods / times)
    
    means = rates * times[:, None] # shape = (len(times), len(rates))
    exponential = np.exp(-means)
    mean_powers = np.array([inline_integer_powers(m, c) for m, c in zip(means, counts)])
    denominators = factorial(counts)
    unscaled = mean_powers * exponential / denominators[:, None]
    scaled = unscaled * likelihoods[:, None]
    composite = np.sum(scaled, axis=0)

    return np.concatenate( (np.array([E]), composite) )

def inline_integer_powers(array, power):
    # np.power is slower for integer powers than inline multiplication.
    # except for squares, which is done differently
    if power == 0:
        return np.ones(len(array))
    elif power == 2:
        return np.power(array, power)
    else:
        expression = 'array*' * power + '1'
        return eval(expression)
