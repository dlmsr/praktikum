# coding=utf-8

from numpy import sum, sqrt, array

def linear_fit(x, y):
    """ Use least squares to fit linear function

    Assumes y = A*x + B

    Parameters
    ----------
    x : N-length sequence
        The independent variable where the data is measured
    y : N-length sequence
        The dependent variable

    Returns
    --------

    popt : array
        Optimal values for the parameters so that the sum of the
        squared errors is minimized.
    perr : array
        The estimated errors of the optimized parameters

    Notes
    -----

    Formulas from 
      TAYLOR, John R.: An Introduction to Error Analysis
    """
    # x and y must be equal in length
    assert(len(x) == len(y))

    N = len(y)
    Delta = N * sum(x**2) - (sum(x))**2

    A = (N * sum(x * y) - sum(x) * sum(y)) / Delta
    B = (sum(x**2) * sum(y) - sum(x) * sum(x * y)) / Delta

    sigma_y = sqrt(sum((y - A * x - B)**2) / (N - 2))

    A_error = sigma_y * sqrt(N / Delta)
    B_error = sigma_y * sqrt(sum(x**2) / Delta)

    # f(x) = A * x + B
    return array([A, A_error, B, B_error])
