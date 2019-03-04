from enum import Enum
import scipy.stats as stat
import numpy as np
import SimPy.InOutFunctions as InOutSupport
import SimPy.StatisticalClasses as StatSupport
import MultiSurvivalModelClasses as SurvivalCls
import CalibrationSettings as CalibSets


class CalibrationColIndex(Enum):
    """ indices of columns in the calibration results cvs file  """
    ID = 0          # cohort ID
    W = 1           # likelihood weight
    MORT_PROB = 2   # mortality probability


class Calibration:
    def __init__(self):
        """ initializes the calibration object"""

        self.cohortIDs = range(CalibSets.POST_N)   # IDs of cohorts to simulate
        self.mortalitySamples = []      # values of mortality probability at which the posterior should be sampled
        self.normalizedWeights = []     # normalized likelihood weights (sums to 1)
        self.csvRows = \
            [['Cohort ID', 'Likelihood Weights', 'Mortality Prob']]  # list containing the calibration results
        self.mortalityResamples = []  # resampled values for constructing posterior estimate and interval

    def sample_posterior(self):
        """ sample the posterior distribution of the mortality probability """

        # specifying the seed of the numpy random number generator

        # find values of mortality probability at which the posterior should be evaluated

        # create a multi cohort

        # simulate the multi cohort

        # calculate the likelihood of each simulated cohort

            # get the average survival time for this cohort

            # construct a gaussian (normal) likelihood
            # with mean calculated from the simulated data and standard deviation from the clinical study.
            # evaluate this pdf (probability density function) at the mean reported in the clinical study.

            # store the weight

        # normalize the likelihood weights

        # re-sample mortality probability (with replacement) according to likelihood weights

        # produce the list to report the results

        # write the calibration result into a csv file

    def get_mortality_estimate_credible_interval(self, alpha):
        """
        :param alpha: the significance level
        :returns tuple (mean, [lower, upper]) of the posterior distribution"""

        # calculate the credible interval

        # estimated mortality probability

        # credible interval

