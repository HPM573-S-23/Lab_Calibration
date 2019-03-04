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

        self.cohortIDs = []             # IDs of cohorts to simulate
        self.mortalitySamples = []      # values of mortality probability at which the posterior should be sampled
        self.normalizedWeights = []     # normalized likelihood weights (sums to 1)
        self.mortalityResamples = []  # resampled values for constructing posterior estimate and interval

    def sample_posterior(self, n_samples):
        """ sample the posterior distribution of the mortality probability,
         :param n_samples: number of samples from the posterior distribution
         """

        # specifying the seed of the numpy random number generator
        np.random.seed(1)

        # cohort ids
        self.cohortIDs = range(n_samples)

        # find values of mortality probability at which the posterior should be evaluated
        self.mortalitySamples = np.random.uniform(
            low=CalibSets.POST_L,
            high=CalibSets.POST_U,
            size=CalibSets.POST_N)

        # create a multi cohort
        multi_cohort = SurvivalCls.MultiCohort(
            ids=self.cohortIDs,
            mortality_probs=self.mortalitySamples,
            pop_sizes=[CalibSets.SIM_POP_SIZE]*CalibSets.POST_N
        )

        # simulate the multi cohort
        multi_cohort.simulate(n_time_steps=CalibSets.TIME_STEPS)

        # calculate the likelihood of each simulated cohort
        weights = []
        for cohort_id in self.cohortIDs:

            # get the average survival time for this cohort
            mean = multi_cohort.multiCohortOutcomes.meanSurvivalTimes[cohort_id]

            # construct a gaussian (normal) likelihood
            # with mean calculated from the simulated data and standard deviation from the clinical study.
            # evaluate this pdf (probability density function) at the mean reported in the clinical study.
            weight = stat.norm.pdf(
                x=CalibSets.OBS_MEAN,
                loc=mean,
                scale=CalibSets.OBS_STDEV)

            # store the weight
            weights.append(weight)

        # normalize the likelihood weights
        sum_weights = np.sum(weights)
        self.normalizedWeights = np.divide(weights, sum_weights)

        # produce the list to report the results
        csv_rows = \
            [['Cohort ID', 'Likelihood Weights', 'Mortality Prob']]  # list containing the calibration results
        for i in range(len(self.mortalitySamples)):
            csv_rows.append(
                [self.cohortIDs[i], self.normalizedWeights[i], self.mortalitySamples[i]])

        # write the calibration result into a csv file
        InOutSupport.write_csv(
            file_name='CalibrationResults.csv',
            rows=csv_rows)

        # re-sample mortality probability (with replacement) according to likelihood weights
        self.mortalityResamples = np.random.choice(
            a=self.mortalitySamples,
            size=n_samples,
            replace=True,
            p=self.normalizedWeights)

    def get_mortality_estimate_credible_interval(self, alpha):
        """
        :param n_samples: number of resamples from parameter values
        :param alpha: the significance level
        :returns tuple (mean, [lower, upper]) of the posterior distribution"""

        # calculate the credible interval
        sum_stat = StatSupport.SummaryStat(name='Posterior samples',
                                           data=self.mortalityResamples)

        estimate = sum_stat.get_mean()  # estimated mortality probability
        credible_interval = sum_stat.get_PI(alpha=alpha)  # credible interval

        return estimate, credible_interval

