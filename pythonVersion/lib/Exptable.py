import numpy as np

class Exptable():
    # Create an exponential table for the lookup calculation of potential.
    def createExpTable(self, dim, sigma, exp_resolution=1000000):
        delta_dis = np.linspace(0.0, dim, num=exp_resolution)  # Range of distance difference.
        tab = self.updateExpTable(sigma, delta_dis)
        return tab, dim

    def updateExpTable(self, sigma, delta_dis):
        sigma2 = sigma * sigma
        # Try 2 time sigma2
        return np.exp(- delta_dis / (2 * sigma2)) / sigma2