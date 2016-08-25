import random as pyrandom
import numpy as np

class DataGen():

    def stddata(self, data, dim):
        for i in range(dim):
            # S1 Standardize data.
            data[:, i] = (data[:, i] - np.mean(data[:, i])) / np.std(data[:, 1])
            # limit range to -.5 ~ .5 in each dimension
            data[:, i] = data[:, i] / np.max(np.absolute(data[:, i]))
            data[:, i] = data[:, i] / 2
        return data

    # Generate data set based of dimension and num_of_cluster.
    def datagen(self, dim, c, sigma=0.2, minnr=50, maxnr=200):
        for i in range(c):
            nr = pyrandom.randrange(minnr, maxnr, 1)
            meanvec = np.random.rand(dim)
            covmat = sigma ** 2 * np.cov(np.random.rand(dim, dim))
            dtmp = np.random.multivariate_normal(meanvec, covmat, nr)
            if (i == 0):
                result = dtmp.copy()
            else:
                result = np.vstack((result, dtmp.copy()))
        return self.stddata(result, dim)