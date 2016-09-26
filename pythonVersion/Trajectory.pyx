cimport numpy as np
import numpy as np
from libc.stdlib cimport rand, malloc, free
from libc.math cimport exp


def potential_ds(np.ndarray[np.float64_t, ndim = 2] data, np.ndarray[np.float64_t, ndim = 1] grid, double sigma = 0.2):
    cdef int N, dim
    N, dim = data.shape[0], data.shape[1]
    cdef double potential, p_sum = 0.0
    cdef int j, i
    cdef double sigma2 = sigma * sigma
    for j in range(N):
        for i in range(dim):
            p_sum += (grid[i] - data[j, i]) * (grid[i] - data[j, i])
        potential += -exp(-0.5 * p_sum / sigma2)
        p_sum = 0
    return potential


def PTSM(np.ndarray[np.float64_t, ndim = 1] pos, np.ndarray[np.float64_t, ndim = 2] X, \
    np.ndarray[np.float64_t, ndim = 1] vel, np.ndarray[np.float64_t, ndim = 1] lookupExp,\
    int table_res, double norm_max, double sigma = 0.25, double dt = 0.01, double r = 0.99,\
    int Nsamp = 5000, double compensation = 0.01):

    cdef int N, dim
    N, dim = X.shape[0], X.shape[1]
    cdef int i, j, steps, lookupIdx  # These 3 are for iterations.
    cdef double sigma2, m
    cdef double d, V, v_sum = 0.0
    cdef double force1p_sum = 0.0
    cdef double vel_sum = 0.0
    cdef double * force = < double * > malloc(dim * sizeof(double))
    cdef double * trj = < double * > malloc(Nsamp * (dim + 1) * sizeof(double))
    cdef double * velocity = < double * > malloc(dim * sizeof(double))
    cdef double * position = < double * > malloc(dim * sizeof(double))
    cdef double * temp = < double * > malloc(dim * sizeof(double))
    cdef double * viewForce = < double * > malloc(Nsamp * dim * sizeof(double))

    for i in range(dim):
        velocity[i] = vel[i]
        position[i] = pos[i]

    sigma2 = sigma * sigma
    m = compensation / sigma2
    # --------------------
    # Force
    for steps in range(Nsamp):
        for i in range(dim):
            force[i] = 0
        for j in range(N):
            for i in range(dim):
                temp[i] = (position[i] - X[j, i])
                force1p_sum += temp[i] * temp[i]
            # Look for exp index

            lookupIdx = int(force1p_sum / norm_max * table_res)


            if (lookupIdx >= table_res):
                lookupIdx = table_res - 1
            for i in range(dim):
                #                 print lookupExp[lookupIdx]
                force[i] += -temp[i] * lookupExp[lookupIdx]

            force1p_sum = 0
            # ------------
            #         Now update pos and vel
            #         This is the main part to get the trajectory information in terms of the new velocity
            #         and position

        for i in range(dim):
            viewForce[steps * dim + i] = force[i]
            # Should be dt /m * force[i], but since m = 1.0
            velocity[i] = r * velocity[i] + dt * force[i] / m
            position[i] = position[i] + dt * velocity[i]
            vel_sum += velocity[i]

            # -----------------------
            #         # Put velocity and new position into the trj array.
        for i in range(dim + 1):
            if i == 0:  # The first column is always given to the velocity
                trj[steps * (dim + 1) + i] = vel_sum
            else:
                trj[steps * (dim + 1) + i] = position[i - 1]
        vel_sum = 0
    # # put trj into a numpy array. Because you can't return a c-array.
    resultTrj = np.zeros(Nsamp * (dim + 1), dtype=np.float64)
    resultVel = np.zeros(dim)
    resultForce = np.zeros(Nsamp * dim)
    # Move c array to numpy array for return
    for i in range(Nsamp * (dim + 1)):
        resultTrj[i] = trj[i]
    for i in range(dim):
        resultVel[i] = velocity[i] * velocity[i]

    for i in range(dim * Nsamp):
        resultForce[i] = viewForce[i]
    # Transpose matrix
    resultTrj = np.reshape(resultTrj, (-1, dim + 1))  # Correct the matrix shape.
    #     # Free up memories.
    #     print testing
    free(trj)
    free(force)
    free(velocity)
    free(position)
    free(viewForce)
    return resultTrj, resultVel, resultForce