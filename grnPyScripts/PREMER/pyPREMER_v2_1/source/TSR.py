import numpy as np
import pandas as pd


class InfoObject:
    def __init__(self, r):
        self.r = r
        self.O = np.where(r == 1)[0]                # indices with existing values
        self.M = np.where(r == 0)[0]                # indices with missing values
        self.nO = sum(r == 1)                       # number of existing values
        self.nM = sum(r == 0)                       # number of missing values


def trimmed_scores_regression(data):
    n_data = data.shape[0]
    n_species = data.shape[1]

    # linear correlation of data
    # (pandas is able to process NaN values)
    B = data.corr()

    # replace remaining NaN value with 0
    B.replace(np.nan, 0, inplace=True)

    # singular value decomposition
    _, s, _ = np.linalg.svd(B)

    # how many components used (90% variance)
    ve = 0                                          # percent variance explained
    n_comp = 0                                      # number of components to use for TSR
    while ve < 0.9:
        ve = ve + s[n_comp] / sum(s)
        n_comp += 1

    # collect info on rows of dataset
    pats = []                                       # list of info objects on each row of dataset
    for i in range(n_data):
        r = np.array(~np.isnan(data.iloc[i, :]))    # binary array with positions of missing values in each row
        pat = InfoObject(r)                         # for details see class above
        pats.append(pat)                            # add info to list

    # prepare output dataset
    x_mean = np.array(data.mean(axis=0))            # array of mean column values
    x_in = np.array(data)                           # switch to arrays (better for indexing)
    x_out = np.array(data)                          # copy original data
    x_missing = np.isnan(x_in)                      # binary array of missing values
    rows, columns = np.where(np.isnan(x_in))        # lists with rows / columns of missing values
    x_out[np.isnan(x_out)] = 0                      # set missing values to 0
    for idx in range(len(rows)):                    # replace missing values with mean of columns values
        x_out[rows[idx], columns[idx]] \
            = x_mean[columns[idx]]

    # some more variables
    max_iterations = 5000
    convergence_criteria = 1e-10
    diff = 100

    # ???
    count = 0
    while count < max_iterations and diff > convergence_criteria:
        count += 1
        data_missing = x_out[x_missing]                     # array of imputed values
        x_mean = x_out.mean(axis=0)                         # array of mean column values
        S = np.cov(x_out.T)                                 # covariance matrix
        x_temp = x_out - np.ones(n_species) * x_mean        # correct each column by mean value

        # singular value decomposition
        if n_data > n_species:
            V, D, U = np.linalg.svd(x_temp.T)               # here python behavior seems to be different from Matlab
        else:
            U, D, V = np.linalg.svd(x_temp)                 # here python behavior seems to be different from Matlab

        V = V[:, :n_comp]                                   # trim matrix V to relevant components

        for i in range(n_data):
            if pats[i].nM > 0:
                L = V[pats[i].O, :min(n_comp, pats[i].nO)]
                S_11 = S[np.ix_(pats[i].O, pats[i].O)]
                S_21 = S[np.ix_(pats[i].M, pats[i].O)]
                z1 = x_temp[i, pats[i].O].T
                part1 = S_21.dot(L)                             # numpy dot product (in python 3.5 can be replaced by @)
                part2 = np.linalg.pinv(L.T.dot(S_11).dot(L))    # numpy dot product (in python 3.5 can be replaced by @)
                part3 = L.T.dot(z1)                             # numpy dot product (in python 3.5 can be replaced by @)
                z2 = part1.dot(part2).dot(part3)                # numpy dot product (in python 3.5 can be replaced by @)
                x_temp[i, pats[i].M] = z2.T
        x_out = x_temp + np.ones((n_data, 1)) * x_mean          # correct missing values
        d = (x_out[x_missing] - data_missing) ** 2              # check convergence
        diff = np.mean(d)                                       # check convergence

    # convert back to dataframe and return
    data_out = pd.DataFrame(x_out, index=data.index, columns=data.columns)
    return data_out
