import numpy as np
import pandas as pd


def outliers(data):
    n_data = data.shape[0]
    n_species = data.shape[1]

    # convert from DataFrame to array
    x_in = np.array(data)

    # transform columns by z-score
    x_as = (x_in - np.tile(x_in.mean(axis=0), (n_data, 1))) / np.tile(x_in.std(axis=0, ddof=1), (n_data, 1))

    # singular value decomposition on covariance matrix
    _, b, _ = np.linalg.svd(np.cov(x_as.T))
    n_comp = sum(b > 1)

    # singular value decomposition on transformed data
    U, S, V = np.linalg.svd(x_as)
    x_rec = U[:, :n_comp].dot(np.diag(S)[:n_comp, :n_comp]).dot(V.T[:, :n_comp].T)

    # matrix E
    E = x_as - x_rec
    EE = E.dot(E.T)
    EE = np.diag(EE)

    # some numbers
    perclim = 0.95
    numb = int(max(1,np.floor(perclim * n_data)))
    rest = n_data - numb

    LIM = []
    for i in range(1000):
        order = np.random.permutation(n_data)
        subset = EE[order[:numb]]
        subset = subset[(-subset).argsort()]            # sort in descending order
        LIM.append(subset[rest] + 0.001)                # choose (rest)th element from list
    LIM = np.median(LIM)

    # Outliers
    OL = []
    index = (-EE).argsort()
    d = EE[index]
    for i, value in enumerate(d):
        if value > LIM:
            OL.append(index[i])
    l = len(OL)

    # 2 times above LIMIT
    OL_2times = []
    i = 0
    stop = 0
    while stop == 0 and i <= l:
        if EE[OL[i]] > 2 * LIM:
            OL_2times = OL[:i + 1]
            i += 1
        else:
            stop = 1

    # 10 times above LAST
    OL_new = []
    EE_new = EE[OL[:l]] - np.ones(l) * LIM
    last = EE_new[-1]
    i = l
    stop = 0
    while stop == 0 and i >= 2:
        i -= 1
        if EE_new[i - 1] < 10 * last:
            OL_new = OL[:i - 1]
        else:
            stop = 1

    OL = np.unique(np.concatenate([OL_new, OL_2times]))
    OL.astype(int)
    l = len(OL)

    OL_out = np.zeros((l, 2))
    OL_out[:, 0] = OL.T
    if l > 0:
        for i in range(l):
            E2_OL = E[OL[i], :] ** 2
            index = (-E2_OL).argsort()
            OL_out[i, 1] = index[0]
    OL_out = OL_out.astype(int)

    # prepare output data
    x_out = x_in
    if len(OL_out) > 0:
        for i in range(l):
            x_out[OL_out[i, 0], OL_out[i, 1]] = np.NaN

    # switch back to dataframe
    data_out = pd.DataFrame(x_out, index=data.index, columns=data.columns)
    return data_out, OL_out


