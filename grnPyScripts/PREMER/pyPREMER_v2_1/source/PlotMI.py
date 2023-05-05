import numpy as np
import pandas as pd
import os
from scipy.stats import zscore
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt


def plot_mi(results, show=0):
    if results.n_species < 50:
        for idx in range(results.n_species):
            fig = plt.figure(idx)
            miarray = results.MIl[idx, :, :]
            miarrayneg = results.MIl[:, idx, :]
            miarraynegorder = np.fliplr(miarrayneg)
            miarraycomplete = np.concatenate((miarraynegorder[:, :-1], miarray), axis=1)
            plt.imshow(miarraycomplete, interpolation='nearest', vmin=0, vmax=1)
            plt.xlabel('Time lags')
            plt.xticks(range(results.taumax * 2 + 1), range(-results.taumax, results.taumax + 1))
            plt.ylabel('Variables')
            plt.yticks(range(results.n_species), results.data.columns)
            cb = plt.colorbar()
            cb.set_label('Normalized Mutual Information')
            plt.title('MI between ' + results.data.columns[idx] + ' and others')
            plt.tight_layout()

            if show == 0:
                plt.show()
            else:
                if not os.path.exists(results.__path__ + 'plots/'):
                    os.mkdir(results.__path__ + 'plots/')
                plt.savefig(results.__path__ + 'plots/' + results.__trunk__ + '_' + results.data.columns[idx] + '_MI.png')
                plt.close()
    else:
        print('Too many plots to plot')


def plot_graph(results, show=0):
    scaled_centered_dist = zscore(results.dist, ddof=1)
    dissimilarities = pdist(scaled_centered_dist, 'euclidean')

    if dissimilarities.min() < 0.5:
        dissimilarities = dissimilarities + dissimilarities.mean()
        print('Some variables are very close. '
              'The distances between variables will be increased to improve visualization.')

    D = squareform(dissimilarities)

    # calculate eigenvalues
    W = np.eye(results.n_species) - np.ones((results.n_species, results.n_species)) * 1.0 / float(results.n_species)
    M = W.dot(-0.5 * D * D).dot(W)
    E, V = np.linalg.eig((M + M.T) / 2.0)

    # sort eigenvalues in descending order
    fi = (-E).argsort()
    fe = E[fi]
    
    # Use only positive eigenvalues
    pe = fe > 0
    Y = -V[:, fi[pe]].dot(np.diag(np.sqrt(fe[pe])))
    Y = pd.DataFrame(Y, index=results.con_array.columns)

    plt.scatter(Y[0], Y[1], color='red', s=50)
    for idx, row in Y.iterrows():
        plt.text(row[0], row[1], idx)
    # plot arrows
    for var1 in results.con_array:
        for var2 in results.con_array:
            if results.con_array.loc[var1, var2] > 0:
                if results.T.loc[var1, var2] > 0:
                    arrow(Y.ix[var1][0], Y.ix[var1][1], Y.ix[var2][0], Y.ix[var2][1])
                elif results.T.loc[var2, var1] > 0:
                    arrow(Y.ix[var2][0], Y.ix[var2][1], Y.ix[var1][0], Y.ix[var1][1])
    plt.xlim(min(Y[0]) - 0.5, max(Y[0]) + 0.5)
    plt.ylim(min(Y[1]) - 0.5, max(Y[1]) + 0.5)
    plt.xticks([], [])
    plt.yticks([], [])

    if show == 0:
        plt.show()
    else:
        if not os.path.exists(results.__path__ + 'plots/'):
            os.mkdir(results.__path__ + 'plots/')
        plt.savefig(results.__path__ + 'plots/' + results.__trunk__ + '_networkgraph.png')
        plt.close()


def arrow(x1, y1, x2, y2):
    # get length of arrow
    l1 = x2 - x1
    l2 = y2 - y1
    # plot
    plt.arrow(x1 + 0.1 * l1, y1 + 0.1 * l2, 0.8 * l1, 0.8 * l2, head_width=0.1, head_length=0.1, fc='k', ec='k')