import numpy as np
import pandas as pd
import os


class ResultsClass():
    """
    :param path: path of results files
    :return: pyPREMER output structure
    """
    def __init__(self, path):
        # get file trunk
        self.__path__ = path
        self.__file_name__ = path + path.rstrip('/').split('/')[-1] + '.csv'
        self.__trunk__ = self.__file_name__.split('/')[-1].split('.')[0]

        # read in data
        self.data = pd.read_csv(self.__file_name__, index_col=0)
        self.n_data = self.data.shape[0]
        self.n_species = self.data.shape[1]

        # read in run info
        info_table = pd.read_csv(self.__path__ + self.__trunk__ + '_runinfo.csv', index_col=0)
        self.taumax = int(info_table.loc['taumax'].value)
        self.ert_crit = int(info_table.loc['ert_crit'].value)

        # con array
        self.con_array = pd.read_csv(self.__path__ + 'out_' + self.__trunk__ + '_con_array.txt', sep='\s+', header=None)
        self.con_array.index = self.data.columns
        self.con_array.columns = self.data.columns

        # prediction (write if doesn't exist)
        if os.path.isfile(self.__path__ + 'out_' + self.__trunk__ + '_prediction.txt'):
            self.prediction = pd.read_csv(self.__path__ + 'out_' + self.__trunk__ + '_prediction.txt', index_col=0)
        else:
            self.prediction = []
            for target, row in self.con_array.iterrows():
                for regulator, item in row.iteritems():
                    self.prediction.append([regulator, target, item])
            self.prediction = pd.DataFrame(self.prediction, columns=['regulator', 'target', 'value'])
            self.prediction.sort_values('value', ascending=False, inplace=True)
            self.prediction.index = range(len(self.prediction))
            self.prediction.to_csv(self.__path__ + 'out_' + self.__trunk__ + '_prediction.txt')

        # cond entr 2
        self.cond_entr_2 = pd.read_csv(self.__path__ + 'out_' + self.__trunk__ + '_cond_entr2.txt', sep='\s+', header=None)
        self.cond_entr_2.index = self.data.columns
        self.cond_entr_2.columns = self.data.columns

        # dist
        self.dist = pd.read_csv(self.__path__ + 'out_' + self.__trunk__ + '_dist.txt', sep='\s+', header=None)
        self.dist.index = self.data.columns
        self.dist.columns = self.data.columns

        # T
        self.T = pd.read_csv(self.__path__ + 'out_' + self.__trunk__ + '_T.txt', sep='\s+', header=None)
        self.T.index = self.data.columns
        self.T.columns = self.data.columns

        # taumin
        self.taumin = pd.read_csv(self.__path__ + 'out_' + self.__trunk__ + '_taumin.txt', sep='\s+', header=None)
        self.taumin.index = self.data.columns
        self.taumin.columns = self.data.columns

        # threshold
        with open(self.__path__ + 'out_' + self.__trunk__ + '_threshold.txt') as fid:
            lines = fid.readlines()
            self.threshold = float(lines[0])

        # H1
        with open(self.__path__ + 'out_' + self.__trunk__ + '_H1.txt') as fid:
            lines = fid.readlines()
            self.H1 = map(float, lines[0].split())

        # mutual information arrays
        MI_matrix = np.loadtxt(self.__path__ + 'out_' + self.__trunk__ + '_MIl_MI_MIm_MIs.txt')
        self.MIl = np.zeros((self.n_species, self.n_species, self.taumax + 1))
        for idx in range(self.taumax + 1):
            self.MIl[:, :, idx] = MI_matrix[idx * self.n_species:(idx + 1) * self.n_species, :]
        offset = self.n_species * (self.taumax + 1)
        self.MI = np.zeros((self.n_species, self.n_species, self.taumax + 1))
        for idx in range(self.taumax + 1):
            self.MI[:, :, idx] = MI_matrix[idx * self.n_species + offset:(idx + 1) * self.n_species + offset, :]
        self.MIm = np.zeros((self.n_species, self.n_species, self.taumax + 1))
        for idx in range(self.taumax + 1):
            self.MIm[:, :, idx] = MI_matrix[idx * self.n_species + 2 * offset:(idx + 1) * self.n_species + 2 * offset, :]
        self.MIs = np.zeros((self.n_species, self.n_species, self.taumax + 1))
        for idx in range(self.taumax + 1):
            self.MIs[:, :, idx] = MI_matrix[idx * self.n_species + 3 * offset:(idx + 1) * self.n_species + 3 * offset, :]

        # H2 array
        H2_matrix = np.loadtxt(self.__path__ + 'out_' + self.__trunk__ + '_H2.txt')
        self.H2 = np.zeros((self.n_species, self.n_species, self.taumax + 1))
        for idx in range(self.taumax + 1):
            self.H2[:, :, idx] = H2_matrix[idx * self.n_species:(idx + 1) * self.n_species, :]

        # Higher order arrays (1)
        if self.ert_crit > 1:

            # H3 arrays
            H3_matrix = np.loadtxt(self.__path__ + 'out_' + self.__trunk__ + '_H3.txt')
            self.H3 = np.zeros((self.n_species, self.n_species, self.n_species))
            for idx in range(self.n_species):
                self.H3[:, :, idx] = H3_matrix[idx * self.n_species:(idx + 1) * self.n_species, :]

            # MI3 arrays
            MI3_matrix = np.loadtxt(self.__path__ + 'out_' + self.__trunk__ + '_MI3.txt')
            self.MI3 = np.zeros((self.n_species, self.n_species, self.n_species))
            for idx in range(self.n_species):
                self.MI3[:, :, idx] = MI3_matrix[idx * self.n_species:(idx + 1) * self.n_species, :]

            # 3D conditional entropy arrays
            CE3_matrix = np.loadtxt(self.__path__ + 'out_' + self.__trunk__ + '_cond_entr3.txt')
            self.cond_entr_3 = np.zeros((self.n_species, self.n_species, self.n_species))
            for idx in range(self.n_species):
                self.cond_entr_3[:, :, idx] = CE3_matrix[idx * self.n_species:(idx + 1) * self.n_species, :]

            # Higher order arrays (2)
            if self.ert_crit > 2:

                # H4 arrays
                H4_matrix = np.loadtxt(self.__path__ + 'out_' + self.__trunk__ + '_H4.txt')
                self.H4 = np.zeros((self.n_species, self.n_species, self.n_species, self.n_species))
                count = 0
                for idx1 in range(self.n_species):
                    for idx2 in range(self.n_species):
                        self.H4[:, :, idx1, idx2] = H4_matrix[count * self.n_species:(count + 1) * self.n_species]
                        count += 1

                # 4D conditional entropy arrays
                CE4_matrix = np.loadtxt(self.__path__ + 'out_' + self.__trunk__ + '_cond_entr4.txt')
                self.cond_entr_4 = np.zeros((self.n_species, self.n_species, self.n_species, self.n_species))
                count = 0
                for idx1 in range(self.n_species):
                    for idx2 in range(self.n_species):
                        self.cond_entr_4[:, :, idx1, idx2] = CE4_matrix[count * self.n_species:(count + 1) * self.n_species]
                        count += 1


