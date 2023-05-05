import shutil
import numpy as np
import pandas as pd
import argparse
import os
import sys
import subprocess
from source.GetResults import ResultsClass
from source.Outliers import outliers
from source.PlotMI import plot_mi, plot_graph
from source.TSR import trimmed_scores_regression

description = """
example usage:
  run pyPREMER ./data/b1_glycolysis.csv
"""


epilog = """
\n\noutput structure:
 pyPREMER produces a structure called \'Output\' containing the fields specified below. The Output
 structure can also be imported by calling the ResultsClass with the name of the results folder
 as argument (eg Output = ResutsClass('./results/example_folder/')). File names in the results
 folder should not be changed for this to work properly.

  -Output.data:         input data as DataFrame
  -Output.n_data:       number of datapoints
  -Output.n_species:    number of variables
  -Output.ert_crit:     Number of entropy reduction rounds carried out
  -Output.taumax:       maximum time lag between X and Y

  -Output.H1:           n-vector of entropies.
  -Output.H2:           n*n*(nlags+1) array, joint entropy of 2 variables.
  -Output.H3:           n*n*n array of joint entropy of 3 variables (calculated only if ert_crit >= 2).
  -Output.H4:           n*n*n*n array of joint entropy of 4 variables (calculated only if ert_crit >= 3).

  -Output.MI:           n*n*(nlags+1) array, mutual information (several lags).
  -Output.MIl:          mutual information normalized as in Linfoot (1957).
  -Output.MIm:          mutual information normalized as in Michaels et al. (1998).
  -Output.MIs:          mutual information normalized as in Studholme et al. (1999).
  -Output.MI3:          n*n*n array of three-way mutual information (calculated only if ert_crit >= 2).

  -Output.cond_entr2:   n*n array of conditional entropies of 2 variables.
  -Output.cond_entr3:   n*n*n array of conditional entropies of 3 variables (only if ert crit >= 2).
  -Output.cond_entr4:   n*n*n*n array of conditional entropies of 4 variables (only if ert crit >= 3).

  -Output.dist:         n*n DataFrame of distance between variables.
  -Output.taumin:       n*n DataFrame of the time lags that minimize the entropic distance.
  -Output.con_array:    n*n DataFrame of connections between variables
  -Output.T:            n*n DataFrame of transfer entropies
  -Output.prediction:   prediction table sorted by value in con_array

  -Output.threshold:    adaptive threshold value
  -Output.Y:            coordinates of the points from multidimensional scaling

citation:
 If you use pyPREMER for your work please use the following citation:
  Villaverde A.F., Becker K., Banga J.R. (2016) PREMER: Parallel Reverse Engineering of Biological
  Networks with Information Theory. In: Bartocci E., Lio P., Paoletti N. (eds) Computational Methods 
  in Systems Biology. CMSB 2016. Lecture Notes in Computer Science, vol 9859. Springer, Cham
"""

#################################################################### command line arguments

ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                             description=description,
                             epilog=epilog)
ap.add_argument('data_file', type=str, nargs='?',
                help='Name of data file. Data must be provided as comma separated text file (csv). Column headers '
                     'should contain name of variables. Row names can refer to experimental condition.')
ap.add_argument('-c', type=str,
                help='Name of connection matrix file. Must be provided as comma separated text file (csv). Column '
                     'headers should contain name of variable as input. Row names should contain name of variables as '
                     'target. A value of 1 indicates that an interaction between input and target is possible, while '
                     'a value of 0 excludes this interaction from the inferred set of interactions. '
                     'If conn matrix is not provided all interactions are considered possible.')
ap.add_argument('-q', type=float, default=1,
                help='Entropic parameter: q = 1 (Boltzmann-Gibbs entropy) | q > 1 (Tsallis entropy) (default = 1). '
                     'Value of the entropic parameter. Choose q = 1 for the classic Shannon entropy (also known as '
                     'Boltzmann-Gibbs), or q > 1 for the generalized Tsallis entropy (Tsallis1988). Hint: in case '
                     'you want to try Tsallis entropy, typical values are 1.5 < q < 3.5.')
ap.add_argument('-MItype', default='MI', type=str,
                help='Normalization of mutual information (MI | MImichaels | MIlinfoot | MIstudholme ), (default=MI). '
                     'Selects the type of normalization of mutual information used to create the distance map. Choose'
                     '\'MI\' for the classic, not normalized value; \'MImichaels\' for the normalization presented in'
                     'Michaels1998; \'MIlinfoot\' for the one in Linfoot1957; or \'MIstudholme\' for the one in '
                     'Studholme1999. Note that PREMER always calculates (and outputs) all the normalizations; this '
                     'option is to select the one used in the distance map.')
ap.add_argument('-taumax', type=int, default=10,
                help='Maximum time lag between X and Y considered in the calculation of I(X,Y). (default=10)')
ap.add_argument('-ert_crit', type=int, default=2,
                help='Number of entropy reduction rounds that will be carried out (0, 1, 2, or 3). (default= 2)')
ap.add_argument('-threshold', type=float, default=1,
                help='Entropy reduction threshold. Enter a number between 0.0 and 0.2 to fix it manually, or 1 to '
                     'use a value obtained from the data. (default=1 (auto))')
ap.add_argument('-numthreads', type=int, default=1,
                help='Number of parallel threads. (default=1)')
ap.add_argument('-plotMI', type=int, default=1,
                help='Show plots mutual information arrays (=0). Save plots mutual information arrays (=1). '
                     'Skip plotting (=2). (default=1)')
ap.add_argument('-useStatistics', action='store_false',
                help='Set to 0 if you wish not to impute missing values. (default = 1)')
ap.add_argument('-correctOutliers', action='store_false',
                help='Set to 0 if you wish not to correct outliers. (default = 1)')
ap.add_argument('-listOutput', action='store_true',
                help='List properties of Output structure')
Args = ap.parse_args()


################################################################### read in data and error checking

# read in data as csv file
assert Args.data_file is not None, 'please provide data file'
data = pd.read_csv(Args.data_file, index_col=0)
n_data = data.shape[0]
n_species = data.shape[1]

# read in conn_matrix, if not provided generate
if not Args.c:
    conn_matrix = pd.DataFrame(np.ones((n_species, n_species)), index=data.columns, columns=data.columns).astype(int)
else:
    conn_matrix = pd.read_csv(Args.conn_matrix, index_col=0)
    assert conn_matrix.shape == (n_species, n_species), \
        'conn_matrix does nat have correct shape (%s, %s)' % (n_species, n_species)

# entropic parameter
assert Args.q >= 1, 'Entropic parameter must be >= 1'
# MItype
assert Args.MItype in ['MI', 'MImichaels', 'MIlinfoot', 'MIstudholme'], \
    'MItype must be one of MI | MImichaels | MIlinfoot | MIstudholme'
# taumax
assert 0 < Args.taumax <= n_data / 2 - 1, 'taumax must be within the range n_data / 2 - 1 (%s)' % (n_data / 2 - 1)
# ert_crit
assert 0 < Args.ert_crit <= 3, 'Number of entropy reduction steps must be 1, 2, or 3'
# Entropy reduction threshold
assert (0.0 <= Args.threshold <= 0.2) or Args.threshold == 1, \
    'Entropy reduction threshold must either be 0.0 < thr < 0.2 or thr == 1'
# Entropy reduction threshold
assert 1 <= Args.numthreads, 'Number of threads must be at least 1'
# plotMI
assert 0 <= Args.plotMI <= 3, 'plotMI must be of 0, 1, 2'

if Args.listOutput:
    print(epilog)

################################################################### data pre-processing

if Args.useStatistics:
    # check if any NaN value exist in data
    if data.isnull().values.any():
        # if yes impute missing values
        data = trimmed_scores_regression(data)

    # correct outliers
    if not Args.correctOutliers:
        data, OL = outliers(data)
        # if outliers removed, impute missing values
        if len(OL) > 0:
            print(len(OL), ' Outliers detected in dataset')
            print(OL)
            data = trimmed_scores_regression(data)
else:
    print('You have chosen not to use statistics. '
          'Therefore outliers will not be detected and missing values not imputed.')


#################################################################### prepare system command and run PREMER executable

# create temporary files (data, conmat)
temp_files = []
temp_file_data = Args.data_file.split('/')[-1].split('.')[0]
temp_file_conmat = 'conmat.txt'
data.to_csv(temp_file_data, sep='\t', index=False, header=None)
conn_matrix.to_csv(temp_file_conmat, sep='\t', index=False, header=None)
temp_files.append(temp_file_data)
temp_files.append(temp_file_conmat)
# run system command
#if sys.platform in ['linux2']:
executable = './PREMERlin64.out'
os.system('chmod +x PREMERlin64.out')
#elif sys.platform in ['win32', 'cygwin']:
#executable = './PREMERx64.exe'
#elif sys.platform in ['darwin', 'os2', 'os2emx']:
#executable = './PREMERosx.out'
#os.system('chmod +x PREMERosx.out')
#else:
#    print(sys.platform)
#    print('unable to detect OS automatically')
#    sys.exit()


subprocess.call([str(executable),
                 str(temp_file_data),
                 str(n_species),
                 str(Args.taumax),
                 str(Args.ert_crit),
                 str(Args.q),
                 str(Args.threshold),
                 str(Args.MItype),
                 str(Args.numthreads)])

    

#################################################################### move generated files to designated sub-folder

# create output folder if needed
file_name = Args.data_file.split('/')[-1].split('.')[0]
folder = './results/' + file_name + '/'
if not os.path.exists(folder):
    os.mkdir(folder)

# copy files and remove old files
shutil.copy(Args.data_file, folder + Args.data_file.split('/')[-1])
shutil.copy(file_name, folder + file_name + '.txt')
shutil.copy('conmat.txt', folder + file_name + '_conmat.txt')
shutil.copy('./results/out_' + file_name + '_H1.txt', folder + 'out_' + file_name + '_H1.txt')
shutil.copy('./results/out_' + file_name + '_H2.txt', folder + 'out_' + file_name + '_H2.txt')
shutil.copy('./results/out_' + file_name + '_MIl_MI_MIm_MIs.txt', folder + 'out_' + file_name + '_MIl_MI_MIm_MIs.txt')
shutil.copy('./results/out_' + file_name + '_T.txt', folder + 'out_' + file_name + '_T.txt')
shutil.copy('./results/out_' + file_name + '_con_array.txt', folder + 'out_' + file_name + '_con_array.txt')
shutil.copy('./results/out_' + file_name + '_cond_entr2.txt', folder + 'out_' + file_name + '_cond_entr2.txt')
shutil.copy('./results/out_' + file_name + '_dist.txt', folder + 'out_' + file_name + '_dist.txt')
shutil.copy('./results/out_' + file_name + '_taumin.txt', folder + 'out_' + file_name + '_taumin.txt')
shutil.copy('./results/out_' + file_name + '_threshold.txt', folder + 'out_' + file_name + '_threshold.txt')
os.remove(file_name)
os.remove('conmat.txt')
os.remove('./results/out_' + file_name + '_H1.txt')
os.remove('./results/out_' + file_name + '_H2.txt')
os.remove('./results/out_' + file_name + '_MIl_MI_MIm_MIs.txt')
os.remove('./results/out_' + file_name + '_T.txt')
os.remove('./results/out_' + file_name + '_con_array.txt')
os.remove('./results/out_' + file_name + '_cond_entr2.txt')
os.remove('./results/out_' + file_name + '_dist.txt')
os.remove('./results/out_' + file_name + '_taumin.txt')
os.remove('./results/out_' + file_name + '_threshold.txt')
if Args.ert_crit > 1:
    shutil.copy('./results/out_' + file_name + '_H3.txt', folder)
    shutil.copy('./results/out_' + file_name + '_MI3.txt', folder)
    shutil.copy('./results/out_' + file_name + '_cond_entr3.txt', folder)
    os.remove('./results/out_' + file_name + '_H3.txt')
    os.remove('./results/out_' + file_name + '_MI3.txt')
    os.remove('./results/out_' + file_name + '_cond_entr3.txt')
    if Args.ert_crit > 2:
        shutil.copy('./results/out_' + file_name + '_H4.txt', folder)
        shutil.copy('./results/out_' + file_name + '_cond_entr4.txt', folder)
        os.remove('./results/out_' + file_name + '_H4.txt')
        os.remove('./results/out_' + file_name + '_cond_entr4.txt')

# generate run info file
df_info = pd.DataFrame()
df_info.loc['taumax', 'value'] = Args.taumax
df_info.loc['ert_crit', 'value'] = Args.ert_crit
df_info.to_csv(folder + file_name + '_runinfo.csv')


#################################################################### read in fortran generated files

# get results from fortran files
Output = ResultsClass(folder)

# write prediction file


# plot results
if Args.plotMI < 2:
    plot_mi(Output, show=Args.plotMI)
    plot_graph(Output, show=Args.plotMI)
