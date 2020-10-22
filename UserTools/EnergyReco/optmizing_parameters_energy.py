##### Script for optimising parameters for Energy BDT
import Store
import sys
import numpy as np
import pandas as pd
import tempfile
import random
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from array import array
from sklearn import datasets
from sklearn import metrics
from sklearn import model_selection
from sklearn import preprocessing
import sklearn
from sklearn.utils import shuffle
from sklearn import linear_model, ensemble
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
import pickle

def Initialise():
    return 1

def Finalise():
    return 1

def optimise_model(model, alg_name, search_method, params, features_train, labels_train):

    #---- GridSearchCV ----
    random_search = GridSearchCV(model, param_grid=params, scoring='neg_mean_squared_error', n_jobs=4, cv=5, verbose=3)
    random_search.fit(features_train, labels_train)

    print("------ Algorithm: " + str(alg_name))
    #print('\n All results:')
    #print(random_search.cv_results_)
    print('\n Best estimator:')
    print(random_search.best_estimator_)
    print('\n Best hyperparameters:')
    print(random_search.best_params_)
    #results = pd.DataFrame(random_search.cv_results_)
    #results.to_csv('xgb-random-grid-search-results-01.csv', index=False)
    print("----------------------------------------------")

def Execute(Toolchain=True, optimizationdatafilename=None, E_threshold=None):
    # Set TF random seed to improve reproducibility
    seed = 170
    np.random.seed(seed)

    #--- load events for training
    if Toolchain:
        optimizationdatafilename = Store.GetStoreVariable('Config','MuonEnergyOptimizationDataFile')
    print( "--- opening training file "+optimization)
    optimizationfile = open(optimization)
    print("evts for optimization in: ",optimizationfile)
    #--- load csv file into pandas DataFrame
    optimizationfiledata=pd.read_csv(optimizationfile)
    optimizationfile.close()
    #--- Extract named columns
    OptimizationDataset=optimizationfiledata[['totalPMTs','totalLAPPDs','TrueTrackLengthInWater','neutrinoE','trueKE','diffDirAbs','TrueTrackLengthInMrd','recoDWallR','recoDWallZ','dirX','dirY','dirZ','vtxX','vtxY','vtxZ','DNNRecoLength']]
    #--- place an upper limit on Muon energies we will train on, filtering only passing rows
    if Toolchain:
        E_threshold = Store.GetStoreVariable('Config','BDT_NuE_threshold')
    #dfsel_train=OptimizationDataset.loc[OptimizationDataset['neutrinoE'] < E_threshold]
    dfsel_train=OptimizationDataset

    #--- print to check:
    print("check training sample: ",dfsel_train.head())
#    print(dfsel_train.iloc[5:10,0:5])
    # check fr NaN values:
    assert(dfsel_train.isnull().any().any()==False)

    # apply normalization scalings
    dfsel_train_normalised = pd.DataFrame([ dfsel_train['DNNRecoLength']/600., dfsel_train['TrueTrackLengthInMrd']/200., dfsel_train['diffDirAbs'], dfsel_train['recoDWallR'], dfsel_train['recoDWallZ'], dfsel_train['totalLAPPDs']/200., dfsel_train['totalPMTs']/200., dfsel_train['vtxX']/150., dfsel_train['vtxY']/200., dfsel_train['vtxZ']/150. ]).T

    print("check optimization sample normalisation: ", dfsel_train_normalised.head())

    #--- convert features and labels DataFrames to numpy arrays:
    features_train = np.array(dfsel_train_normalised[['DNNRecoLength','TrueTrackLengthInMrd','diffDirAbs','recoDWallR','recoDWallZ','totalLAPPDs','totalPMTs','vtxX','vtxY','vtxZ']])
    labels_train = np.array(dfsel_train[['trueKE']])
    print('events for training: ',len(labels_train))
    print("features_train.shape: ",features_train.shape)
    
    ########### BDTG ############
    params = {'n_estimators':[100, 200, 500, 600, 1000], 
              'max_depth': [10, 50, 100],
              'learning_rate': [0.01, 0.025, 0.05, 0.1], 
              'loss': ['lad']}
    
    # optimize
    print("optimizing BDTG...")
    net_hi_E = ensemble.GradientBoostingRegressor(**params)
    optimise_model(net_hi_E, " BDTG ", "grid", params, features_train, labels_train)
    #---------------------------#

    return 1

if __name__ == "__main__":
    # Make the script runnable as a standalone python script too?
    optimizationdatafilename =  '../TrackLengthReconstruction/vars_Ereco_train_05202019.csv'
    E_threshold=2.
    Execute(False, optimizationdatafilename, E_threshold)
