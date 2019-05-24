##### Script for Muon Energy Reconstruction in the water tank
import Store
import sys
import numpy as np
import pandas as pd
import tensorflow as tf
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
import pickle

def Initialise():
    return 1

def Finalise():
    return 1

def Execute(Toolchain=True, trainingdatafilename=None, E_threshold=None, modelfilename=None, testingdatafilename=None):
    # Set TF random seed to improve reproducibility
    seed = 180
    np.random.seed(seed)

    #--- load events for training
    if Toolchain:
        trainingdatafilename = Store.GetStoreVariable('EnergyReco','NeutrinoEnergyTrainingDataFile')
    print( "--- opening training file "+trainingdatafilename)
    trainingfile = open(trainingdatafilename)
    print("evts for training in: ",trainingfile)
    #--- load csv file into pandas DataFrame
    trainingfiledata=pd.read_csv(trainingfile)
    trainingfile.close()
    #--- Extract named columns
    TrainingDataset=trainingfiledata[['totalPMTs','totalLAPPDs','TrueTrackLengthInWater','neutrinoE','trueKE','diffDirAbs','TrueTrackLengthInMrd','recoDWallR','recoDWallZ','dirX','dirY','dirZ','vtxX','vtxY','vtxZ','DNNRecoLength']]
    #--- place an upper limit on Muon energies we will train on, filtering only passing rows
    if Toolchain:
        E_threshold = Store.GetStoreVariable('EnergyReco','BDT_NuE_threshold')
    dfsel_train=TrainingDataset.loc[TrainingDataset['neutrinoE'] < E_threshold]

    #--- print to check:
    print("check training sample: ",dfsel_train.head())
#    print(dfsel_train.iloc[5:10,0:5])
    # check fr NaN values:
    assert(dfsel_train.isnull().any().any()==False)

    # apply normalization scalings
    dfsel_train_normalised = pd.DataFrame([ dfsel_train['DNNRecoLength']/600., dfsel_train['TrueTrackLengthInMrd'], dfsel_train['diffDirAbs'], dfsel_train['recoDWallR'], dfsel_train['recoDWallZ'], dfsel_train['totalLAPPDs']/1000., dfsel_train['totalPMTs']/1000., dfsel_train['vtxX']/150., dfsel_train['vtxY']/200., dfsel_train['vtxZ']/150. ]).T

    print("check train sample normalisation: ", dfsel_train_normalised.head())

    #--- convert features and labels DataFrames to numpy arrays:
    features_train = np.array(dfsel_train_normalised[['DNNRecoLength','TrueTrackLengthInMrd','diffDirAbs','recoDWallR','recoDWallZ','totalLAPPDs','totalPMTs','vtxX','vtxY','vtxZ']])
    labels_train = 1000.*np.array(dfsel_train[['NeutrinoEnergy']])
    print('events for training: ',len(labels_train))
    print("features_train.shape: ",features_train.shape)
    
    ########### BDTG ############
    # build model
    n_estimators=1000
    params = {'n_estimators':n_estimators, 'max_depth': 50,
              'learning_rate': 0.01, 'loss': 'lad'}
    
    # train
    print("training BDTG...")
    net_hi_E = ensemble.GradientBoostingRegressor(**params)
    model = net_hi_E.fit(features_train, labels_train)
    net_hi_E

    # save the model to disk
    if Toolchain:
        modelfilename = Store.GetStoreVariable('EnergyReco','BDTNeutrinoModelFile')
    pickle.dump(model, open(modelfilename, 'wb'))
    
    #############################

    # Measure model metrics
    # load testing dataset: same process as for the training sample
    if Toolchain:
        testingdatafilename = Store.GetStoreVariable('EnergyReco','NeutrinoEnergyTestingDataFile')
    if testingdatafilename == 'NA':
        return 1    # if we have no testing sample, we're done
    
    print( "--- Opening test sample file" + testingdatafilename)
    testingfile = open(testingdatafilename)
    testingfiledata=pd.read_csv(testingfile)
    testingfile.close()
    TestingDataset=testingfiledata[['totalPMTs','totalLAPPDs','TrueTrackLengthInWater','neutrinoE','trueKE','diffDirAbs','TrueTrackLengthInMrd','recoDWallR','recoDWallZ','dirX','dirY','dirZ','vtxX','vtxY','vtxZ','DNNRecoLength']]
    dfsel_test=TestingDataset.loc[TestingDataset['neutrinoE'] < E_threshold]
    print("check testing sample: ",dfsel_test.head())
    assert(dfsel_test.isnull().any().any()==False)

    #--- normalise the sample parameters:
    dfsel_test_normalised = pd.DataFrame([ dfsel_test['DNNRecoLength']/600., dfsel_test['TrueTrackLengthInMrd'], dfsel_test['diffDirAbs'], dfsel_test['recoDWallR'], dfsel_test['recoDWallZ'], dfsel_test['totalLAPPDs']/1000., dfsel_test['totalPMTs']/1000., dfsel_test['vtxX']/150., dfsel_test['vtxY']/200., dfsel_test['vtxZ']/150. ]).T
    print("check test sample normalisation: ", dfsel_test_normalised.head())

    #--- convert features and labels to numpy arrays:
    features_test = np.array(dfsel_test_normalised[['DNNRecoLength','TrueTrackLengthInMrd','diffDirAbs','recoDWallR','recoDWallZ','totalLAPPDs','totalPMTs','vtxX','vtxY','vtxZ']])
    labels_test = 1000.*np.array(dfsel_test[['NeutrinoEnergy']])
    print('events for testing: ',len(labels_test))
    print("features_test.shape: ",features_test.shape)

    #---- Predict on the test set and measure MSE
    mse = mean_squared_error(labels_test, net_hi_E.predict(features_test))
    print("MSE: %.4f" % mse)
    if Toolchain:
        Store.SetStoreVariable('EnergyReco','NeutrinoEnergyMSE',mse)
 
    test_score = np.zeros((params['n_estimators'],), dtype=np.float64)
    for i, y_pred in enumerate(net_hi_E.staged_predict(features_test)):
        test_score[i] = net_hi_E.loss_(labels_test, y_pred)
    if Toolchain:
        Store.SetStoreVariable('EnergyReco','NeutrinoEnergyLosses',test_score)

#    fig,ax=plt.subplots(ncols=1, sharey=True)
#    ax.plot(np.arange(params['n_estimators']) + 1, net_hi_E.train_score_, 'b-',
#             label='Training Set Deviance')
#    ax.plot(np.arange(params['n_estimators']) + 1, test_score, 'r-',
#             label='Test Set Deviance')
#     ax.set_ylim(0.,500.)
#     ax.set_xlim(0.,n_estimators)
#     ax.legend(loc='upper right')
#     ax.set_ylabel('Least Absolute Deviations [MeV]')
#     ax.set_xlabel('Number of Estimators')
#     ax.yaxis.set_label_coords(-0.1, 0.6)
#     ax.xaxis.set_label_coords(0.85, -0.08)
#     plt.savefig("deviation_train_test.png")

    return 1

if __name__ == "__main__":
    # Make the script runnable as a standalone python script too?
    trainingdatafilename =  '../LocalFolder/vars_Ereco.csv'
    testingdatafilename = '../LocalFolder/vars_Ereco.csv'
    modelfilename = '../LocalFolder/finalized_BDTmodel_forNeutrinoEnergy.sav'
    E_threshold=2.
    Execute(False, trainingdatafilename, E_threshold, modelfilename, testingdatafilename)
