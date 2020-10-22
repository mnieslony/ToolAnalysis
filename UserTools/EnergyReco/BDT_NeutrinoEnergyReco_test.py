##### Script for Neutrino Energy Reconstruction in the water tank
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
import pickle

def Initialise():
    return 1

def Finalise():
    return 1

def Execute(Toolchain=True, testingdatafilename=None, E_threshold=None, modelfilename=None, predictionsdatafilename=None):

    #--- get testing data file name
    if Toolchain:
        testingdatafilename = Store.GetStoreVariable('Config','NeutrinoEnergyTestingDataFile')
    
    #--- load events for testing
    print( "--- Opening test sample file" + testingdatafilename)
    testingfile = open(testingdatafilename)
    testingfiledata=pd.read_csv(testingfile)
    testingfile.close()
    TestingDataset=testingfiledata[['totalPMTs', 'totalLAPPDs', 'TrueTrackLengthInWater', 'neutrinoE', 'trueKE', 'diffDirAbs', 'TrueTrackLengthInMrd', 'recoDWallR', 'recoDWallZ', 'dirX', 'dirY', 'dirZ', 'vtxX', 'vtxY', 'vtxZ', 'DNNRecoLength']]
    dfsel_test=TestingDataset.loc[TestingDataset['neutrinoE'] < E_threshold]
    print("check testing sample: ",dfsel_test.head())
    assert(dfsel_test.isnull().any().any()==False)

    #--- normalise the sample parameters:
    dfsel_test_normalised = pd.DataFrame([ dfsel_test['DNNRecoLength']/600., dfsel_test['TrueTrackLengthInMrd'], dfsel_test['diffDirAbs'], dfsel_test['recoDWallR'], dfsel_test['recoDWallZ'], dfsel_test['totalLAPPDs']/1000., dfsel_test['totalPMTs']/1000., dfsel_test['vtxX']/150., dfsel_test['vtxY']/200., dfsel_test['vtxZ']/150. ]).T
    print("check test sample normalisation: ", dfsel_test_normalised.head())

    #--- convert features and labels to numpy arrays:
    features_test = np.array(dfsel_test_normalised[['DNNRecoLength', 'TrueTrackLengthInMrd', 'diffDirAbs', 'recoDWallR', 'recoDWallZ', 'totalLAPPDs', 'totalPMTs', 'vtxX', 'vtxY', 'vtxZ']])
    labels_test = np.array(dfsel_test[['trueKE']])
    print('events for testing: ',len(labels_test))
    print("features_test.shape: ",features_test.shape)

    ########### BDTG ############
    # read model from the disk
    modelfilename = Store.GetStoreVariable('Config','BDTNeutrinoModelFile')
    loaded_model = pickle.load(open(modelfilename, 'rb'))
    
    #############################

    #---- Predict on the test set
    BDTGoutput_E = loaded_model.predict(features_test)
    
    #---- print MSE
    mse = mean_squared_error(labels_test, BDTGoutput_E)
    print("MSE: %.4f" % mse)

    #---- Make histogram of fractional error
    Y=[0 for j in range (0,len(test_data_trueKE_hi_E))]
    for i in range(len(test_data_trueKE_hi_E)):
        Y[i] = 100.*(test_data_trueKE_hi_E[i]-BDTGoutput_E[i])/(1.*test_data_trueKE_hi_E[i])
    nbins=np.arange(-100,100,2)
    fig,ax0=plt.subplots(ncols=1, sharey=True)#, figsize=(8, 6))
    cmap = sns.light_palette('b',as_cmap=True)
    f=ax0.hist(np.array(Y), nbins, histtype='step', fill=True, color='gold',alpha=0.75)
    ax0.set_xlim(-100.,100.)
    ax0.set_xlabel('$\Delta E/E$ [%]')
    ax0.set_ylabel('Number of Entries')
    ax0.xaxis.set_label_coords(0.95, -0.08)
    ax0.yaxis.set_label_coords(-0.1, 0.71)
    title = "mean = %.2f, std = %.2f " % (np.array(Y).mean(), np.array(Y).std())
    plt.title(title)
    plt.savefig("BDT_Neutrino_DE_E.png")
 
#    #---- Score and plot loss
#    test_score = np.zeros((params['n_estimators'],), dtype=np.float64)
#    for i, y_pred in enumerate(loaded_model.staged_predict(features_test)):
#        test_score[i] = loaded_model.loss_(labels_test, y_pred)
#    fig,ax=plt.subplots(ncols=1, sharey=True)
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

    # write predictions to the old-style csv file, for validation while we migrate
    if Toolchain:
        predictionsdatafilename = Store.GetStoreVariable('Config','MuonEnergyPredictionsFile')
    if predictionsdatafilename == 'NA':
        return 1  # if not saving to legacy file, just return
    # build a pandas DataFrame from the predicted result
    df1 = pd.DataFrame(labels,columns=['NeutrinoEnergy'])
    df2 = pd.DataFrame(BDTGoutput_E,columns=['RecoE'])
    df_final = pd.concat([df1,df2],axis=1)
    df_final.to_csv(predictionsdatafilename, float_format = '%.3f')

    return 1

if __name__ == "__main__":
    # Make the script runnable as a standalone python script too?
    testingdatafilename = '../LocalFolder/BDT_testing_input.csv'
    predictionsdatafilename = '../LocalFolder/E_Nu_reco_results.csv'
    modelfilename = '../LocalFolder/finalized_BDTmodel_forNeutrinoEnergy.sav'
    E_threshold=2.
    Execute(False, testingdatafilename, E_threshold, modelfilename, predictionsdatafilename)
