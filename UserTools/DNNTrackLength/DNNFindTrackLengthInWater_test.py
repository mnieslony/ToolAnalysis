##### Script To Validate DNN for Track Length Reconstruction in the water tank
import Store
import sys
import glob
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
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.wrappers.scikit_learn import KerasRegressor

def Initialise():
    return 1

def Finalise():
    return 1

def Execute(Toolchain=True, testingdatafilename=None, weightsfilename=None, predctionsdatafilename=None):
    
    # Load Data
    #-----------------------------
    if Toolchain:
        testingdatafilename = Store.GetStoreVariable('Config','TrackLengthTestingDataFile')
    # open the file
    testfile = open(testdatafilename)
    print("evts for testing in: ",testfile)
    # read into a pandas structure
    testfiledata = pd.read_csv(testfile)
    testfile.close()
    # convert to 2D numpy array
    TestingDataset = np.array(testfiledata)
    # split the numpy array up into sub-arrays
    testfeatures, testlambdamax, testlabels, testrest = np.split(TestingDataset,[2203,2204,2205],axis=1)
    # scale the features
    testfeatures_transformed = scaler.transform(testfeatures)

    # print info
    print( "lambdamax ", lambdamax[:2], labels[:2])
    print(features[0])
    num_events, num_pixels = features.shape
    print(num_events, num_pixels)
    
    # Preprocess data and load model
    #-----------------------------
    
    # rename variables for obfuscation
    test_x = features
    test_y = labels
    print("test sample features shape: ", test_x.shape," test sample label shape: ", test_y.shape)

    # Scale data to 0 mean and unit standard deviation.
    scaler = preprocessing.StandardScaler()
    x_transformed = scaler.transform(test_x)
    
    # define keras model, loading weight from weights file
    model = Sequential()
    model.add(Dense(50, input_dim=2203, kernel_initializer='he_normal', activation='relu'))
    model.add(Dense(5, kernel_initializer='he_normal', activation='relu'))
    model.add(Dense(1, kernel_initializer='he_normal', activation='relu'))

    # load weights
    if Toolchain:
        weightsfilename = Store.GetStoreVariable('Config','TrackLengthWeightsFile')
    model.load_weights(weightsfilename)

    # Compile model
    model.compile(loss='mean_squared_error', optimizer='Adamax', metrics=['accuracy'])
    print("Created model and loaded weights from file"+weightsfilename)

    # Score accuracy / Make predictions
    #----------------------------------
    print('predicting...')
    y_predicted = model.predict(x_transformed)
    
    # estimate accuracy on dataset using loaded weights
    scores = model.evaluate(test_x, test_y, verbose=0)
    print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))
    
    # Score with sklearn.
    score_sklearn = metrics.mean_squared_error(y_predicted, test_y)
    print('MSE (sklearn): {0:f}'.format(score_sklearn))

    # Write to the output csv file
    #-----------------------------
    
    # firstly, maybe we don't want to save the predictions at all. See if we've been given at least one file:
    if Toolchain:
        predctionsdatafilename = Store.GetStoreVariable('Config','TrackLengthPredictionsDataFile')
    if (predctionsdatafilename = None) or (predctionsdatafilename = ''):
        # no output files today
        return 1
    
    # build a dataframe from the true and predicted track lengths
    outputdataarray = np.concatenate((test_y, y_predicted),axis=1)
    outputdataframe=pd.DataFrame(outputdataarray, columns=['TrueTrackLengthInWater','DNNRecoLength'])

    # append as additional columns to the input dataframe
    testfiledata.insert(2217, 'TrueTrackLengthInWater', outputdataframe['TrueTrackLengthInWater'].values, allow_duplicates="True")
    testfiledata.insert(2218, 'DNNRecoLength', outputdataframe['DNNRecoLength'].values, allow_duplicates="True")

    # check if we're splitting the output into two files (for training/testing the BDTs)
    firstfilesentries = Store.GetStoreVariable('Config','FirstFileEntries')
    predctionsdatafilename2 = Store.GetStoreVariable('Config','TrackLengthPredictionsDataFile2')

    # write to csv file(s)
    if (firstfilesentries = None) or (firstfilesentries = 0) or (predctionsdatafilename2 = None) or (predctionsdatafilename2 = ''):
        testfiledata.to_csv(predctionsdatafilename, float_format = '%.3f')
    else:
        testfiledata[firstfilesentries:].to_csv(predctionsdatafilename, float_format = '%.3f')
        testfiledata[:firstfilesentries].to_csv(predctionsdatafilename2, float_format = '%.3f')

    return 1

if __name__ == "__main__":
    # Make the script runnable as a standalone python script too?
    testingdatafilename = '../LocalFolder/data_forRecoLength_05202019.csv'
    predctionsdatafilename = '../LocalFolder/vars_Ereco_train_06082019.csv'
    predctionsdatafilename2 = '../LocalFolder/vars_Ereco_test_06082019.csv'
    weightsfilename = 'weights_bets.hdf5'
    Execute(False, testingdatafilename, weightsfilename, predctionsdatafilename, predctionsdatafilename2)
