##### Script Track Length Reconstruction in the water tank 
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

def Initialise():
    return 1

def Finalise():
    return 1

def Execute():
    # Set TF random seed to improve reproducibility
    seed = 150
    np.random.seed(seed)

    #--- events for training - MC events
    # get training data file path from store
    trainingdatafilename = Store.GetStoreVariable('EnergyReco','TrackLengthTrainingDataFile')
    # open the file
    trainingfile = open(trainingdatafilename)
    print("evts for training in: ",trainingfile)
    # read into a pandas structure
    trainingfiledata = pd.read_csv(trainingfile)
    trainingfile.close()
    # convert to 2D numpy array
    TrainingDataset = np.array(trainingfiledata)
    # split the numpy array up into sub-arrays
    features, lambdamax, labels, rest = np.split(TrainingDataset,[2203,2204,2205],axis=1)
    # This puts splits the arrays column-wise as follows:
    # 0-2202 into 'features', element 2203 into 'lambdamax', 2204 into 'labels' and 2205+ into 'rest'
    # csv file columns are:
    # 0-1099: hit lambda values, 1100-2199: hit times, 2200: lambda_max, 2201: Num PMT hits,
    # 2202: Num LAPPD hits, 2203: lambda_max (again), 2204: TrackLengthInWater, 2205+: nuE, muE ... etc
    
    # print info, initialize seed
    print( "lambdamax ", lambdamax[:2], labels[:2])
    print(features[0])
    num_events, num_pixels = features.shape
    print(num_events, num_pixels)
    np.random.seed(0)
    
    # rename variables for obfuscation
    train_x = features
    train_y = labels
    print("train sample features shape: ", train_x.shape," train sample label shape: ", train_y.shape)

    # Scale the training set to 0 mean and unit standard deviation.
    scaler = preprocessing.StandardScaler()
    train_x = scaler.fit_transform(train_x)
 
    #Build 2 layer fully connected DNN with 10, 10 units respectively.
    feature_columns = [
       tf.feature_column.numeric_column('x', shape=np.array(train_x).shape[1:])]
    regressor = tf.estimator.DNNRegressor(
       feature_columns=feature_columns, hidden_units=[70, 20])

    # Train.
    print('training....')
    batch_size = 1#2
    epochs_no= 2000
    n_batches = int(np.ceil(num_events / batch_size))
    train_input_fn = tf.estimator.inputs.numpy_input_fn(
          x={'x': train_x}, y=train_y, batch_size=batch_size, num_epochs=epochs_no, shuffle=False,num_threads=1)
    regressor.train(input_fn=train_input_fn,steps=1000) #1000)
    
    
    # Score accuracy
    #-----------------------------
    # if we want to score the trained model, we need a test set
    testingdatafilename = Store.GetStoreVariable('EnergyReco','TrackLengthTestingDataFile')
    if testingdatafilename != "":
        # open the file
        testfile = open(testdatafilename)
        print("evts for testing in: ",trainingfile)
        # read into a pandas structure
        testfiledata = pd.read_csv(testfile)
        trainingfile.close()
        # convert to 2D numpy array
        TestingDataset = np.array(testfiledata)
        # split the numpy array up into sub-arrays
        testfeatures, testlambdamax, testlabels, testrest = np.split(TestingDataset,[2203,2204,2205],axis=1)
        # scale the features
        testfeatures_transformed = scaler.transform(testfeatures)
        # construct the test data structure
        test_input_fn = tf.estimator.inputs.numpy_input_fn(
             x={'x': testfeatures_transformed}, y=testlabels, shuffle=False)
        
        # Score with tensorflow.
        scores = regressor.evaluate(input_fn=test_input_fn)
        print('MSE (tensorflow): {0:f}'.format(scores['average_loss']))
    
    # Save model
    #-----------------------------
    # combined weights + model file
    combinedmodelfilename = Store.GetStoreVariable('EnergyReco','CombinedTrackLengthModelFile')
    regressor.save(combinedmodelfilename) #  XXX this section is keras based
    
    # alternatively split model and weights files
    #modelfilename = Store.GetStoreVariable('EnergyReco','TrackLengthModelFile')
    #weightsfilename = Store.GetStoreVariable('EnergyReco','TrackLengthWeightsFile')
    #model_as_json_string = regressor.to_json()
    #modelfile = open(modelfilename, 'w')
    #modelfile.write(model_as_json_string)
    #modelfile.close()
    #regressor.save_weights(weightsfilename, overwrite=True)
    
    # N.B. by saving the model as separate json and weights files, we could potentially load
    # this model using lwtnn and run the prediction step with a c++ tool, rather than calling
    # the predict python script. 
    # BUT: to do this we also need an intermediate step where the model json file and weights file
    # are combined with a variable description json file* to create a single lwtnn json file
    # *created manually? the combination can be performed by a keras2json.py script
    

    return 1

