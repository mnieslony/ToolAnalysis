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
    # train the model
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
    
    # The easiest way to save a model during training and restore it later for prediction
    # is via checkpoints. Simply pass a desired output directory for the checkpoint files
    # and train calls will write checkpoints, while subsequent predict calls will retrieve them.
    # In fact subsequent train calls will also attempt to use checkpoints to pre-load a starting
    # point, so that training can be split into steps.
    # XXX     We assume each training should be an independent one-shot process     XXX
    # XXX          and wipe the checkpoint directory before training.               XXX
    # The limitation of this method is that the Estimator defined before 'predict' calls
    # must match that defined before 'train' calls, for the checkpoint to load correctly
    # XXX i.e. if you change the DNNRegressor Estimator here (e.g. # hidden layers) XXX
    # XXX  you must change it in DNNFindTrackLenghInWater_pred.py to match as well  XXX
    # you could of course have multiple checkpoint directories for multiple trained models
    # Retrieve the user's desired output location to store the model checkpoint
    checkpointdir = Store.GetStoreVariable('EnergyReco','TrackLengthCheckpointDir')
    evtnum = Store.GetStoreVariable('ANNIEEvent','EventNumber')
    if evtnum == 0:
        print('Clearing any existing checkpoints in...'+checkpointdir)
        shutil.rmtree(checkpointdir, True)  # try to remove checkpoint dir, ignore errors
 
    # Build a fully connected DNN with 2 hidden layers, with 70 and 20 nodes respectively.
    feature_columns = [
       tf.feature_column.numeric_column('x', shape=np.array(train_x).shape[1:])]
    regressor = tf.estimator.DNNRegressor(
       feature_columns=feature_columns, hidden_units=[70, 20],model_dir=checkpointdir)

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
    if testingdatafilename != 'NA':
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
        # construct the test data structure
        test_input_fn = tf.estimator.inputs.numpy_input_fn(
             x={'x': testfeatures_transformed}, y=testlabels, shuffle=False)
        
        # Score with tensorflow.
        scores = regressor.evaluate(input_fn=test_input_fn)
        print('MSE (tensorflow): {0:f}'.format(scores['average_loss']))
        Store.SetStoreVariable('EnergyReco','DNNTrackLengthTrainScore',scores['average_loss'])

    return 1

