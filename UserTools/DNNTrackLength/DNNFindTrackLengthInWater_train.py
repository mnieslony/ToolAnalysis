##### Script to Train DNN for Track Length Reconstruction in the water tank 
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

def create_model():
    # create model
    model = Sequential()
    model.add(Dense(50, input_dim=2203, kernel_initializer='he_normal', activation='relu'))
    model.add(Dense(5, kernel_initializer='he_normal', activation='relu'))
    model.add(Dense(1, kernel_initializer='he_normal', activation='relu'))
    # Compile model
    model.compile(loss='mean_squared_error', optimizer='Adamax', metrics=['accuracy'])
    return model

def Execute(Toolchain=True, trainingdatafilename=None, weightsfilename=None):
    # train the model
    # Set TF random seed to improve reproducibility
    seed = 150
    np.random.seed(seed)

    #--- events for training - MC events
    # get training data file path from store
    if Toolchain:
        trainingdatafilename = Store.GetStoreVariable('Config','TrackLengthTrainingDataFile')
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
    # 2202: Num LAPPD hits, 2203: lambda_max (again), 2204: TrueTrackLengthInWater, 2205+: nuE, muE ... etc
    
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
    
    # Construct the DNN model
    estimator = KerasRegressor(build_fn=create_model, epochs=10, batch_size=2, verbose=0)

    # load weights
    if Toolchain:
        weightsfilename = Store.GetStoreVariable('Config','TrackLengthWeightsFile')
    checkpoint = ModelCheckpoint(weightsfilename, monitor='val_loss', verbose=1, save_best_only=True, save_weights_only=True, mode='auto')
    callbacks_list = [checkpoint]
    # Run the model
    print('training....')
    history = estimator.fit(train_x, train_y, validation_split=0.33, epochs=10, batch_size=2, callbacks=callbacks_list, verbose=0)

    # summarize history for loss
    f, ax2 = plt.subplots(1,1)
    ax2.plot(history.history['loss'])
    ax2.plot(history.history['val_loss'])
    ax2.set_title('model loss')
    ax2.set_ylabel('Performance')
    ax2.set_xlabel('Epochs')
    #ax2.set_xlim(0.,10.)
    ax2.legend(['loss', 'val_loss'], loc='upper left')
    plt.savefig("../LocalFolder/keras_DNN_training_loss.pdf")

    return 1

if __name__ == "__main__":
    # Make the script runnable as a standalone python script too?
    trainingdatafilename = '../LocalFolder/DNN_training_input.csv'
    weightsfilename = '../LocalFolder/weights_bets.hdf5'
    Execute(False, trainingdatafilename, weightsfilename)
