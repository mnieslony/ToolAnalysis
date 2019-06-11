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
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.callbacks import ModelCheckpoint
from sklearn.model_selection import GridSearchCV
from tensorflow.keras.wrappers.scikit_learn import KerasRegressor

def Initialise():
    return 1

def Finalise():
    return 1

def create_model(optimizer='Adamax', init_mode='he_normal', activation='relu', neurons1=50, neurons2=15):
    # create model
    model = Sequential()
    model.add(Dense(neurons1, input_dim=2203, kernel_initializer=init_mode, activation=activation))
    model.add(Dense(neurons2, kernel_initializer=init_mode, activation=activation))
    model.add(Dense(1, kernel_initializer=init_mode, activation=activation))
    # Compile model
    model.compile(loss='mean_squared_error', optimizer=optimizer, metrics=['accuracy'])
    return model

def Execute(Toolchain=True, optimizationdatafilename=None):
    # train the model
    # Set TF random seed to improve reproducibility
    seed = 150
    np.random.seed(seed)

    #--- events for optimization - MC events
    # get optimization data file path from store
    if Toolchain:
        optimizationdatafilename = Store.GetStoreVariable('Config','TrackLengthOptimizationDataFile')
    # open the file
    optimizationfile = open(optimizationdatafilename)
    print("evts for optimization in: ",optimizationfile)
    # read into a pandas structure
    optimizationfiledata = pd.read_csv(optimizationfile)
    optimizationfile.close()
    # convert to 2D numpy array
    OptimizationDataset = np.array(optimizationfiledata)
    # split the numpy array up into sub-arrays
    features, lambdamax, labels, rest = np.split(OptimizationDataset,[2203,2204,2205],axis=1)
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

    # Scale the optimization set to 0 mean and unit standard deviation.
    scaler = preprocessing.StandardScaler()
    train_x = scaler.fit_transform(train_x)
    
    # Construct the DNN model
    estimator = KerasRegressor(build_fn=create_model, epochs=10, batch_size=2, verbose=0)

    #--- define the grid search parameters: comment out one option at a time to reduce CPU time!
    #1) tune batch_size, epochs
    batch_size = [1, 2, 5]
    epochs = [10, 50, 100, 500]
    param_grid = dict(batch_size=batch_size, epochs=epochs)

    #2) tune the optimization algorithm, the weight initialisation, 
    optimizer = ['SGD', 'RMSprop', 'Adagrad', 'Adadelta', 'Adam', 'Adamax', 'Nadam']
    init_mode = ['uniform', 'lecun_uniform', 'normal', 'zero', 'glorot_normal', 'glorot_uniform', 'he_normal', 'he_uniform']
    #param_grid = dict(optimizer=optimizer, init_mode=init_mode)

    #3) tune the activation function, the number of neurons for the 1st layer
    activation = ['softmax', 'softplus', 'softsign', 'relu', 'tanh', 'sigmoid', 'hard_sigmoid', 'linear']
    neurons = [25, 50, 70, 90, 100]
    #param_grid = dict(activation=activation, neurons1=neurons)

    #4) tune the number of neurons for the 2nd layer
    neurons2 = [5, 10, 15, 20, 25, 30]
    #param_grid = dict(neurons2=neurons2)

    # search the grid parameters:
    grid = GridSearchCV(estimator=estimator, param_grid=param_grid, n_jobs=-1)
    grid_result = grid.fit(train_x, train_y)
    # summarize results
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
    means = grid_result.cv_results_['mean_test_score']
    stds = grid_result.cv_results_['std_test_score']
    params = grid_result.cv_results_['params']
    for mean, stdev, param in zip(means, stds, params):
        print("%f (%f) with: %r" % (mean, stdev, param))

    return 1

if __name__ == "__main__":
    # Make the script runnable as a standalone python script too?
    optimizationdatafilename = '../LocalFolder/data_forRecoLength_05202019.csv'
    Execute(False, optimizationdatafilename)
