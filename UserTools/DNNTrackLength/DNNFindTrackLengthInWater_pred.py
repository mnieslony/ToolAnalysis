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
    # since we're predicting on an event-wise basis, we need to accumulate the square errors
    # so that we can return the MSE in finalise
    sum_square_errors = Store.GetStoreVariable('EnergyReco','TrackLengthInWaterSumSquaredErrors')
    num_events_processed = Store.GetStoreVariable('EnergyReco','EventsProcessed')
    mean_sq_err = sum_square_errors / num_events_processed
    print('MSE (sklearn): {0:f}'.format(mean_sq_err))

    

    

    
    return 1

def Execute():
    
    # Load Data
    #-----------------------------
    
    #--- events for predicting - MC events
    #filein = open(str(infile))
    #print("evts for predicting in: ",filein)
    #filedata = pd.read_csv(filein)
    #Dataset = np.array(filedata)
    #features, lambdamax, labels, rest = np.split(Dataset,[2203,2204,2205],axis=1)
    # The above imports the csv file into one continguous array, then splits it up, putting elements
    # 0-2202 into 'features', element 2203 into 'lambdamax', 2204 into 'labels' and 2205+ into 'rest'
    ## csv file columns are:
    ## 0-1099: hit lambda values, 1100-2199: hit times, 2200: lambda_max, 2201: Num PMT hits,
    ## 2202: Num LAPPD hits, 2203: lambda_max (again), 2204: TrackLengthInWater, 2205+: nuE, muE ... etc
    
    # equivalent using BoostStore variables
    print( "--- loading input variables from store!")
    hit_lambdas = Store.GetStoreVariable('EnergyReco','lambda_vec')  # std::vectors in the Store are returned
    hit_times = Store.GetStoreVariable('EnergyReco','digit_ts_vec')  # as python lists
    lambdamax = Store.GetStoreVariable('EnergyReco','lambda_max')
    num_pmt_hits = Store.GetStoreVariable('EnergyReco','num_pmt_hits')
    num_lappd_hits = Store.GetStoreVariable('EnergyReco','num_lappd_hits')
    features = np.concatenate([hit_lambdas,hit_times,lambdamax,num_pmt_hits,num_lappd_hits], axis=None)
    labels = Store.GetStoreVariable('EnergyReco','TrackLengthInWater')
    
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
    
    # load BDT model from file XXX XXX XXX XXX XXX this section of code is keras based...
    # combined weights + model file
    combinedmodelfilename = Store.GetStoreVariable('EnergyReco','CombinedTrackLengthModelFile')
    regressor = load_model(combinedmodelfilename)
    
    # alternatively split model and weights files
    #modelfilename = Store.GetStoreVariable('EnergyReco','TrackLengthModelFile')
    #weightsfilename = Store.GetStoreVariable('EnergyReco','TrackLengthWeightsFile')
    #modelfile = open(modelfilename,"r")
    #model_as_json_string = modelfile.read()
    #modelfile.close()
    #regressor = model_from_json(model_as_json_string)
    #regressor.load_weights(weightsfilename)

    # Scale data to 0 mean and unit standard deviation.
    scaler = preprocessing.StandardScaler()
    x_transformed = scaler.transform(test_x)
    
    # Do prediction
    #-----------------------------
    print('predicting...')
    test_input_fn = tf.estimator.inputs.numpy_input_fn(
         x={'x': x_transformed}, y=test_y, shuffle=False)
    predictions = regressor.predict(input_fn=test_input_fn)
    y_predicted = np.array(list(p['predictions'] for p in predictions))
    y_predicted = y_predicted.reshape(np.array(test_y).shape)
    
    
    # Score accuracy
    #-----------------------------
    # Since we only predict one event at a time here, we can't really score here
    # the 'sklearn' metric was a MSE, which we can do ourselves: update the MSE calculation
    sum_square_errors = Store.GetStoreVariable('EnergyReco','TrackLengthInWaterSumSquaredErrors')
    sum_square_errors += ((y_predicted-test_y)**2.)
    Store.SetStoreVariable('EnergyReco','TrackLengthInWaterSumSquaredErrors',sum_square_errors)
    # the final result will be calculated and printed in 'finalise'
    
    
    # Store outputs
    #-----------------------------
    
    # put the predicted lengths into the BoostStore
    Store.SetStoreVariable('EnergyReco','DNNRecoLength',y_predicted)
    
    # old method to append predicted track length to the input file and write to new csv
    #print("shapes: ", test_y.shape, ", ", y_predicted.shape)
    #print(" saving .csv file with energy variables..")
    data = np.concatenate((test_y, y_predicted),axis=1)
    df=pd.DataFrame(data, columns=['TrueTrackLengthInWater','DNNRecoLength'])
    #df_final = pd.concat([filedata,df], axis=1).drop(['lambda_max.1'], axis=1)
    
    
    # Backward Compatibility
    #-----------------------------
    # append this entry to the old-style csv file, for validation while we migrate
    ## filedata here was effectively: [features, lambdamax, labels, rest]
    ## where 'rest' was [nuE, muE, diffDirAbs2, TrueTrackLengthInMrd2, recoDWallR2, recoDWallZ2,
    ## dirX, dirY, dirZ, vtxX, vtxY, vtxZ]. These variables aren't normally loaded,
    ## since they're not used in this script. Load them here to produce a complete legacy csv file...
    truenue = Store.GetStoreVariable('EnergyReco','trueNeuE')
    truemue = Store.GetStoreVariable('EnergyReco','trueEnergy')
    diffdirabs = Store.GetStoreVariable('EnergyReco','diffDirAbs2')
    truemrdtracklen = Store.GetStoreVariable('EnergyReco','TrueTrackLengthInMrd2')
    recodwallr2 = Store.GetStoreVariable('EnergyReco','recoDWallR2')
    recodwallz2 = Store.GetStoreVariable('EnergyReco','recoDWallZ2')
    dirvec = Store.GetStoreVariable('EnergyReco','dirVec')  # these are 3-element lists
    vtxvec = Store.GetStoreVariable('EnergyReco','vtxVec')  # 
    df_final = pd.concat([features,lambdamax,labels,truenue,truemue,diffdirabs,truemrdtracklen,recodwallr2,recodwallz2,dirvec,vtxvec,df], axis=1).drop(['lambda_max.1'], axis=1)

    #-logical tests:
    print("checking..."," df0.shape[0]: ",df0.shape[0]," len(y_predicted): ", len(y_predicted))
    assert(df0.shape[0]==len(y_predicted))
    assert(df_final.shape[0]==df.shape[0])

    df_final.to_csv("../LocalFolder/vars_Ereco.csv", float_format = '%.3f', mode='a', header=False) ## append

    #---if asserts fails check dimensions with these print outs:
    #print("df: ",df.head())
    #print(df.iloc[:,2200:])
    #print(df0.head())
    #print(df0.shape)
    #print(df0.iloc[:,2200:])
    #print(df_final.shape)
    #print(df_final.iloc[:,2200:])


    return 1

