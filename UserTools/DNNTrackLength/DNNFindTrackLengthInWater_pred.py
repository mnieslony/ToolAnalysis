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
    # retrieve the accumulated square errors and calculate the MSE
    sum_square_errors = Store.GetStoreVariable('EnergyReco','TrackLengthInWaterSumSquaredErrors')
    evtnum = Store.GetStoreVariable('ANNIEEvent','EventNumber')
    mean_sq_err = sum_square_errors / evtnum
    print('MSE (sklearn): {0:f}'.format(mean_sq_err))
    # This is, of course, meaningless unless we have truth data
    
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
    ## 2202: Num LAPPD hits, 2203: lambda_max (again), 2204: TrueTrackLengthInWater, 2205+: nuE, muE ... etc
    
    # equivalent using BoostStore variables
    print( "--- loading input variables from store!")
    hit_lambdas = Store.GetStoreVariable('EnergyReco','lambda_vec')  # std::vectors in the Store are returned
    hit_times = Store.GetStoreVariable('EnergyReco','digit_ts_vec')  # as python lists
    lambdamax = Store.GetStoreVariable('EnergyReco','lambda_max')
    num_pmt_hits = Store.GetStoreVariable('EnergyReco','num_pmt_hits')
    num_lappd_hits = Store.GetStoreVariable('EnergyReco','num_lappd_hits')
    features = np.concatenate([hit_lambdas,hit_times,lambdamax,num_pmt_hits,num_lappd_hits], axis=None)
    labels = Store.GetStoreVariable('EnergyReco','TrueTrackLengthInWater')
    
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

    # Scale data to 0 mean and unit standard deviation.
    scaler = preprocessing.StandardScaler()
    x_transformed = scaler.transform(test_x)
    
    # define BDT model, loading weight from checkpointdir
    # XXX the Estimator definition here must match the Estimator definition in   XXX
    # XXX the train script used to generate the checkpoints storing the weights! XXX
    checkpointdir = Store.GetStoreVariable('EnergyReco','TrackLengthCheckpointDir')
    print('Loading model weights from '+checkpointdir)
    feature_columns = [
       tf.feature_column.numeric_column('x', shape=np.array(test_x).shape[1:])]
    regressor = tf.estimator.DNNRegressor(
       feature_columns=feature_columns, hidden_units=[70, 20],model_dir=checkpointdir)
    
    # Do prediction
    #-----------------------------
    print('predicting...')
    test_input_fn = tf.estimator.inputs.numpy_input_fn(
         x={'x': x_transformed}, y=test_y, shuffle=False)
    predictions = regressor.predict(input_fn=test_input_fn)
    y_predicted = np.array(list(p['predictions'] for p in predictions))
    y_predicted = y_predicted.reshape(np.array(test_y).shape)
    
    
    # Score accuracy (if the truth label is meaningful)
    #-----------------------------
    # Update the accumulated square error for the MSE calculation (DNN accuracy score)
    evtnum = Store.GetStoreVariable('ANNIEEvent','EventNumber')
    if evtnum==0:
        sum_square_errors = 0
    else:
        sum_square_errors = Store.GetStoreVariable('EnergyReco','TrackLengthInWaterSumSquaredErrors')
    sum_square_errors += ((y_predicted-test_y)**2.)
    # the final result will be calculated and printed in 'finalise'
    
    
    # Store outputs
    #-----------------------------
    # put the predicted lengths and errors if we have them into the BoostStore
    Store.SetStoreVariable('EnergyReco','DNNRecoLength',y_predicted)
    Store.SetStoreVariable('EnergyReco','TrackLengthInWaterSumSquaredErrors',sum_square_errors)
    #print("shapes: ", test_y.shape, ", ", y_predicted.shape)
    
    #-----------------------------
    # Backward Compatibility
    #-----------------------------
    # append this entry to the old-style csv file, for validation while we migrate
    # Append the predicted track length to the input file and write to new csv file
    #print(" saving .csv file with energy variables..")
    outputfilepath = Store.GetStoreVariable('EnergyReco','TrackLengthPredictionsFile')
    if outputfilepath == 'NA':
        return 1  # if not saving to legacy file, just return
    
    # check if this is the first execute iteration. if so, we'll create a header row first.
    if evtnum==0:
        headers=['']
        for i in range(len(hit_lambdas)):
            headers.append('l_'+str(i))
        for i in range(len(hit_lambdas)):
            headers.append('T_'+str(i))
        headers = headers+['lambda_max', 'totalPMTs', 'totalLAPPDs', 'TrueTrackLengthInWater', 'neutrinoE', 'trueKE', 'diffDirAbs', 'TrueTrackLengthInMrd', 'recoDWallR', 'recoDWallZ', 'dirX', 'dirY', 'dirZ', 'vtxX', 'vtxY', 'vtxZ']
        # convert to pandas dataframe
        headersframe = pd.DataFrame(headers)
        # write the headers to file
        headersframe.T.to_csv(outputfilepath,header=False,index=False)

    # combine true and predicted track lengths into a numpy array
    #data = np.concatenate((test_y, y_predicted),axis=1)
    # convert to pandas dataframe with headers
    #df=pd.DataFrame(data, columns=['TrueTrackLengthInWater','DNNRecoLength'])
    # combine with loaded file data ...
    #df_final = pd.concat([filedata,df], axis=1).drop(['lambda_max.1'], axis=1)  # '.1' added by pandas
    # ...
    # but we no longer have the 'filedata' array. This array previously contained:
    # [features, lambdamax, labels, rest]
    # where 'rest' included
    # [nuE, muE, diffDirAbs2, TrueTrackLengthInMrd2, recoDWallR2, recoDWallZ2, dirX, dirY, dirZ, vtxX, vtxY, vtxZ]
    # but many of these variables are no longer loaded as they're not needed in this script.
    # For legacy/testing, we can load them here to produce a complete csv file...
    truenue = Store.GetStoreVariable('EnergyReco','trueNeuE')
    truemue = Store.GetStoreVariable('EnergyReco','trueEnergy')
    diffdirabs = Store.GetStoreVariable('EnergyReco','diffDirAbs2')
    truemrdtracklen = Store.GetStoreVariable('EnergyReco','TrueTrackLengthInMrd2')
    recodwallr2 = Store.GetStoreVariable('EnergyReco','recoDWallR2')
    recodwallz2 = Store.GetStoreVariable('EnergyReco','recoDWallZ2')
    dirvec = Store.GetStoreVariable('EnergyReco','dirVec')  # these are 3-element lists
    vtxvec = Store.GetStoreVariable('EnergyReco','vtxVec')  #
    
    # combine everything into a numpy array
    data = np.concatenate((evtnum,features,labels,truenue,truemue,diffdirabs,truemrdtracklen,recodwallr2,recodwallz2,dirvec,vtxvec,test_y, y_predicted),axis=1)
    # convert to pandas dataframe
    filedata=pd.DataFrame(data)
    
    # then append this event's data row to the output file
    filedata.to_csv(outputfilepath, float_format = '%.3f', mode='a', header=False, index=False)

    # some quick validation tests:
    print("checking..."," df0.shape[0]: ",df0.shape[0]," len(y_predicted): ", len(y_predicted))
    assert(df0.shape[0]==len(y_predicted))
    assert(df_final.shape[0]==df.shape[0])

    #---if asserts fails check dimensions with these print outs:
    #print("df: ",df.head())
    #print(df.iloc[:,2200:])
    #print(df0.head())
    #print(df0.shape)
    #print(df0.iloc[:,2200:])
    #print(df_final.shape)
    #print(df_final.iloc[:,2200:])


    return 1

