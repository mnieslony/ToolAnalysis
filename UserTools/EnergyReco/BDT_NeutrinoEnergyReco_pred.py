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

def Execute():
    # Set TF random seed to improve reproducibility
    seed = 180
    np.random.seed(seed)

    # Retrieve features for predicting
    # --------------------------------
    print( "--- loading input variables from store!")
    num_pmt_hits = Store.GetStoreVariable('EnergyReco','num_pmt_hits')
    num_lappd_hits = Store.GetStoreVariable('EnergyReco','num_lappd_hits')
    neutrinoE = Store.GetStoreVariable('EnergyReco','neutrinoE')
    trueKE = Store.GetStoreVariable('EnergyReco','trueKE')
    diffDirAbs = Store.GetStoreVariable('EnergyReco','diffDirAbs')
    TrueTrackLengthInMrd = Store.GetStoreVariable('EnergyReco','TrueTrackLengthInMrd')
    recoDWallR = Store.GetStoreVariable('EnergyReco','recoDWallR')
    recoDWallZ = Store.GetStoreVariable('EnergyReco','recoDWallZ')
    dirX = Store.GetStoreVariable('EnergyReco','dirX')
    dirY = Store.GetStoreVariable('EnergyReco','dirY')
    dirZ = Store.GetStoreVariable('EnergyReco','dirZ')
    vtxX = Store.GetStoreVariable('EnergyReco','vtxX')
    vtxY = Store.GetStoreVariable('EnergyReco','vtxY')
    vtxZ = Store.GetStoreVariable('EnergyReco','vtxZ')
    DNNRecoLength = Store.GetStoreVariable('EnergyReco','DNNRecoLength')
    
    # XXX I don't think we can apply this here...
    #if trueKE<E_threshold:
    #    return 1
    
    # Normalize them
    # --------------
    DNNRecoLength =/ 600.
    num_lappd_hits =/ 1000
    num_pmt_hits =/ 1000
    vtxX =/ 150.
    vtxY =/ 200.
    vtxZ =/ 150.
    #diffDirAbs                    # no scaling
    #TrueTrackLengthInMrd=/200.    # already scaled in FindTrackLengthInWater
    #recoDWallR=/tank_radius       # already scaled in FindTrackLengthInWater
    #recoDWallZ=/2*tank_halfheight # already scaled in FindTrackLengthInWater

    # build numpy arrays of features and labels
    # -----------------------------------------
    features = np.array([DNNRecoLength,TrueTrackLengthInMrd,diffDirAbs,recoDWallR,recoDWallZ,num_lappd_hits,num_pmt_hits,vtxX,vtxY,vtxZ])
    labels = 1000.*np.array([NeutrinoEnergy])

    ########### BDTG ############
    # read model from the disk
    modelfilename = Store.GetStoreVariable('EnergyReco','BDTNeutrinoModelFile')
    #pickle.dump(model, open(filename, 'wb'))
 
    # load the model from disk
    loaded_model = pickle.load(open(modelfilename, 'rb'))
    
    #############################

    # predicting...
    # -------------
    BDTGoutput_E = loaded_model.predict(features)

    # measure accuracy
    # ----------------
    Y = 100.*(trueKE-BDTGoutput_E[0])/(1.*trueE)
#   print("MC Energy: ", labels[0]," Reco Energy: ",BDTGoutput_E[0]," DE/E[%]: ",Y)

    # Pass to the BoostStore
    # ----------------------
    Store.SetStoreVariable('EnergyReco','NeutrinoEnergyReco',BDTGoutput_E)
    Store.SetStoreVariable('EnergyReco','NeutrinoEnergyAccuracy',Y)

    #-----------------------------
    # Backward Compatibility
    #-----------------------------
    # append this entry to the old-style csv file, for validation while we migrate
    outputfilepath = Store.GetStoreVariable('EnergyReco','NeutrinoEnergyPredictionsFile')
    if outputfilepath == 'NA':
        return 1  # if not saving to legacy file, just return
    
    # check if this is the first execute iteration. if so, we'll create a header row first.
    evtnum = Store.GetStoreVariable('ANNIEEvent','EventNumber')
    if evtnum==0:
        headers=['','NeutrinoEnergy', 'RecoE']
        # convert to pandas dataframe
        headersframe = pd.DataFrame(headers)
        # write the headers to file
        headersframe.T.to_csv(outputfilepath,header=False,index=False)

    # build a pandas DataFrame from the predicted result
    df0 = pd.DataFrame([evtnum],columns=[''])
    df1 = pd.DataFrame(labels,columns=['NeutrinoEnergy'])
    df2 = pd.DataFrame(BDTGoutput_E,columns=['RecoE'])
    df_final = pd.concat([df0,df1,df2],axis=1)
    
    #save results to .csv:
    df_final.to_csv(outputfilepath, float_format = '%.3f', mode='a', header=False, index=False)

#    nbins=np.arange(-100,100,2)
#    fig,ax0=plt.subplots(ncols=1, sharey=True)#, figsize=(8, 6))
#    cmap = sns.light_palette('b',as_cmap=True)
#    f=ax0.hist(np.array(Y), nbins, histtype='step', fill=True, color='gold',alpha=0.75)
#    ax0.set_xlim(-100.,100.)
#    ax0.set_xlabel('$\Delta E/E$ [%]')
#    ax0.set_ylabel('Number of Entries')
#    ax0.xaxis.set_label_coords(0.95, -0.08)
#    ax0.yaxis.set_label_coords(-0.1, 0.71)
#    title = "mean = %.2f, std = %.2f " % (np.array(Y).mean(), np.array(Y).std())
#    plt.title(title)
#    plt.savefig("DE_E.png")
 
    return 1


