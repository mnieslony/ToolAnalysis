#!/usr/bin/env python
# coding: utf-8

# In[95]:


##### Tool CSV -> pickle

import numpy as np
import pickle
import itertools
import os


import tensorflow as tf
from sklearn.metrics import confusion_matrix



#Hello,
#Please check for yourself the paths and labels up to line 35. 

Label1="SR" #or Electron Y=[1,0]
Label2="MR" # or Muon Y=[0,1]
#Electron/muon labels
#Label1="Electron"
#Label2="Muon"
MC=True


########### Load MODEL ###########
print("Loading Model")
model = tf.keras.models.load_model("CNN-(3, 3)-filter_size-2-double_conv-100-nodes-1-dense0.77.model")
print("Model successfully loaded")
##################################



########### Load Data to Check/Test ###########
#X=Trainingdata, Y=Labels

print("Loading Data")
X= pickle.load(open("cnn_pickle/X_Test_12.13.2020.pickle","rb"))
if MC == True:
	Y= pickle.load(open("cnn_pickle/Y_Test_12.13.2020.pickle","rb"))
print("Data successfully loaded")
##################################





##########################################################################################


##### Prediction


y_prob = np.array(model.predict(X,batch_size=32, verbose=0))
y_classes = y_prob.argmax(axis=-1)
print("Prediction of Data successfull")

file1 = open("CNN_RC_Beam_Prediction.txt","w")
#file1 = open("CNN_PID_Beam_Prediction.txt","w")
if MC == True:
	rounded_labels =np.argmax(Y, axis=1)
	cm = confusion_matrix(rounded_labels, y_classes)
	print("Creation of Confusion Matrix successfull")
	file1.write("True {} correctly identifed: {} \n".format(Label1, cm[0,0]))
	file1.write("True {} wrongly identifed as {} : {} \n".format(Label1,Label2, cm[0,1]))
	file1.write("True {} correctly identifed: {} \n".format(Label2, cm[1,1]))
	file1.write("True {} wrongly identifed as {} : {} \n".format(Label2,Label1, cm[1,0]))
counter = 0
for row in y_classes:
    if row == 0:
        file1.write("{},{},{} \n".format(counter,Label1,y_prob[counter]))
    else:
        file1.write("{},{},{} \n".format(counter,Label2,y_prob[counter]))
    counter = counter+1
    
file1.close()
print("Saving successfull")
print("Fin")

