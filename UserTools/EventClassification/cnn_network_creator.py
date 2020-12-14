#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import tensorflow as tf
import time
from datetime import datetime
import pickle
import numpy as np
import matplotlib.pyplot as plt
import random
from tensorflow import keras
from tensorflow.keras import regularizers, layers
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout,LeakyReLU, Activation, Flatten, Conv2D, MaxPooling2D, BatchNormalization, AveragePooling2D, GlobalAveragePooling2D
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from sklearn.metrics import confusion_matrix
import itertools

########### Load Data ###########
#X=Trainingdata, Y=Labels
#RingCounting dataset
X= pickle.load(open("/pnfs/annie/persistent/users/dmaksimo/ml_tutorial/pickle/X_RC_PMT_110k_Files_mitTopBottom.pickle","rb"))
Y= pickle.load(open("/pnfs/annie/persistent/users/dmaksimo/ml_tutorial/pickle/Y_RC_PMT_110k_Files_mitTopBottom.pickle","rb"))

#Exemplary muon/electron PID files can be found at the following paths
#X= pickle.load(open("/pnfs/annie/persistent/users/dmaksimo/ml_tutorial/pickle/X_PID_PMT_120k_Files_gn_tf.pickle","rb"))
#Y= pickle.load(open("/pnfs/annie/persistent/users/dmaksimo/ml_tutorial/pickle/Y_PID_PMT_120k_Files_gn_tf.pickle","rb"))
##################################

#CNN Network Paramters
dense_layer = 1        # How many Dense layers?
nodes=512              # Nodes per dense layer
conv_layer = 2 #4      # How many double Conv layers? 
layer_size = 100 #400  # Conv filters per layer
filter_size = (3,3)    # Filter sizes

b_s= 100               # batch_size
ep = 2                 # epochs

#RingCounting labels
Label_1 = "Single-ring" #electron
Label_2 = "Multi-ring"  #muon
#Electron/Muon labels
#Label_1 = "Electron"
#Label_2 = "Muon"

NAME ="CNN-{}-filter_size-{}-double_conv-{}-nodes-{}-dense".format(filter_size,conv_layer, layer_size, dense_layer)

#Data will be shuffeld. Safe shuffeld sets?
safe_suffeld = True    #will be safed with date and time
path = 'cnn_pickle/'



################################# 
#Should work from here on automatacally
#################################









unique, counts = np.unique(Y, return_counts=True, axis=0)
print("How much from each label? ",counts)
########### Shuffle data #########
training_data = list(zip(X, Y))
random.shuffle(training_data)
########### Seperate data ########
a=len(training_data)
X1 =[]
Y1 =[]
for x in training_data[:int(a*0.7)]:  
    X1.append(x[0])
    Y1.append(x[1])  
XTraining = np.array(X1)
YTraining = np.array(Y1)
X2 =[]
Y2 =[]
for x in training_data[int(a*0.7):int(a*0.85)]: 
    X2.append(x[0])
    Y2.append(x[1])  
XVal = np.array(X2)
YVal = np.array(Y2)
X3 =[]
Y3 =[]
for x in training_data[int(a*0.85):]:
    X3.append(x[0])
    Y3.append(x[1])
XTest = np.array(X3)
YTest = np.array(Y3)
print(XTraining.shape,XVal.shape,XTest.shape)
del X,X1,X2,X3,Y,Y1,Y2,Y3

if safe_suffeld == True:
    now = datetime.now()
    pickle_out = open(path+"X_Training_{}.pickle".format(now.strftime("%m.%d.%Y")),"wb")
    pickle.dump(XTraining,pickle_out,protocol=4)
    pickle_out.close()
    pickle_out = open(path+"X_Val_{}.pickle".format(now.strftime("%m.%d.%Y")),"wb")
    pickle.dump(XVal,pickle_out,protocol=4)
    pickle_out.close()
    pickle_out = open(path+"X_Test_{}.pickle".format(now.strftime("%m.%d.%Y")),"wb")
    pickle.dump(XTest,pickle_out,protocol=4)
    pickle_out.close()

if safe_suffeld == True:
    pickle_out = open(path+"Y_Training_{}.pickle".format(now.strftime("%m.%d.%Y")),"wb")
    pickle.dump(YTraining,pickle_out,protocol=4)
    pickle_out.close()
    pickle_out = open(path+"Y_Val_{}.pickle".format(now.strftime("%m.%d.%Y")),"wb")
    pickle.dump(YVal,pickle_out,protocol=4)
    pickle_out.close()
    pickle_out = open(path+"Y_Test_{}.pickle".format(now.strftime("%m.%d.%Y")),"wb")
    pickle.dump(YTest,pickle_out,protocol=4)
    pickle_out.close()

########## Network ################

model = Sequential()
model.add(Conv2D(layer_size,filter_size,strides=1, input_shape= XTraining.shape[1:],activation="relu", padding='same'))                                               
model.add(Conv2D(layer_size,filter_size,padding='same',activation="relu"))    
model.add(MaxPooling2D(pool_size=(2,2),padding='same'))
model.add(BatchNormalization())
model.add(Dropout(0.2))
for l in range(conv_layer-1):                   
    model.add(Conv2D(layer_size,filter_size,padding='same',activation="relu"))
    model.add(Conv2D(layer_size,filter_size,padding='same',activation="relu"))
    model.add(MaxPooling2D(pool_size=(2,2),padding='same'))
    model.add(BatchNormalization())
    model.add(Dropout(0.2))            
#model.add(GlobalAveragePooling2D())
model.add(Flatten())
for l in range(dense_layer-1):
    model.add(Dense(nodes,activation="relu" ))
    model.add(BatchNormalization())
    model.add(Dropout(0.2))
#model.add(Dense(32,activation="relu"))
model.add(BatchNormalization())
model.add(Dropout(0.2))
model.add(Dense(2))
model.add(Activation('softmax'))
#adam = tf.keras.optimizers.Adam(learning_rate=0.01, beta_1=0.9, beta_2=0.999, amsgrad=True, epsilon = 0.001)
model.compile(loss="binary_crossentropy",
             optimizer="adam",
              metrics=['accuracy']
             )   
filepath=NAME+"{val_acc:.2f}.model"
checkpoint = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, mode='max')
#monitor = EarlyStopping(monitor='val_loss', min_delta=1e-3, patience=10, verbose=1, mode='auto', restore_best_weights=False)
model.summary()
history=model.fit(
    XTraining,YTraining,
 validation_data=(XVal,YVal)   
,batch_size=b_s,
shuffle=True,
class_weight='balanced',
callbacks=[
            #monitor,
            checkpoint,
            #tensorboard 
],
epochs= ep)
print(history.history.keys())
# summarize history for accuracy
plt.plot(history.history['acc'])
plt.plot(history.history['val_acc'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()
# summarize history for loss
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()

print("################################### \n \n Finished epochs \n \n ###################################")
print("Testing")
#model = tf.keras.models.load_model(NAME+".model")
score = model.evaluate(XTest, YTest, verbose=False) 
model.metrics_names
print('Test score (loss): ', score[0])    #Loss on test
print('Test accuracy: ', score[1])



# In[ ]:




