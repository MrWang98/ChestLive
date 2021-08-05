"""
    1. The index of python list is 0.
    2. Load data function is at .py of Func_ReadFromMatlab
    3. Training model is at .py of BuiltNN_Keras
"""

import numpy as np
from keras import metrics, optimizers
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from keras.utils import to_categorical
from sklearn.utils import shuffle
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import cross_val_score

# Self-Define function.
import confusion_matrix
import Chest_Model_Metrics
import Func_ReadAllData
import Func_ReadFromMatlab
import Func_ReadAllTrainData
import Func_ReadPartData
# import Func_ReadFromSinglePerson
import Func_ReadFromMimic
import Func_ReadSiglePersonData
import BuiltNN_Keras
import BuiltNN_Keras_bak


# Load Data Set.
# 1. In two method multi-person & single person.
# data_cup, cupLabels = Func_ReadFromMatlab.LoadData()
# data_cup, cupLabels = Func_ReadFromMimic.LoadData()
data_cup, cupLabels = Func_ReadAllTrainData.LoadData()
# data_cup_2, cupLabels_2 = Func_ReadAllTestData.LoadData()
# data_cup, cupLabels = Func_ReadPartData.LoadData()
# data_cup, cupLabels = Func_ReadSiglePersonData.LoadData()
# Split Training Data. to_categorical is used for binary classification.
# 1. Rename Data
# 2. Label data in one-hot method.
# 3. Shuffle data randomly. Another method is use np.random.shuffle & np.random.seed
#
# # Train Data
# train_x = data_cup_1
# train_y = np.array(LabelEncoder().fit_transform(cupLabels_1)) # LabelEncoder LabelBinarizer()
# train_y = to_categorical(train_y, num_classes=11)# 15z
# train_x, train_y = shuffle(train_x, train_y, random_state=0)
# # Test Data
# test_x = data_cup_2
# test_y = np.array(LabelEncoder().fit_transform(cupLabels_2)) # LabelEncoder LabelBinarizer()
# test_y = to_categorical(test_y, num_classes=11)# 15z
# test_x, test_y = shuffle(test_x, test_y, random_state=0)

# print(train_y)
# print('xuemeng')
# print(test_y)

# For all
train_x = data_cup
train_y = np.array(LabelEncoder().fit_transform(cupLabels)) # LabelEncoder LabelBinarizer()
train_y = to_categorical(train_y, num_classes=15)# 15z
# train_y = OneHotEncoder(cupLabels , num_classes=15)

print("Data shape", train_x.shape, "Data label", train_y.shape)
print(train_x, train_y)
train_x, train_y = shuffle(train_x, train_y, random_state=0)
x_train, x_test, y_train, y_test = train_test_split(train_x, train_y, test_size=0.3, random_state=0)

# Create Instantiation of Training Model.
VoiceChestModel = BuiltNN_Keras_bak.VoiceChestModel_Keras(x_train.shape[1:])
# Compile Model.
sgd = optimizers.SGD(lr=0.0004, decay=1e-7, momentum=0.88, nesterov=True)
VoiceChestModel.compile(optimizer='adam',
                        loss='categorical_crossentropy',
                        metrics=['accuracy',
                                 Chest_Model_Metrics.acc,
                                 Chest_Model_Metrics.precision,
                                 Chest_Model_Metrics.sensitivity,
                                 Chest_Model_Metrics.f1_socre,
                                 Chest_Model_Metrics.Modle_TP, # True Accept Rate(TAR)
                                 Chest_Model_Metrics.Modle_TN, # True Reject Rate(TRR)
                                 Chest_Model_Metrics.Modle_FP, # False Accept Rate(FAR)
                                 Chest_Model_Metrics.Modle_FN, # False Reject Rate(FRR)
                                 ]) #class_mode='binary_crossentropy' categorical_crossentropy
# Training Model.
VoiceChestModel.fit(x_train, y_train, epochs=50, batch_size=30) #50
# scores = cross_val_score(VoiceChestModel, x_test, y_test, cv=5)
# print(scores)
# Prediction
preds = VoiceChestModel.evaluate(x_test, y_test, batch_size=30, verbose=1, sample_weight=None)
print("Loss = " + str(preds[0]))
print("Test Accuracy = " + str(preds[1]))
print("Test Accuracy = " + str(preds[2]))
print("Precision = " + str(preds[3]), "Recall = " + str(preds[4]), "F1-Score = " + str(preds[5]))
print("TP = " + str(preds[6]), "TN = " + str(preds[7]), "FP = " + str(preds[8]), "FN = " + str(preds[9]))
print("Classification over.")
confusion_matrix.plot_confuse(VoiceChestModel, x_test, y_test)

