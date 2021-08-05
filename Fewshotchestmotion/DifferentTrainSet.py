"""
    1. The index of python list is 0.
    2. Load data function is at .py of Func_ReadFromMatlab
    3. Training model is at .py of BuiltNN_Keras
"""

from keras.models import load_model
import numpy as np
from keras import metrics, optimizers
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from keras.utils import to_categorical
from sklearn.utils import shuffle
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

# Self-Define function.
import confusion_matrix
import Chest_Model_Metrics
import DifferentTrainSet_Train
import DifferentTrainSet_Test
import BuiltNN_Keras

# Load Data Set.
# 1. In two method multi-person & single person.
data_cup_train, cupLabels_train = DifferentTrainSet_Train.LoadData()
data_cup_test, cupLabels_test = DifferentTrainSet_Test.LoadData()
# 1. Rename Data
# 2. Label data in one-hot method.
# 3. Shuffle data randomly. Another method is use np.random.shuffle & np.random.seed
#
# # Train Data
train_x = data_cup_train
train_y = np.array(LabelEncoder().fit_transform(cupLabels_train)) # LabelEncoder LabelBinarizer()
train_y = to_categorical(train_y, num_classes=2)# 15z
x_train, y_train = shuffle(train_x, train_y, random_state=0)
# Test Data
test_x = data_cup_test
y_test = np.array(LabelEncoder().fit_transform(cupLabels_test)) # LabelEncoder LabelBinarizer()
y_test = to_categorical(y_test, num_classes=2)# 15z
x_test, y_test = shuffle(test_x, y_test, random_state=0)

print(y_train)
print('xuemeng')
print(y_test)

print("Data shape", x_train.shape, "Data label", y_train.shape)
print(x_train, y_train)
# Create Instantiation of Training Model.
VoiceChestModel = BuiltNN_Keras.VoiceChestModel_Keras(train_x.shape[1:])
print(train_x.shape[1:])
# Compile Model.
sgd = optimizers.SGD(lr=0.0001, decay=1e-7, momentum=0.9, nesterov=True)
VoiceChestModel.compile(optimizer='sgd',
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
VoiceChestModel.fit(x_train, y_train, epochs=100, batch_size=100)

# Prediction
preds = VoiceChestModel.evaluate(x_test, y_test, batch_size=40, verbose=1, sample_weight=None)
print("Loss = " + str(preds[0]))
print("Test Accuracy = " + str(preds[1]))
print("Test Accuracy = " + str(preds[2]))
print("Precision = " + str(preds[3]), "Recall = " + str(preds[4]), "F1-Score = " + str(preds[5]))
print("TP = " + str(preds[6]), "TN = " + str(preds[7]), "FP = " + str(preds[8]), "FN = " + str(preds[9]))
print("Classification over.")
confusion_matrix.plot_confuse(VoiceChestModel, x_test, y_test)
# VoiceChestModel.save('my_model.h5')

