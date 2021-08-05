import tensorflow as tf
from tensorflow import keras
import os
import Func_ReadAllData
import numpy as np

import os
import keras.backend.tensorflow_backend as KTF
os.environ["CUDA_VISIBLE_DEVICES"]="0"
config = tf.compat.v1.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 0.4
session = tf.compat.v1.Session(config=config)
KTF.set_session(session)

class TestDataset:
    # This class will facilitate the creation of a few-shot dataset
    # from Collected data that can be sampled from quickly while also
    # allowing to create new labels at the same time.
    def __init__(self, training):
        self.data = {}
        self.data_train, self.labels_train = Func_ReadAllData.LoadData(training)
        for data_tmp, label_tmp in zip(self.data_train, self.labels_train):
            data_tmp = np.expand_dims(data_tmp, axis=2)
            if label_tmp not in self.data:
                self.data[label_tmp] = []
            self.data[label_tmp].append(data_tmp)
            self.labels = list(self.data.keys())
        self.length = self.data_train.shape[0]//len(self.labels)

    def get_dataset(
        self, batch_size, repetitions, idx
    ):
        dim1,dim2,dim3 = self.data_train.shape
        temp_labels = np.zeros(shape=(self.length,dim2,dim3,1))
        label_subset = self.labels[idx]
        temp_labels[0: self.length] = idx
        temp_images = self.data[label_subset]
        temp_images = np.array(temp_images)

        dataset = tf.data.Dataset.from_tensor_slices(
            (temp_images.astype(np.float32), temp_labels.astype(np.int32))
        )
        dataset = dataset.shuffle(100).batch(batch_size).repeat(repetitions)

        return dataset

testSet = TestDataset(False)

# currentDir = "."
# for file in os.listdir(currentDir):
#     if(file.split(".")[1]=="h5")
modelName = "model_1900.h5"
model = keras.models.load_model(modelName)

# finalDir = "../Final"
# fileList = os.listdir(finalDir)
fileList = list(testSet.data.keys())

test_batch_size = 25
test_iters = 5
print("开始测试：")
for i in range(len(fileList)):
    accAll = 0
    print("{}…………".format(fileList[i]))
    testSample = testSet.get_dataset(test_batch_size, test_iters, i)
    for images, labels in testSample:
        preds = model.predict(images)
        preds = np.array(preds)
        labels = np.array(labels)
        accAll+=(preds==labels).sum()
    print("The acc of {} : {:5.2f}%".format(fileList[i],((accAll/5)/testSet.length)*100))


