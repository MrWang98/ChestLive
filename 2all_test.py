"""
    解决：如果有100个类，如何划分，训练集、验证集和测试集。可以60个类训练、20个验证，最后20个测试。
"""
# import keras
# from keras_preprocessing import image
#import matplotlib.pyplot as plt
import numpy as np
import random
import tensorflow as tf
#import keras
#from keras import layers
# from keras_preprocessing import image
from tensorflow import keras
from tensorflow.keras import layers
#import tensorflow_datasets as tfds
import Func_ReadAllData
import ReadData
from tensorflow.keras.layers import Input, Conv1D, Activation, Dropout, MaxPooling1D, Conv2D, MaxPooling2D
import copy

import os

os.environ['CUDA_VISIBLE_DEVICES']='1'
gpus= tf.config.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(gpus[0], True)

"""
## Define the Hyperparameters
"""
learning_rate = 0.003 # 0.003
meta_step_size = 0.1 # 0.25

inner_batch_size = 25
eval_batch_size = 25

meta_iters = 2000  # 8000
eval_iters = 5
inner_iters = 4

eval_interval = 1
train_shots = 20  # 20 5
shots = 20  # 5

classes = 2  # 5

"""
## Change the model
"""
# X_input = Input(input_shape)
inputs = layers.Input(shape=(2000, 39, 1)) #, 1
X = Conv2D(filters=32, kernel_size=4, strides=1, name='conv0')(inputs) # 32
X = Activation('relu')(X)
# X = Dropout(0.3)(X)
X = MaxPooling2D(pool_size=2, strides=2, name='max_pool0')(X) #pool_size=2, strides=2 (4)

X = Conv2D(256, 4, strides=1, name='conv1')(X)  # filter setting. 64
X = Activation('relu')(X)
# X = Dropout(0.3)(X)
X = MaxPooling2D(pool_size=2, strides=2, name='max_pool01')(X)

X = Conv2D(128, 4, strides=1, name='conv11')(X)  # filter setting. #128
X = Activation('relu')(X)
# X = Dropout(0.3)(X)
X = MaxPooling2D(pool_size=2, strides=2, name='max_pool1')(X)

X = Conv2D(64, 4, strides=1, name='conv2', padding="same")(X)  # filter setting. # 256
X = Activation('relu')(X)
# X = Dropout(0.3)(X)

X = MaxPooling2D(pool_size=2, strides=2, name='max_pool')(X) #2
X = layers.Flatten()(X)
X = layers.Dense(classes, activation='softmax', name='fc')(X)
model = keras.Model(inputs=inputs, outputs=X, name='VoiceChestModel_Keras')

model.compile(metrics=['accuracy'])
optimizer = keras.optimizers.SGD(learning_rate=learning_rate)

import copy

class testDataset:
    # This class will facilitate the creation of a few-shot dataset
    # from Collected data that can be sampled from quickly while also
    # allowing to create new labels at the same time.
    def __init__(self, training):
        self.data = {}
        self.test_data = {}
        data_train, labels_train = ReadData.loadData(startPath="split2",training=training)
        data_test, labels_test = ReadData.loadData(startPath=os.path.join("Set"),training=training)
        for data_tmp, label_tmp in zip(data_train, labels_train):
            data_tmp = np.expand_dims(data_tmp, axis=2)
            if label_tmp not in self.data:
                self.data[label_tmp] = []
            self.data[label_tmp].append(data_tmp)
            self.labels = list(self.data.keys())
        
        for data_tmp, label_tmp in zip(data_test, labels_test):
            data_tmp = np.expand_dims(data_tmp, axis=2)
            if label_tmp not in self.test_data:
                self.test_data[label_tmp] = []
            self.test_data[label_tmp].append(data_tmp)
            self.test_labels = list(self.test_data.keys())

            

    def get_mini_dataset(
        self, batch_size, repetitions, shots, num_classes, split=False, test=False, test_size=20
    ):
        testImage = {}
        testLabels = {}
        temp_labels = np.zeros(shape=(num_classes * shots))
        temp_images = np.zeros(shape=(num_classes * shots, 2000, 39, 1))
        if split:
            keys = [0 for i in range(num_classes)]
            test_labels = np.zeros(shape=(num_classes))
            test_images = np.zeros(shape=(num_classes, 2000, 39, 1))

        # Get a random subset of labels from the entire label set.
        # label_subset = random.choices(self.labels, k=num_classes)  # k 次有放回采样
        label_subset =[]
        label_subset = np.array(label_subset)
        tmpLabels = copy.deepcopy(self.labels)
        for i in range(num_classes):                                # num_classes 次不放回采样
            tmp = random.choices(tmpLabels,k=1)
            label_subset = np.concatenate((label_subset,tmp),axis=0)
            tmpLabels.remove(tmp[0])
        for class_idx, class_obj in enumerate(label_subset):
            if class_obj not in testImage:
              # testImage[class_obj] = []
              # testLabels[class_obj] = []
              testImage[class_obj] = np.zeros(shape=(test_size, 2000, 39, 1))
              testLabels[class_obj] = np.zeros(shape=(test_size))
            # Use enumerated index value as a temporary label for mini-batch in
            # few shot learning.
            temp_labels[class_idx * shots : (class_idx + 1) * shots] = class_idx

            # If creating a split dataset for testing, select an extra sample from each label to create the test dataset.
            if split:
                # test_labels[class_idx] = class_idx
                test_labels[class_idx] = class_idx
                # images_to_split = random.choices(
                #     self.data[label_subset[class_idx]], k=shots + 1   # 选出shot + 1个data
                # )
                images_to_split = self.data[label_subset[class_idx]][0:shots+1]
                
                if test:
                  keys[class_idx] = class_obj
                  print("keys:{}".format(keys))

                  if(label_subset[class_idx]=="WZY"):
                    print(label_subset[class_idx]+"_hungry")
                    testImage[class_obj] = random.choices(self.test_data[label_subset[class_idx]+"_hungry"],k=test_size)
                    testLabels[class_obj][0:] = class_idx
                  else:
                    testImage[class_obj] = random.choices(self.data[label_subset[class_idx]][shots+1:],k=test_size)
                    testLabels[class_obj][0:] = class_idx

                  # testImage[class_obj] = random.choices(self.data[label_subset[class_idx]][shots+1:],k=test_size)
                  # testLabels[class_obj][0:] = class_idx

                test_images[class_idx] = images_to_split[-1]
                temp_images[
                    class_idx * shots : (class_idx + 1) * shots
                ] = images_to_split[:-1]
            else:
                # For each index in the randomly selected label_subset, sample the
                # necessary number of images.
                temp_images[
                    class_idx * shots : (class_idx + 1) * shots
                ] = random.choices(self.data[label_subset[class_idx]], k=shots)

        # print(keys)
        dataset = tf.data.Dataset.from_tensor_slices(
            (temp_images.astype(np.float32), temp_labels.astype(np.int32))
        )
        dataset = dataset.shuffle(100).batch(batch_size).repeat(repetitions)  #打乱，喂入batch size，迭代repetitions次
        if test:
          tmpImage = []
          tmpLabels = []
          # print(len(list(testLabels.values())[0]))
          l = len(list(testLabels.values())[0])
          # print(l)
          for i in range(l):
            for k in label_subset:
              # print("i:{},lemgth:{}".format(i,len(list(testImage[k]))))
              # print("length:{},idx:{}".format(len(list(testImage[k])),i))
              tmpImage.append(list(testImage[k])[i])
              tmpLabels.append(list(testLabels[k])[i])
          # print(len(tmpImage))
          tmpImage = np.array(tmpImage)
          tmpLabels = np.array(tmpLabels)
          testset = tf.data.Dataset.from_tensor_slices((tmpImage.astype(np.float32), tmpLabels.astype(np.float32)))
          testset = testset.batch(num_classes)
          keys = np.array(keys)
          return dataset, test_images, test_labels, keys, testset
        if split:
            return dataset, test_images, test_labels
        return dataset

testSet = testDataset(False)

model = tf.keras.models.load_model("checkpoint/model_2way.h5")

number = 1000
test_size = 5
acc = {}
# for shots in range(1,21,1):
#   print("shots:{}".format(shots))
shots=20
# test_size = 45


all_count = {}  #保存各个类被抽取出的总次数
label_dic = {}  #保存各个类被抽取出来后预测正确的次数，与上面的总次数相除暂定为预测的准确度
true_list=[]
score_list=[]
person = {}
for i in range(number):
  testing = []
  train_set, test_images, test_labels, keys, test_set = testSet.get_mini_dataset(
        eval_batch_size, eval_iters, shots, classes, split=True, test=True, test_size=test_size
    )
  # print(test_set)
  for l in keys:
    if l not in label_dic:
      label_dic[l]=0
    if l not in all_count:
      all_count[l]=0
    if l not in person:
        person[l]=[]
  old_vars1 = model.get_weights()
  # train
  for images, labels in train_set:
    # print(labels)
    with tf.GradientTape() as tape:
        preds = model(images)
        loss = keras.losses.sparse_categorical_crossentropy(labels, preds)
    # print(loss)
    # print("labels:{}, preds:{}".format(labels,preds))
    grads = tape.gradient(loss, model.trainable_weights)
    optimizer.apply_gradients(zip(grads, model.trainable_weights))
  test_preds = model.predict(test_images)
  test_preds = tf.argmax(test_preds).numpy()
  # if(i%50==0):
  #   print(i)
  print(i)
  # print("labels:{}, preds:{}".format(test_labels,test_preds))
  # test
  for images, labels in test_set:
      old_vars2 = model.get_weights()
      for label in labels:
        # print(label.numpy())
        all_count[keys[label.numpy().astype(int)]]+=1 
      preds = model.predict(images)
      


      person[keys[0]].append([labels[0],preds[0][1]])
      person[keys[1]].append([labels[1],preds[1][1]])
      
      
      true_list.append(labels[0])
      score_list.append(preds[0])
      true_list.append(labels[1])
      score_list.append(preds[1])
      preds = tf.argmax(preds).numpy()
      print("labels:{}, preds:{}".format(labels,preds))
      flags = (labels==preds)
      num_correct = flags.numpy().sum()
      testing.append(num_correct/classes)
      for idx, flag in enumerate(flags):
        if flag:
          label_dic[keys[idx]]+=1
      model.set_weights(old_vars2)
  model.set_weights(old_vars1)
filelist=os.listdir("./split/test/")
person_number=len(filelist)
result = {}
for key in label_dic.keys():
  result[key] = (label_dic[key]/all_count[key])*100
  print("label:{},  抽取出来的总次数:{},  其中预测正确的次数:{},  acc:{:.2f}%".format(key,all_count[key],label_dic[key],(label_dic[key]/all_count[key])*100))
with open("npy/r2.txt","w") as f:
    f.write(str(result))
with open("npy/result_{}.txt".format(person_number),"w") as f:
    f.write(str(person))
#np.save("npy/true.npy",true_list)
#np.save("npy/score.npy",score_list)
