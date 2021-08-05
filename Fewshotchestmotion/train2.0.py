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
from tensorflow.keras.layers import Input, Conv1D, Activation, Dropout, MaxPooling1D, Conv2D, MaxPooling2D

import os
#import keras.backend.tensorflow_backend as KTF
from scipy.io import savemat
file_name = 'TrainAcc.mat'


os.environ["CUDA_VISIBLE_DEVICES"]="0"
#config = tf.compat.v1.ConfigProto()
#config.gpu_options.per_process_gpu_memory_fraction = 0.4
#session = tf.compat.v1.Session(config=config)
#KTF.set_session(session)

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
train_shots = 70  # 20 5
shots = 20  # 5
classes = 2  # 5

class Dataset:
    # This class will facilitate the creation of a few-shot dataset
    # from Collected data that can be sampled from quickly while also
    # allowing to create new labels at the same time.
    def __init__(self, training):
        self.data = {}
        data_train, labels_train = Func_ReadAllData.LoadData(training=training)
        for data_tmp, label_tmp in zip(data_train, labels_train):
            data_tmp = np.expand_dims(data_tmp, axis=2)
            if label_tmp not in self.data:
                self.data[label_tmp] = []
            self.data[label_tmp].append(data_tmp)
            self.labels = list(self.data.keys())

    def get_mini_dataset(
        self, batch_size, repetitions, shots, num_classes, split=False
    ):
        temp_labels = np.zeros(shape=(num_classes * shots))
        temp_images = np.zeros(shape=(num_classes * shots, 2000, 39, 1))
        if split:
            test_labels = np.zeros(shape=(num_classes))
            test_images = np.zeros(shape=(num_classes, 2000, 39, 1))

        # Get a random subset of labels from the entire label set.
        label_subset = random.choices(self.labels, k=num_classes)  # k 次有放回采样
        for class_idx, class_obj in enumerate(label_subset):
            # Use enumerated index value as a temporary label for mini-batch in
            # few shot learning.
            temp_labels[class_idx * shots : (class_idx + 1) * shots] = class_idx
            # If creating a split dataset for testing, select an extra sample from each
            # label to create the test dataset.
            if split:
                test_labels[class_idx] = class_idx
                images_to_split = random.choices(
                    self.data[label_subset[class_idx]], k=shots + 1
                )
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

        dataset = tf.data.Dataset.from_tensor_slices(
            (temp_images.astype(np.float32), temp_labels.astype(np.int32))
        )
        dataset = dataset.shuffle(100).batch(batch_size).repeat(repetitions)
        if split:
            return dataset, test_images, test_labels
        return dataset

train_dataset = Dataset(training=True)
test_dataset = Dataset(training=False)

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

model.compile()
# optimizer = keras.optimizers.SGD(learning_rate=learning_rate)
optimizer = tf.optimizers.SGD()

"""
## Train the model
"""

training = []
testing = []
for meta_iter in range(meta_iters):
    frac_done = meta_iter / meta_iters
    cur_meta_step_size = (1 - frac_done) * meta_step_size
    # Temporarily save the weights from the model.
    old_vars = model.get_weights()
    # Get a sample from the full dataset. Train the few-shot.
    mini_dataset = train_dataset.get_mini_dataset(
        inner_batch_size, inner_iters, train_shots, classes
    )

    for images, labels in mini_dataset:
        with tf.GradientTape() as tape:
            preds = model(images)
            loss = keras.losses.sparse_categorical_crossentropy(labels, preds)
        grads = tape.gradient(loss, model.trainable_weights)
        optimizer.apply_gradients(zip(grads, model.trainable_weights))
    new_vars = model.get_weights()
    # Perform SGD for the meta step.
    for var in range(len(new_vars)):
        new_vars[var] = old_vars[var] + (
            (new_vars[var] - old_vars[var]) * cur_meta_step_size
        )
    # After the meta-learning step, reload the newly-trained weights into the model.
    model.set_weights(new_vars)

    # Evaluation loop
    if meta_iter % eval_interval == 0:
        accuracies = []
        for dataset in (train_dataset, test_dataset):
            # Sample a mini dataset from the full dataset.
            train_set, test_images, test_labels = dataset.get_mini_dataset(
                eval_batch_size, eval_iters, shots, classes, split=True
            )
            old_vars = model.get_weights()
            # Train on the samples and get the resulting accuracies.
            for images, labels in train_set:
                with tf.GradientTape() as tape:
                    preds = model(images)
                    loss = keras.losses.sparse_categorical_crossentropy(labels, preds)
                grads = tape.gradient(loss, model.trainable_weights)
                optimizer.apply_gradients(zip(grads, model.trainable_weights))
            test_preds = model.predict(test_images)
            test_preds = tf.argmax(test_preds).numpy()
            num_correct = (test_preds == test_labels).sum()
            # Reset the weights after getting the evaluation accuracies.
            model.set_weights(old_vars)
            accuracies.append(num_correct / classes)

        training.append(accuracies[0])
        testing.append(accuracies[1])



        if meta_iter % 10 == 0:
            print(
                "batch %d: train=%f test=%f" % (meta_iter, accuracies[0], accuracies[1])
            )

        if meta_iter % 100 == 0:
            preModelName = "model_{}.h5".format(meta_iter-100)
            modelName = "model_{}.h5".format(meta_iter)
            model.save(modelName)
            if(os.path.exists(preModelName)):
                os.remove(preModelName)

        if meta_iter % 1999 == 0:
            preModelName = "model_{}.h5".format(meta_iter-99)
            modelName = "model_{}.h5".format(meta_iter)
            model.save(modelName)
            if(os.path.exists(preModelName)):
                os.remove(preModelName)

# First, some preprocessing to smooth the training and testing arrays for display.
window_length = 100
train_s = np.r_[
    training[window_length - 1: 0: -1], training, training[-1:-window_length:-1]
]
test_s = np.r_[
    testing[window_length - 1: 0: -1], testing, testing[-1:-window_length:-1]
]
w = np.hamming(window_length)
train_y = np.convolve(w / w.sum(), train_s, mode="valid")  # Numpy 中的卷积函数
test_y = np.convolve(w / w.sum(), test_s, mode="valid")

# Display the training accuracies.
x = np.arange(0, len(test_y), 1)
savemat(file_name, {'x': x, 'test_y': test_y, 'train_y': train_y})

train_set, test_images, test_labels = dataset.get_mini_dataset(
    eval_batch_size, eval_iters, shots, classes, split=True
)

# Fetch a mini-set to train
for images, labels in train_set:
    with tf.GradientTape() as tape:
        preds = model(images)
        loss = keras.losses.sparse_categorical_crossentropy(labels, preds)
    grads = tape.gradient(loss, model.trainable_weights)
    optimizer.apply_gradients(zip(grads, model.trainable_weights))
test_preds = model.predict(test_images)
test_preds = tf.argmax(test_preds).numpy()

sample_keys = list(train_dataset.data.keys())

for i in range(len(test_preds)):
    print(
        "Label : {}, Prediction : {}".format(int(test_labels[i]), test_preds[i])
    )

