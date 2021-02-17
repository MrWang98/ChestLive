from keras import layers
from keras.layers import Input, Dense, Activation, ZeroPadding2D, BatchNormalization, Flatten, Conv2D, Conv1D, ZeroPadding1D
from keras.layers import AveragePooling2D, MaxPooling2D, Dropout, GlobalMaxPooling2D, GlobalAveragePooling2D, MaxPooling1D
from keras.models import Model
from keras.preprocessing import image
from keras.utils import layer_utils
from keras.utils.data_utils import get_file
from keras.applications.imagenet_utils import preprocess_input
import pydot
from IPython.display import SVG
from keras.utils.vis_utils import model_to_dot
from keras.utils import plot_model
# from kt_utils import *
import keras.backend as K
K.set_image_data_format('channels_last')
import matplotlib.pyplot as plt
from matplotlib.pyplot import imshow


def VoiceChestModel_Keras(input_shape):
    X_input = Input(input_shape)

    # X = ZeroPadding1D(3)(X_input)

    # X = Conv2D(32, (7, 7), strides=(1, 1), name='conv0')(X_input)
    X = Conv1D(32, 4, strides=1, name='conv0')(X_input) #filter setting. 3
    # X = BatchNormalization(axis=1, name='bn0')(X)
    X = Activation('relu')(X)
    X = Dropout(0.3)(X)
    X = MaxPooling1D((4), name='max_pool0')(X)

    X = Conv1D(64, 4, strides=1, name='conv1')(X)  # filter setting.
    X = Activation('relu')(X)
    X = Dropout(0.3)(X)
    X = MaxPooling1D((4), name='max_pool01')(X)

    X = Conv1D(128, 4, strides=1, name='conv11')(X)  # filter setting.
    X = Activation('relu')(X)
    X = Dropout(0.3)(X)
    X = MaxPooling1D((4), name='max_pool1')(X)

    X = Conv1D(256, 4, strides=1, name='conv2')(X)  # filter setting.
    X = Activation('relu')(X)
    X = Dropout(0.3)(X)
    # X = MaxPooling1D((4), name='max_pool11')(X)

    # X = Conv1D(256, 4, strides=1, name='conv12')(X)  # filter setting.
    # X = Activation('relu')(X)
    # X = Dropout(0.3)(X)
    X = MaxPooling1D((2), name='max_pool')(X)

    X = Flatten()(X)
    # X = Dropout(0.2)(X)
    X = Dense(15, activation='softmax', name='fc')(X)  # 1 softmax sigmoid

    model = Model(inputs=X_input, outputs=X, name='VoiceChestModel_Keras')

    return model


