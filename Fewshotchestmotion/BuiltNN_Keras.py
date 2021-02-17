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

    X = Conv1D(32, 3, strides=1, name='conv0')(X_input) #filter setting.
    X = Activation('relu')(X)
    X = Dropout(0.2)(X)

    X = MaxPooling1D((2), name='max_pool')(X)
    X = Flatten()(X)
    X = Dropout(0.2)(X)
    X = Dense(2, activation='softmax', name='fc')(X)  # 1 softmax sigmoid 2

    model = Model(inputs=X_input, outputs=X, name='VoiceChestModel_Keras')

    return model


