import tensorflow as tf
import os
os.environ['TF_XLA_FLAGS'] = '--tf_xla_enable_xla_devices'
print(tf.test.is_gpu_available())
