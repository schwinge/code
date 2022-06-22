import tensorflow as tf
from tensorflow.keras.layers import (
  Conv2D, MaxPool2D, Flatten, Dense)
from tensorflow.keras import Model
from matplotlib import pyplot as plt
import numpy as np
#from tensorflow.python import keras as ks
#from keras.applications import VGG16
#tf v2.8.2, keras: 2.8.0
#print('keras: %s' % tf.keras.__version__)

#pixel value 0-255/8-bit norm to 0-1, 0 is black 
mnist = tf.keras.datasets.mnist
(x_train, y_train), (x_test, y_test) = mnist.load_data()
x_train, x_test = x_train / 255.0, x_test / 255.0
#print(x_train.shape)


#matrix to NHWC
x_train = x_train[..., tf.newaxis].astype("float32")
x_test = x_test[..., tf.newaxis].astype("float32")
train_ds = tf.data.Dataset.from_tensor_slices(
    (x_train, y_train)).shuffle(10000).batch(32)
test_ds = tf.data.Dataset.from_tensor_slices((x_test, y_test)).batch(32)


#class MyModel(tf.keras.Model), MyModel is self
#initialize state & data
class VGG16(Model):
#  def __init__(self,input_shape=(None,28,28,1), **kwargs):
  def __init__(self):
    super(VGG16, self).__init__()
    self.conv11 = Conv2D(32, (3,3), activation='relu')
    self.conv12 = Conv2D(64, (3,3), activation='relu')
    self.pool1 = MaxPool2D(pool_size=[2, 2], strides=2)

    self.conv21 = Conv2D(256, (3,3), activation='relu')
    self.conv22 = Conv2D(256, (3,3), activation='relu')
    self.pool2 = MaxPool2D(strides=1)

#  self.conv31 = Conv2D(32, (3,3), padding='same')
    self.flatten = Flatten()
    self.d1 = Dense(128, activation='relu')
    self.d2 = Dense(32, activation='relu')
    self.d3 = Dense(10)

#  def call(self, input_tensor):
#    x = self.conv1(input_tensor)
  def call(self, x):
    x = self.conv11(x)
    x = self.conv12(x)
    x = self.pool1(x)

    x = self.conv21(x)
    x = self.conv22(x)
    x = self.pool2(x)

    x = self.flatten(x)
    x = self.d1(x)
    x = self.d2(x)
    return self.d3(x)
#    output = self.d2(x)
#    return output
#return x the tensor from d2 the layer


# Create an oop instance, not ec2
vgg = VGG16()
optimizer = tf.keras.optimizers.Adam()
#optimizer = tf.keras.optimizers.Adam(lr)
loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
train_loss = tf.keras.metrics.Mean(name='train_loss')
train_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(name='train_accuracy')
test_loss = tf.keras.metrics.Mean(name='test_loss')
test_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(name='test_accuracy')


#autograph: eager exe env speed up, & graph format exportable
#model.compile
@tf.function
def train_step(images, labels):
  with tf.GradientTape() as tape:
    predictions = vgg(images, training=True)
    loss = loss_object(labels, predictions)
  gradients = tape.gradient(loss, vgg.trainable_variables)
  optimizer.apply_gradients(zip(gradients, vgg.trainable_variables))
  train_loss(loss)
  train_accuracy(labels, predictions)

@tf.function
def test_step(images, labels):
  predictions = vgg(images, training=False)
  t_loss = loss_object(labels, predictions)
  test_loss(t_loss)
  test_accuracy(labels, predictions)

EPOCHS = 5
acc = np.zeros((EPOCHS, 2))
lss = np.zeros((EPOCHS, 2))


for epoch in range(EPOCHS):
  # Reset the metrics at the start of the next epoch
  train_loss.reset_states()
  train_accuracy.reset_states()
  test_loss.reset_states()
  test_accuracy.reset_states()

  for images, labels in train_ds:
    train_step(images, labels)

  for test_images, test_labels in test_ds:
    test_step(test_images, test_labels)

  #acc = np.append(acc, np.array([[train_accuracy.result(), test_accuracy.result()]]), axis=0) 
  #lss = np.vstack((loss, [train_loss.result(), test_loss.result()]))
  acc[epoch] = np.array([[train_accuracy.result(), test_accuracy.result()]])
  lss[epoch] = np.array([[train_loss.result(), test_loss.result()]])


print(acc)
#plt.plot(acc[EPOCHS:, 0])
plt.plot(acc[:, 0], "k")
plt.plot(acc[:, 1], "r--")
plt.legend(["train accuracy","test_accuracy"])
plt.show()
plt.plot(lss[:, 0], "m")
plt.plot(lss[:, 1], "c--")
plt.legend(["train_loss","test_loss"])
plt.show()


vgg.summary()
#print(vgg.layers[0].trainable_weights)
#plot_ tfjs

if False: '''
references:
https://colab.research.google.com/github/tensorflow/docs/blob/master/site/en/tutorials/quickstart/advanced.ipynb#scrollTo=OZACiVqA8KQV
vgg16
'''
