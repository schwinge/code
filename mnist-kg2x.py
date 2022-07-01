import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)

import os
for dirname, _, filenames in os.walk('/kaggle/input'):
    for filename in filenames:
        print(os.path.join(dirname, filename))  #pwd input

import tensorflow as tf
from tensorflow.keras.layers import (
  Conv2D, MaxPool2D, Flatten, Dense)
from tensorflow.keras import Model
from matplotlib import pyplot as plt
import numpy as np


train = pd.read_csv("/kaggle/input/digit-recognizer/train.csv")
test = pd.read_csv("/kaggle/input/digit-recognizer/test.csv")
x_tmp = train.drop("label",1)
y_tmp = train.label


#from sklearn.model_selection import train_test_split
#x_train, x_valid, y_train, y_valid = train_test_split(x_train , y_train , test_size=0.2 , random_state=40 )
#simply use fixed epoch then
#x_train, x_valid, test = x_train / 255.0, x_valid / 255.0, test / 255.0
x_train, test = x_tmp / 255.0, test / 255.0
x_train.shape  #2D matrix
y_train = y_tmp
img_size = 28

#x_train = x_train[:1000].reshape(1000,28,28,1)
#2 NHWC, need drop y first
x_train = x_train.values.reshape(-1,img_size,img_size,1)
#x_valid = x_valid.values.reshape(-1,img_size,img_size,1)
train_ds = tf.data.Dataset.from_tensor_slices(
    (x_train, y_train)).shuffle(10000).batch(128)
plt.imshow(np.array(x_train[1]))
plt.show()


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

