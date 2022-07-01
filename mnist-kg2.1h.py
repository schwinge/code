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

#cross-entropy losses to categorical y 
train = pd.read_csv('/kaggle/input/digit-recognizer/train.csv')
test = pd.read_csv('/kaggle/input/digit-recognizer/test.csv')
x_tmp = train.drop('label',1)
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
test = test.values.reshape(-1,img_size,img_size,1)
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
    self.pool1 = MaxPool2D(strides=1)

    self.conv21 = Conv2D(256, (3,3), activation='relu')
    self.conv22 = Conv2D(256, (3,3), activation='relu', padding = 'same')
    self.pool2 = MaxPool2D(strides=1)

#  self.conv31 = Conv2D(32, (3,3), padding='same')
    self.flatten = Flatten()
    self.d1 = Dense(128, activation='relu')
    self.d2 = Dense(32, activation='relu')
    self.d3 = Dense(10)

#  @property
#  def metrics(self):
#  return [loss_tracker, mae_metric]


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

#optimizer help minimize loss function 
vgg = VGG16()
optimizer = tf.keras.optimizers.Adam()
loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)

#autograph: eager exe env speed up, & graph format exportable
#model.compile
@tf.function
def train_step(images, labels):
  with tf.GradientTape() as tape:
    predictions = vgg(images, training=True)
    loss = loss_object(labels, predictions)
  gradients = tape.gradient(loss, vgg.trainable_variables)
  optimizer.apply_gradients(zip(gradients, vgg.trainable_variables))


EPOCHS = 5
#acc = np.zeros((EPOCHS, 2))
#lss = np.zeros((EPOCHS, 2))


for epoch in range(EPOCHS):
  for images, labels in train_ds:
    train_step(images, labels)


#vgg.compile(optimizer='adam', loss='mse', metrics=['mae'])
predictions = vgg.predict(test)
y = np.argmax(predictions,axis = 1)
#2 yeah! <3
sub = pd.read_csv('/kaggle/input/digit-recognizer/sample_submission.csv')
sub['Label'] = y
sub.to_csv("vggcon3.csv", index=False)
sub.head()

for i in range(10): 
  print(y[i])
  plt.imshow(np.array(test[i]))
  plt.show() 

vgg.summary()
#print(vgg.layers[0].trainable_weights)
#plot_ tfjs


if False: '''
for i in range(10): 
  print(predictions[i])
  plt.imshow(np.array(test[i]))
  plt.show() 
#2 yeah! <3
'''

