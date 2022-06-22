if False: '''
class VGG16(Model):
  def __init__(self):
    super(VGG16, self).__init__()
    self.conv11 = Conv2D(32, (3,3), activation='relu')
    self.conv12 = Conv2D(64, (3,3), activation='relu')
    self.pool1 = MaxPool2D(strides=1)

    self.conv21 = Conv2D(256, (3,3), activation='relu')
    self.conv22 = Conv2D(256, (3,3), activation='relu')
    self.pool2 = MaxPool2D(strides=1)

    self.conv31 = Conv2D(256, (3,3), activation='relu')
    self.conv32 = Conv2D(256, (3,3), activation='relu')
    self.conv33 = Conv2D(256, (3,3), activation='relu')
    self.pool3 = MaxPool2D(strides=1)

    self.conv41 = Conv2D(512, (3,3), activation='relu')
    self.conv42 = Conv2D(512, (3,3), activation='relu')
    self.conv43 = Conv2D(512, (3,3), activation='relu')
    self.pool4 = MaxPool2D(strides=1)

    self.conv51 = Conv2D(1, (3,3), activation='relu')
    self.conv52 = Conv2D(1, (3,3), activation='relu')
    self.conv53 = Conv2D(1, (3,3), activation='relu')
    self.pool5 = MaxPool2D(strides=1)

    self.flatten = Flatten()
    self.d1 = Dense(128, activation='relu')
    self.d2 = Dense(64, activation='relu')
    self.d3 = Dense(10)


  def call(self, x):
    x = self.conv11(x)
    x = self.conv12(x)
    x = self.pool1(x)

    x = self.conv21(x)
    x = self.conv22(x)
    x = self.pool2(x)

    x = self.conv31(x)
    x = self.conv32(x)
    x = self.conv33(x)
    x = self.pool3(x)

    x = self.conv41(x)
    x = self.conv42(x)
    x = self.conv43(x)
    x = self.pool4(x)

    x = self.conv51(x)
    x = self.conv52(x)
    x = self.conv53(x)
    x = self.pool5(x)

    x = self.flatten(x)
    x = self.d1(x)
    x = self.d2(x)
    return self.d3(x)
'''