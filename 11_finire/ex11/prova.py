# NOTA: USO PYTHON 3.9.12
import numpy as np
import matplotlib.pyplot as plt

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Activation
from tensorflow.keras import backend as K
from tensorflow.keras.utils import get_custom_objects

# parametri con cui giocare
N_epochs = 30
N_train = 3000
Sigma = 0.1

# target parameters of f(x) = ax^3 + bx^2 + cx + d
a = 3 
b = -2
c = -3
d = 4

# generate training inputs
np.random.seed(0)
x_train = np.random.uniform(-1, 1, N_train) 
x_valid = np.random.uniform(-1, 1, 50)
x_valid.sort()
y_target = a * x_valid*x_valid*x_valid + b*x_valid*x_valid + c*x_valid + d # ideal (target) polinomial function

sigma = Sigma # noise standard deviation
y_train = np.random.normal(a * x_train*x_train*x_train + b*x_train*x_train + c*x_train + d, sigma) 
y_valid = np.random.normal(a * x_valid*x_valid*x_valid + b*x_valid*x_valid + c*x_valid + d, sigma)

# plot validation and target dataset
#plt.plot(x_valid, y_target, label='target')
#plt.scatter(x_valid, y_valid, color='r', label='validation data')
#plt.legend()
#plt.grid(True)
#plt.show()

neurons = [30, 100, 500, 1000]
colors = ['r', 'b', 'g', 'green']

for i in range(0,4):
    K.clear_session()
    # compose the NN model
    model = tf.keras.Sequential()
    # provare ad allargare e poi stringere

    model.add(Dense(neurons[i], input_shape=(1,), activation = 'relu')) # <<<<< numero di neuroni e funzione di attivazione
    model.add(Dense(70, activation = 'relu'))
    model.add(Dense(50, activation = 'relu'))
    model.add(Dense(30, activation = 'relu'))
    model.add(Dense(1, activation = 'relu'))

    # compile the model
    model.compile(optimizer='sgd', loss='mse', metrics=['mse']) # <<<<<<< ottimizzatore, loss, metrics

    # train our model: fit the model using training dataset over 10 epochs of 32 batch size each
    history = model.fit(x=x_train, y=y_train, 
                        batch_size=32, epochs=N_epochs,
                        shuffle=True, # a good idea is to shuffle input before (at each epoch)
                        validation_data=(x_valid, y_valid)
                        )

    # Plot training & validation loss values
    #plt.plot(history.history['loss'])
    #plt.plot(history.history['val_loss'])
    #plt.title('Model loss')
    #plt.ylabel('Loss')
    #plt.xlabel('Epoch')
    #plt.legend(['Train', 'Test'], loc='best')
    #plt.show()
    #plt.savefig("imgs/sgm_"+str(sigma)+".png")

    # plot predictions
    x_predicted = np.random.uniform(-1, 2, 100)
    y_predicted = model.predict(x_predicted) # perchè ha dimensione doppia di x_predicted? perchè devo finire con un solo nodo
    plt.scatter(x_predicted, y_predicted, marker = ".", color=colors[i], label = str(neurons[i]))

plt.grid(True)
plt.plot(x_valid, y_target)
plt.legend()
plt.show()