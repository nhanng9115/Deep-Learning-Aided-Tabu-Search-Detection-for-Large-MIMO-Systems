# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 22:17:21 2020

@author: nhannguyen
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 11:58:21 2019

@author: Nguyen Nhan
"""


import mimo_simple
import numpy as np
import tensorflow as tf


#import matplotlib.pyplot as plot
import scipy.io as spio

tf.reset_default_graph()
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                                Parameters
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# System parameters
Nt = 2
Nr = 4
L = 4
n_epoch_train = 100
mod_scheme = "QPSK"
SNR_vec = [0, 2, 4, 6, 8, 10, 12] # dB
    
N = 2*Nt

# Training parameters
start_learning_rate = 0.0001
decay_factor = 0.97
decay_step_size = 100
train_batch_size = 100



'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                        The architecture of DetNet
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
def neuron_layer(x, n_neuron):
    n_input = int(x.get_shape()[1])
    W = tf.matmul(tf.Variable(tf.random_normal([n_input, n_neuron], stddev=0.01)), tf.eye(n_neuron))
    b = tf.Variable(tf.random_normal([1, n_neuron], stddev=0.01))
    y = tf.matmul(x, W) + b
    return y
    
def network(x, hy, hh):
    batch_size = tf.shape(hy)[0]
    LOSS = []

    x_est = tf.zeros([batch_size, N])
        
    t = tf.Variable(0.5)
    alpha = tf.constant(1.0)
    
    for i in range(1,L):
        hhx = tf.squeeze(tf.matmul(tf.expand_dims(x_est, 1), hh), 1) 
        
        Z = neuron_layer(hy - hhx, N) + neuron_layer(x_est, N)
        beta = tf.constant(0.5)


        if mod_scheme == "QPSK":
            x_est = -1 + tf.nn.relu(Z+t)/tf.abs(t) - tf.nn.relu(Z-t)/tf.abs(t)

        else:        
            x_est = -3 + (tf.nn.relu(Z+2+t) - tf.nn.relu(Z+2-t) \
                      + tf.nn.relu(Z+t) - tf.nn.relu(Z-t) \
                      + tf.nn.relu(Z-2+t) - tf.nn.relu(Z-2-t))/tf.abs(t)
                
        dis = tf.reduce_mean(tf.reduce_mean(tf.square(x - x_est), 1))
        x_est_norm = tf.nn.l2_normalize(x_est, 1)        
        x_norm = tf.nn.l2_normalize(x, 1)
        cos = 1 - tf.abs(tf.reduce_mean(tf.reduce_mean(tf.multiply(x_est_norm, x_norm), 1)))
        LOSS.append(np.log(i)*( alpha*dis + beta*cos ))        
        
    return x_est, LOSS#, BER
 
with tf.device('/gpu:0'):
    HY = tf.placeholder(tf.float32, shape=[None,N], name="HY")
    X = tf.placeholder(tf.float32, shape=[None,N], name="X")
    HH = tf.placeholder(tf.float32, shape=[None,N,N], name="HH")
    
    x_hat, LOSS = network(X, HY, HH)
    x_hat = tf.multiply(x_hat,1,name="x_hat")

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                        Optimization
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
TOTAL_LOSS = tf.add_n(LOSS)
global_step = tf.Variable(0, trainable=False)
learning_rate = tf.train.exponential_decay(start_learning_rate, global_step, decay_step_size, decay_factor, staircase=True)
optimizer = tf.train.AdamOptimizer(learning_rate)
training_op = optimizer.minimize(TOTAL_LOSS)
    
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                            Training DetNet
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
init = tf.global_variables_initializer()
saver = tf.train.Saver(max_to_keep=20)

with tf.Session(config=tf.ConfigProto(allow_soft_placement=True, log_device_placement=True)) as sess:
  
  sess.run(init)
  for ss in range(len(SNR_vec)):
      SNR_dB = SNR_vec[ss]
      DELTA_SNR = 1
      print('SNR_dB = ', SNR_dB, 'dB')
      
      for epoch in range(n_epoch_train + 1):
          with tf.device('/gpu:0'):
              batch_HY, batch_HH, batch_X, _, _, _ = mimo_simple.gen_data(train_batch_size, Nt, Nr, SNR_dB, mod_scheme, DELTA_SNR, 0)
              feed_dict = {HY: batch_HY, HH: batch_HH, X: batch_X}
          
              sess.run(training_op, feed_dict)
          
          if epoch == 10:
              
              with tf.device('/gpu:0'):
                  loss, s_hat = sess.run([LOSS[-1], x_hat], feed_dict)

              BER = mimo_simple.compute_BER(s_hat, batch_X, Nt, mod_scheme)
              print('Training epoch : ', epoch, ', loss : ', loss, '=====> BER = ', BER)
