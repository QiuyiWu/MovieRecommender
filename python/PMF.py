# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 10:12:09 2022
Origin: https://github.com/fuhailin/Probabilistic-Matrix-Factorization

@author: mingy
"""

import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

#%% PMF
class PMF(object):
    def __init__(self, num_feat=10, epsilon=1, _lambda=0.01, momentum=0.8, maxepoch=50, num_batches=9, batch_size=1000):
        #******** Assign Parameters ******
        self.num_feat = num_feat  # Number of latent features,
        self.epsilon = epsilon  # learning rate,
        self._lambda = _lambda  # L2 regularization,
        self.momentum = momentum  # momentum of the gradient,
        self.maxepoch = maxepoch  # Number of epoch before stop,
        self.num_batches = num_batches  # Number of batches in each epoch (for SGD optimization),
        self.batch_size = batch_size  # Number of training samples used in each batches (for SGD optimization)

        self.w_Item = None  # Item feature vectors
        self.w_User = None  # User feature vectors

        self.rmse_train = []
        self.rmse_test = []

    # ***Fit the model with train_tuple and evaluate RMSE on both train and test data.  ***********#
    # ***************** train_vec=TrainData, test_vec=TestData*************#
    def fit(self, train_vec, test_vec):
        # mean subtraction
        self.mean_inv = np.mean(train_vec[:, 2])  # Mean Ratings

        pairs_train = train_vec.shape[0]
        pairs_test = test_vec.shape[0]

        # 1-p-i, 2-m-c
        num_user = int(max(np.amax(train_vec[:, 0]), np.amax(test_vec[:, 0]))) + 1  # number of users
        num_item = int(max(np.amax(train_vec[:, 1]), np.amax(test_vec[:, 1]))) + 1  # number of movies

        incremental = False
        if ((not incremental) or (self.w_Item is None)): # run when there is incrememntal? run nayways by default
            # initialize
            self.epoch = 0
            self.w_Item = 0.1 * np.random.randn(num_item, self.num_feat)  # numpy.random.randn uniform M x D movies
            self.w_User = 0.1 * np.random.randn(num_user, self.num_feat)  # numpy.random.randn uniform N x D users

            self.w_Item_inc = np.zeros((num_item, self.num_feat))  # M x D
            self.w_User_inc = np.zeros((num_user, self.num_feat))  # N x D

        while self.epoch < self.maxepoch:  # Loop through 1:maxepoch
            self.epoch += 1

            # Shuffle training triples
            shuffled_order = np.arange(train_vec.shape[0])  # creaste ordered and evenly distributed id for traning date
            np.random.shuffle(shuffled_order)  # Randomly permute a sequence for traning data id

            # Batch update
            for batch in range(self.num_batches):  # for each batch
                print("epoch %d batch %d " %(self.epoch, batch+1))

                test = np.arange(self.batch_size * batch, self.batch_size * (batch + 1))
                batch_idx = np.mod(test, shuffled_order.shape[0])  # 本次迭代要使用的索引下标

                batch_UserID = np.array(train_vec[shuffled_order[batch_idx], 0], dtype='int32')
                batch_ItemID = np.array(train_vec[shuffled_order[batch_idx], 1], dtype='int32')

                # Compute Objective Function
                pred_out = np.sum(np.multiply(self.w_User[batch_UserID, :],
                                              self.w_Item[batch_ItemID, :]),
                                  axis=1)  # mean_inv subtracted # np.multiply对应位置元素相乘

                rawErr = pred_out - train_vec[shuffled_order[batch_idx], 2] + self.mean_inv

                # Compute gradients
                Ix_User = 2 * np.multiply(rawErr[:, np.newaxis], self.w_Item[batch_ItemID, :]) \
                       + self._lambda * self.w_User[batch_UserID, :]
                Ix_Item = 2 * np.multiply(rawErr[:, np.newaxis], self.w_User[batch_UserID, :]) \
                       + self._lambda * (self.w_Item[batch_ItemID, :])  # np.newaxis :increase the dimension

                dw_Item = np.zeros((num_item, self.num_feat))
                dw_User = np.zeros((num_user, self.num_feat))

                # loop to aggreate the gradients of the same element
                for i in range(self.batch_size):
                    dw_Item[batch_ItemID[i], :] += Ix_Item[i, :]
                    dw_User[batch_UserID[i], :] += Ix_User[i, :]

                # Update with momentum
                self.w_Item_inc = self.momentum * self.w_Item_inc + self.epsilon * dw_Item / self.batch_size
                self.w_User_inc = self.momentum * self.w_User_inc + self.epsilon * dw_User / self.batch_size

                self.w_Item = self.w_Item - self.w_Item_inc
                self.w_User = self.w_User - self.w_User_inc

                # Compute Objective Function after
                if batch == self.num_batches - 1:
                    pred_out = np.sum(np.multiply(self.w_User[np.array(train_vec[:, 0], dtype='int32'), :],
                                                  self.w_Item[np.array(train_vec[:, 1], dtype='int32'), :]),
                                      axis=1)  # mean_inv subtracted
                    rawErr = pred_out - train_vec[:, 2] + self.mean_inv
                    obj = np.linalg.norm(rawErr) ** 2 \
                          + 0.5 * self._lambda * (np.linalg.norm(self.w_User) ** 2 + np.linalg.norm(self.w_Item) ** 2)

                    self.rmse_train.append(np.sqrt(obj / pairs_train))

                # Compute validation error
                if batch == self.num_batches - 1:
                    pred_out = np.sum(np.multiply(self.w_User[np.array(test_vec[:, 0], dtype='int32'), :],
                                                  self.w_Item[np.array(test_vec[:, 1], dtype='int32'), :]),
                                      axis=1)  # mean_inv subtracted
                    rawErr = pred_out - test_vec[:, 2] + self.mean_inv
                    self.rmse_test.append(np.linalg.norm(rawErr) / np.sqrt(pairs_test))

                    # Print info
                    if batch == self.num_batches - 1:
                        print('Training RMSE: %f, Test RMSE %f' % (self.rmse_train[-1], self.rmse_test[-1]))

    def predict(self, invID):
        return np.dot(self.w_Item, self.w_User[int(invID), :]) + self.mean_inv  # numpy.dot 点乘

    # ****************Set parameters by providing a parameter dictionary.  ***********#
    def set_params(self, parameters):
        if isinstance(parameters, dict):
            self.num_feat = parameters.get("num_feat", 10)
            self.epsilon = parameters.get("epsilon", 1)
            self._lambda = parameters.get("_lambda", 0.01)
            self.momentum = parameters.get("momentum", 0.8)
            self.maxepoch = parameters.get("maxepoch", 50)
            self.num_batches = parameters.get("num_batches", 10)
            self.batch_size = parameters.get("batch_size", 1000)

    def topK(self, test_vec, k=10):
        inv_lst = np.unique(test_vec[:, 0])
        pred = {}
        for inv in inv_lst:
            if pred.get(inv, None) is None:
                pred[inv] = np.argsort(self.predict(inv))[-k:]  # numpy.argsort索引排序

        intersection_cnt = {}
        for i in range(test_vec.shape[0]):
            if test_vec[i, 1] in pred[test_vec[i, 0]]:
                intersection_cnt[test_vec[i, 0]] = intersection_cnt.get(test_vec[i, 0], 0) + 1
        invPairs_cnt = np.bincount(np.array(test_vec[:, 0], dtype='int32'))

        precision_acc = 0.0
        recall_acc = 0.0
        for inv in inv_lst:
            precision_acc += intersection_cnt.get(inv, 0) / float(k)
            recall_acc += intersection_cnt.get(inv, 0) / float(invPairs_cnt[int(inv)])

        return precision_acc / len(inv_lst), recall_acc / len(inv_lst)

#%% Load Data Function
def load_rating_data(folder_path='D:/PhD/Recommender System/PMF/ml-25m/'):
    """
    load movie lens 25m ratings from original rating file.
    """
    mlist = pd.read_csv(folder_path + 'genome-scores.csv', usecols = ['movieId']).squeeze().unique() # List of unique movies that have genome scores
    prefer = []
    for line in open(folder_path + 'ratings.csv', 'r').readlines()[1:]:  # opening "ratings file"
        (userid, movieid, rating, ts) = line.split(',') 
        uid = int(userid)
        mid = int(movieid)
        rat = float(rating)
        if mid in mlist:
            prefer.append([uid, mid, rat])
    data = np.array(prefer)
    return data

#%% Spliting Training and Testing
def spilt_rating_dat(data, size=0.2):
    train_data = []
    test_data = []
    for line in data:
        rand = random.random()
        if rand < size:
            test_data.append(line)
        else:
            train_data.append(line)
    train_data = np.array(train_data)
    test_data = np.array(test_data)
    return train_data, test_data

#%% Run Selected Dataset with 24,674,113 entries
if __name__ == "__main__":
    folder_path='D:/PhD/Recommender System/PMF/ml-25m/'
    pmf = PMF()
    pmf.set_params({"num_feat": 10, "epsilon": 1, "_lambda": 0.01, "momentum": 0.8, "maxepoch": 50, "num_batches": 2500,
                    "batch_size": 1000})
    ratings = load_rating_data(folder_path)
    print(len(np.unique(ratings[:, 0])), len(np.unique(ratings[:, 1])), pmf.num_feat) # 162540 13816 10
    train, test = train_test_split(ratings, test_size=0.2)  # spilt_rating_dat(ratings)
    pmf.fit(train, test)

    # Check performance by plotting train and test errors
    plt.plot(range(pmf.maxepoch), pmf.rmse_train, marker='o', label='Training Data')
    plt.plot(range(pmf.maxepoch), pmf.rmse_test, marker='v', label='Test Data')
    plt.title('The MovieLens Dataset Learning Curve')
    plt.xlabel('Number of Epochs')
    plt.ylabel('RMSE')
    plt.legend()
    plt.grid()
    plt.show()
    print("precision_acc,recall_acc:" + str(pmf.topK(test)))

#%% 