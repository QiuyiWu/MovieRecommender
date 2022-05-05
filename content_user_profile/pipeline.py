"""
Module that handls the reading
"""
import pandas as pd
import numpy as np

from functools import lru_cache

class ContentBasedModelPipeline:
    def __init__(self, random_seed=None):
        self.rng = np.random.RandomState(random_seed)
    
    @lru_cache(maxsize=1)
    def movie_features(self):
        """get above movies and their genome tags & scores"""
        return pd.read_feather('./data/movie_profiles.feather').set_index('movieId')
    
    @lru_cache(maxsize=1)
    def user_hist(self):
        return pd.read_feather('./data/user_hist.feather').set_index(['userId','movieId'])
    
    @lru_cache(maxsize=1)
    def user_profiles(self):
        return pd.read_feather('./data/user_profiles.feather').set_index('userId')

    def get_user_train_test(self, user_id):
        def calc_feature(user_profile, movie_profile):
            if user_profile is not None:
                x = np.r_[user_profile, movie_profile]
                x /= np.linalg.norm(x)
            else:
                x = None
            return x

        user_profiles = self.user_profiles().loc[user_id]
        user_hist = self.user_hist().loc[user_id]
        data = self.movie_features().copy()

        # add in avaliable ratings
        data['rating'] = user_hist['rating']

        # construct features for train and test
        train, test = {}, {}
        for movieId, genome in data.genome.iteritems():
            train[movieId] = calc_feature(user_profiles.train, genome)
            test[movieId] = calc_feature(user_profiles.test, genome)

        data['train'] = pd.Series(train)
        data['test'] = pd.Series(test)
        return data

    def get_batch_train_test(self, batch_size=5000):
        """Constructing the finalized dataset for model"""
        hist = self.user_hist()
        user_profiles = self.user_profiles()
        movie_features = self.movie_features()

        data = hist.sample(batch_size, random_state=self.rng)
        data = data.join(user_profiles).join(movie_features)

        # concate user profile and movie genome features to one
        # normalize each row, so that euclidean distance emulate cosine distance ordering
        has_train, has_test = True, True

        if data.train.isna().sum() != data.shape[0]:
            data['train'] = (
                data[~data.train.isna()]
                .apply(lambda x: np.r_[x.train, x.genome], axis=1)
                .apply(lambda x: x / np.linalg.norm(x))
            )
        else:
            has_train = False
        
        if data.test.isna().sum() != data.shape[0]:
            data['test'] = (
                data.apply(lambda x: np.r_[x.test, x.genome], axis=1)
                .apply(lambda x: x / np.linalg.norm(x))
            )
        else:
            has_test = False

        # cleanup
        data.drop(['genome'], axis=1, inplace=True)

        # get train & test subsets
        if has_train:
            X_train = np.r_[data[data.is_train].train.to_list()]
            y_train = data[data.is_train].rating.to_numpy()
        else:
            X_train, y_train = None, None

        if has_test:
            X_test = np.r_[data[~data.is_train].test.to_list()]
            y_test = data[~data.is_train].rating.to_numpy()
        else:
            X_test, y_test = None, None

        return X_train, y_train, X_test, y_test
