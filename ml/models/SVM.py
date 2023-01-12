# Author: Brenton Graham
# Description: SVM classifier with hyperparameter tuning
# Last updated: 12/29/2022


import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold, GridSearchCV


class SVM():
    '''
    Class to tune, train and test an SVM classifier
    '''
    
    def __init__(self, kfold=5):
        
        # Hyperparameters to tune
        self.grid = [
            {'C': [0.01, 0.1, 1, 10, 100], 'kernel': ['linear'], 'class_weight': ['balanced', None]}, 
            {'C': [0.01, 0.1, 1, 10, 100], 'kernel': ['rbf'], 'gamma': [0.001, 0.0001], 'class_weight': ['balanced', None]}]
        
        # Specify model with randomized grid search
        self.model = GridSearchCV(
            estimator=SVC(probability=True, cache_size=10000), param_grid=self.grid, 
            cv=StratifiedKFold(kfold, shuffle=True), scoring='roc_auc', error_score='raise', n_jobs=-1)
        
        # Set model attributes
        self.best_params_ = {}
        
        
    def train(self, trainingData: list):
        '''
        Function to train SVM classifier. Each training instance will perform a grid search of
        hyperparameters using the provided training data.
        
        Input
        - training_data: list in the form of [train_x, train_y]
        
        '''
        # Fit model on training data
        x, y = trainingData
        self.model.fit(x, y)
        
        # Choose best model and update attributes
        self.best_params_ = self.model.best_params_
        self.model = self.model.best_estimator_
        
        
    def predict(self, x, proba: bool):
        '''
        Function to test trained SVM classifier on a test set.
        
        Input
        - x: feature set of new observations to classify
        - proba: flag to return probabilities rather than class prediction
        '''
        
        # Make predictions on test set
        preds = self.model.predict(x)
        probs = self.model.predict_proba(x)[:,1]
        
        # Return predictions (or probabilities)
        return probs if proba else preds
    
    