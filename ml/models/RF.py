# Author: Brenton Graham
# Description: Random forest classifier with hyperparameter tuning
# Last updated: 12/29/2022


import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV


class RF():
    '''
    Class to tune, train and test a Random Forest classifier
    
    TO DO:
    - Add regressor functionality
    - We could introduce an internal feature selection step based on feature importance
    '''
    
    def __init__(self, kfold=5):
        
        # Hyperparameters to tune
        self.grid = [{
            'n_estimators': [100, 500, 1000], 
            'max_depth': [5, 10, 50, 100], 
            'max_features': ['sqrt', 'log2'], 
            'min_samples_leaf': [1, 2, 4], 
            'class_weight': ['balanced', None]}]
        
        # Specify model with randomized grid search
        self.model = RandomizedSearchCV(
            estimator=RandomForestClassifier(n_jobs=-1), param_distributions=self.grid, n_iter=10,
            cv=StratifiedKFold(kfold, shuffle=True), scoring='roc_auc', error_score='raise', n_jobs=-1)
        
        # Set model attributes
        self.best_params_ = {}
        self.feature_importances_ = []
        
        
    def train(self, trainingData: list):
        '''
        Function to train RF classifier. Each training instance will perform a randomized grid search of
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
        self.feature_importances_ = self.model.feature_importances_
        
        
    def predict(self, x, proba: bool):
        '''
        Function to test trained RF classifier on a test set.
        
        Input
        - x: feature set of new observations to classify
        - proba: flag to return probabilities rather than class prediction
        '''
        
        # Make predictions on test set
        preds = self.model.predict(x)
        probs = self.model.predict_proba(x)[:,1]
        
        # Return predictions (or probabilities)
        return probs if proba else preds
    
    