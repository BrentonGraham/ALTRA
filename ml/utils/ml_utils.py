# Author: Brenton Graham
# Description: Common functions used in ML scripts
# Last updated: 01/11/2023


import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, matthews_corrcoef, precision_score, recall_score, f1_score


def xy_split(df, target: str, to_numpy: bool):
    '''
    Function to split df into features (x) and target (y)
    
    Input
    - target: string corresponding to the target column name
    - to_numpy: boolean flag that can be used to convert x and y outputs to numpy arrays
    '''
    
    # Split data frame
    x = df.loc[:, df.columns != target]
    y = df.loc[:, target]
    
    # Return x, y in specified format
    return (x.to_numpy(), y.to_numpy()) if to_numpy else (x, y)


def flatten_list(nested_list):
    '''
    Function to flatten a list of lists into one list
    e.g. [[1, 0, 0], [0, 0, 1]] -> [1, 0, 0, 0, 0, 1]
    Source: https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
    '''
    
    return [item for sublist in nested_list for item in sublist]


def numpy_to_DF(numpyMatrix, columns):
    '''
    Function to convert numpy matrix to data frame
    
    Input
    - numpyMatrix
    - columns: column names, usually extracted from original df (e.g. columns = x_df.columns)
    '''
    
    return pd.DataFrame(numpyMatrix, columns=columns)


def getClfStats(y_true, y_pred, y_prob):
    '''
    Function to evaluate classifier metrics, including AUC, Precision, Recall, and F1 score
    '''
    
    statsDict = {}
    statsDict["AUC"] = round(roc_auc_score(y_true, y_prob), 3)
    statsDict["Precision"] = round(precision_score(y_true, y_pred), 3)
    statsDict["Recall"] = round(recall_score(y_true, y_pred), 3)
    statsDict["F1"] = round(f1_score(y_true, y_pred), 3)
    return statsDict
