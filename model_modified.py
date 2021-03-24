# FROM: https://alex.miller.im/posts/linear-model-custom-loss-function-regularization-python/
#%matplotlib inline
#from matplotlib import pyplot as plt
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.preprocessing import StandardScaler
from scipy.optimize import minimize
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

# def objective_function(beta, X, Y):
#     error = loss_function(np.matmul(X,beta), Y)
#     return(error)


def loss_function(y_true, y_pred, sample_weight=None, loss_type='mse'):
    if loss_type == 'mse':
        return(mean_squared_error(y_true, y_pred, sample_weight=None))
    if loss_type == 'rmse':
        return(mean_squared_error(y_true, y_pred, sample_weight=None, squared=False))
    if loss_type == 'mae':
        return(mean_absolute_error(y_true, y_pred, sample_weight=None))

    # y_true = np.array(y_true)
    # y_pred = np.array(y_pred)
    # assert len(y_true) == len(y_pred)
    
    # if np.any(y_true==0):
    #     print("Found zeroes in y_true. MAPE undefined. Removing from set...")
    #     idx = np.where(y_true==0)
    #     y_true = np.delete(y_true, idx)
    #     y_pred = np.delete(y_pred, idx)
    #     if type(sample_weight) != type(None):
    #         sample_weight = np.array(sample_weight)
    #         sample_weight = np.delete(sample_weight, idx)
        
    # if type(sample_weight) == type(None):
    #     return(np.mean(np.abs((y_true - y_pred) / y_true)) * 100)
    # else:
    #     sample_weight = np.array(sample_weight)
    #     assert len(sample_weight) == len(y_true)
    #     return(100/np.sum(sample_weight)*np.dot(
    #             sample_weight, (np.abs((y_true - y_pred) / y_true))
    #     ))
#loss_function = mean_absolute_percentage_error


class CustomLinearModel:
    """
    Linear model: Y = XB, fit by minimizing the provided loss_function
    with L2 regularization
    """
    def __init__(self, loss_function=loss_function, loss_type='mse', 
                 X=None, Y=None, sample_weight=None, beta_init=None, 
                 regularization=0.00012, regularization_type='l2', 
                 include_extra_penalty=False, verbose=False):
        self.beta = None
        self.beta_init = beta_init
        self.sample_weight = sample_weight

        self.loss_type = loss_type
        self.loss_function = loss_function
        self.regularization = regularization
        self.regularization_type = regularization_type
        self.include_extra_penalty = include_extra_penalty

        self.X = X
        self.Y = Y

        self.verbose = verbose
    
    
    def predict(self, X):
        prediction = np.matmul(X, self.beta)
        return(prediction)


    def model_error(self):
        error = self.loss_function(
            y_true=self.Y, 
            y_pred=self.predict(self.X), 
            sample_weight=self.sample_weight,
            loss_type=self.loss_type
        )
        return(error)
    

    # L1, L2: Description: https://towardsdatascience.com/l1-and-l2-regularization-methods-ce25e7fc831c
    def l1_regularized_loss(self, beta):
        self.beta = beta
        
        model_error_value = self.model_error()
        regularization_value = self.regularization*np.sum(np.abs(np.array(self.beta)))
        extra_penalty = 25

        if self.include_extra_penalty:
            return_val = model_error_value + extra_penalty + regularization_value
        else: 
            return_val = model_error_value + regularization_value

        if self.verbose: 
            print(f'REG: {str(regularization_value)}; RETURN: {str(return_val)}')

        return(return_val)


    def l2_regularized_loss(self, beta):
        self.beta = beta

        model_error_value = self.model_error()

        # Original code value
        #regularization_value = sum(self.regularization*np.array(self.beta)**2)
        # NOTE: Reason for no square root outside: https://stats.stackexchange.com/questions/449748/why-does-l-2-norm-regularization-not-have-a-square-root
        regularization_value = self.regularization*np.sum(np.array(self.beta)**2)
        extra_penalty = 25

        if self.include_extra_penalty:
            return_val = model_error_value + extra_penalty + regularization_value
        else: 
            return_val = model_error_value + regularization_value

        if self.verbose: 
            print(f'REG: {str(regularization_value)}; RETURN: {str(return_val)}')

        return(return_val)
    

    def fit(self, maxiter=250):        
        # Initialize beta estimates (you may need to normalize
        # your data and choose smarter initialization values
        # depending on the shape of your loss function)
        if type(self.beta_init) == type(None):
            # set beta_init = 1 for every feature
            self.beta_init = np.array([1]*self.X.shape[1])
        else: 
            # Use provided initial values
            pass
            
        if self.beta != None and all(self.beta_init == self.beta):
            print("Model already fit once; continuing fit with more iterations.")
            

        if self.regularization_type == 'l1':
            res = minimize(self.l1_regularized_loss, self.beta_init,
                           method='BFGS', options={'maxiter': 500})
        elif self.regularization_type == 'l2': 
            res = minimize(self.l2_regularized_loss, self.beta_init,
                           method='BFGS', options={'maxiter': 500})
        else: 
            sys.exit(f'ERROR: Unknown regularization type: {self.regularization_type}')

        self.beta = res.x
        self.beta_init = self.beta

# Set random seed
np.random.seed(seed=123)

# Generate predictors
X_raw = np.random.random(100*9)
X_raw = np.reshape(X_raw, (100, 9))

# Standardize the predictors
scaler = StandardScaler().fit(X_raw)
X = scaler.transform(X_raw)

# Add an intercept column to the model.
X = np.abs(np.concatenate((np.ones((X.shape[0],1)), X), axis=1))

# Define my "true" beta coefficients
beta = np.array([2,6,7,3,5,7,1,2,2,8])

# Y = Xb
Y_true = np.matmul(X,beta)

# Observed data with noise
Y = Y_true*np.exp(np.random.normal(loc=0.0, scale=0.2, size=100))

regularization_type = 'l1'
loss_type = 'mae'
include_extra_penalty = False
verbose = True 

l2_model = CustomLinearModel(
    loss_function=loss_function,
    loss_type=loss_type,
    X=X, Y=Y, 
    regularization=0.01,
    regularization_type=regularization_type,
    include_extra_penalty=include_extra_penalty,
    verbose=verbose
)
l2_model.fit()
l2_model.beta

model = CustomLinearModel(
    loss_function=loss_function,
    loss_type=loss_type,
    X=X, Y=Y, 
    regularization=0,
    regularization_type=regularization_type,
    include_extra_penalty=include_extra_penalty,
    verbose=verbose
)
model.fit()
model.beta

pd.DataFrame({
    "true_beta": beta, 
    "estimated_beta": model.beta,
    "regularized_beta": l2_model.beta
})[["true_beta", "estimated_beta", "regularized_beta"]]

stats.pearsonr(beta, l2_model.beta)
stats.pearsonr(model.beta, l2_model.beta)

# OTHER
oddsratio, pvalue = stats.fisher_exact([[8, 2], [1, 5]])
pvals = np.array([0.002,0.06,0.07,0.003,0.05,0.07,0.001,0.02,0.002,0.08])
reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(pvals, method='fdr_bh')

