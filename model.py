# FROM: https://alex.miller.im/posts/linear-model-custom-loss-function-regularization-python/
#%matplotlib inline
#from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from scipy.optimize import minimize

# def objective_function(beta, X, Y):
#     error = loss_function(np.matmul(X,beta), Y)
#     return(error)


def loss_function(y_pred, y_true, sample_weights=None):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    assert len(y_true) == len(y_pred)
    
    if np.any(y_true==0):
        print("Found zeroes in y_true. MAPE undefined. Removing from set...")
        idx = np.where(y_true==0)
        y_true = np.delete(y_true, idx)
        y_pred = np.delete(y_pred, idx)
        if type(sample_weights) != type(None):
            sample_weights = np.array(sample_weights)
            sample_weights = np.delete(sample_weights, idx)
        
    if type(sample_weights) == type(None):
        return(np.mean(np.abs((y_true - y_pred) / y_true)) * 100)
    else:
        sample_weights = np.array(sample_weights)
        assert len(sample_weights) == len(y_true)
        return(100/np.sum(sample_weights)*np.dot(
                sample_weights, (np.abs((y_true - y_pred) / y_true))
        ))
#loss_function = mean_absolute_percentage_error


class CustomLinearModel:
    """
    Linear model: Y = XB, fit by minimizing the provided loss_function
    with L2 regularization
    """
    def __init__(self, loss_function=loss_function, 
                 X=None, Y=None, sample_weights=None, beta_init=None, 
                 regularization=0.00012, regularization_type='l2'):
        self.regularization = regularization
        self.regularization_type = regularization_type
        self.beta = None
        self.loss_function = loss_function
        self.sample_weights = sample_weights
        self.beta_init = beta_init
        
        self.X = X
        self.Y = Y
    
    
    def predict(self, X):
        prediction = np.matmul(X, self.beta)
        return(prediction)


    def model_error(self):
        error = self.loss_function(
            y_pred=self.predict(self.X), 
            y_true=self.Y, 
            sample_weights=self.sample_weights
        )
        return(error)
    

    # L1, L2: Description: https://towardsdatascience.com/l1-and-l2-regularization-methods-ce25e7fc831c
    def l1_regularized_loss(self, beta):
        self.beta = beta
        return(self.model_error() + self.regularization*np.sum(np.abs(np.array(self.beta))))


    def l2_regularized_loss(self, beta):
        self.beta = beta
        #return(self.model_error() + sum(self.regularization*np.array(self.beta)**2))
        # NOTE: Reason for no square root outside: https://stats.stackexchange.com/questions/449748/why-does-l-2-norm-regularization-not-have-a-square-root
        return(self.model_error() + self.regularization*np.sum(np.array(self.beta)**2))
    

    def fit(self, maxiter=250):        
        # Initialize beta estimates (you may need to normalize
        # your data and choose smarter initialization values
        # depending on the shape of your loss function)
        if type(self.beta_init)==type(None):
            # set beta_init = 1 for every feature
            self.beta_init = np.array([1]*self.X.shape[1])
        else: 
            # Use provided initial values
            pass
            
        if self.beta!=None and all(self.beta_init == self.beta):
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

regularization_type = 'l2'

l2_mape_model = CustomLinearModel(
    loss_function=loss_function,
    X=X, Y=Y, 
    regularization=0.01,
    regularization_type=regularization_type
)
l2_mape_model.fit()
l2_mape_model.beta

mape_model = CustomLinearModel(
    loss_function=loss_function,
    X=X, Y=Y, 
    regularization=0,
    regularization_type=regularization_type
)
mape_model.fit()
mape_model.beta

pd.DataFrame({
    "true_beta": beta, 
    "estimated_beta": mape_model.beta,
    "regularized_beta": l2_mape_model.beta
})[["true_beta", "estimated_beta", "regularized_beta"]]
