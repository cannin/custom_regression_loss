from sklearn.model_selection import KFold

# Used to cross-validate models and identify optimal lambda
class CustomCrossValidator:
    """
    Cross validates arbitrary model using MAPE criterion on list of lambdas.
    """
    def __init__(self, X, Y, ModelClass,
                 sample_weights=None,
                 loss_function=mean_absolute_percentage_error):
        self.X = X
        self.Y = Y
        self.ModelClass = ModelClass
        self.loss_function = loss_function
        self.sample_weights = sample_weights


    def cross_validate(self, lambdas, num_folds=10):
        """
        lambdas: set of regularization parameters to try
        num_folds: number of folds to cross-validate against
        """
        self.lambdas = lambdas
        self.cv_scores = []
        X = self.X
        Y = self.Y 
        
        # Beta values are not likely to differ dramatically
        # between differnt folds. Keeping track of the estimated
        # beta coefficients and passing them as starting values
        # to the .fit() operator on our model class can significantly
        # lower the time it takes for the minimize() function to run
        beta_init = None
        
        for lam in self.lambdas:
            print("Lambda: {}".format(lam))
            
            # Split data into training/holdout sets
            kf = KFold(n_splits=num_folds, shuffle=True)
            kf.get_n_splits(X)
            
            # Keep track of the error for each holdout fold
            k_fold_scores = []
            
            # Iterate over folds, using k-1 folds for training
            # and the k-th fold for validation
            f = 1
            for train_index, test_index in kf.split(X):
                # Training data
                CV_X = X[train_index,:]
                CV_Y = Y[train_index]
                CV_weights = None
                if type(self.sample_weights) != type(None):
                    CV_weights = self.sample_weights[train_index]
                
                # Holdout data
                holdout_X = X[test_index,:]
                holdout_Y = Y[test_index]
                holdout_weights = None
                if type(self.sample_weights) != type(None):
                    holdout_weights = self.sample_weights[test_index]
                
                # Fit model to training sample
                lambda_fold_model = self.ModelClass(
                    regularization=lam,
                    X=CV_X,
                    Y=CV_Y,
                    sample_weights=CV_weights,
                    beta_init=beta_init,
                    loss_function=self.loss_function
                )
                lambda_fold_model.fit()
                
                # Extract beta values to pass as beta_init 
                # to speed up estimation of the next fold
                beta_init = lambda_fold_model.beta
                
                # Calculate holdout error
                fold_preds = lambda_fold_model.predict(holdout_X)
                fold_mape = mean_absolute_percentage_error(
                    holdout_Y, fold_preds, sample_weights=holdout_weights
                )
                k_fold_scores.append(fold_mape)
                print("Fold: {}. Error: {}".format( f, fold_mape))
                f += 1
            
            # Error associated with each lambda is the average
            # of the errors across the k folds
            lambda_scores = np.mean(k_fold_scores)
            print("LAMBDA AVERAGE: {}".format(lambda_scores))
            self.cv_scores.append(lambda_scores)
        
        # Optimal lambda is that which minimizes the cross-validation error
        self.lambda_star_index = np.argmin(self.cv_scores)
        self.lambda_star = self.lambdas[self.lambda_star_index]
        print("\n\n**OPTIMAL LAMBDA: {}**".format(self.lambda_star))

# User must specify lambdas over which to search
lambdas = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001]

cross_validator = CustomCrossValidator(
    X, Y, CustomLinearModel,
    loss_function=mean_absolute_percentage_error
)
cross_validator.cross_validate(lambdas, num_folds=5)