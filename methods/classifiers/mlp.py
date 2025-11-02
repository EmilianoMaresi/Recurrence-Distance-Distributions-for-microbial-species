#import joblib
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import LabelEncoder

class MLPModel:
    
    def __init__(self, **kwargs):
    	self.model = MLPClassifier(**kwargs)
    	self.label_encoder = None  # will be set during training

    def train(self, X_train, y_train):
        # Encode labels if they are strings
        self.label_encoder = LabelEncoder()
        y_train_encoded = self.label_encoder.fit_transform(y_train)
        self.model.fit(X_train, y_train_encoded)

    def predict(self, X_test):
        y_pred_encoded = self.model.predict(X_test)
        # Decode labels back to original form if encoder exists
        if self.label_encoder:
            return self.label_encoder.inverse_transform(y_pred_encoded)
        return y_pred_encoded

    def score(self, X_test, y_test):
        if self.label_encoder:
            y_test_encoded = self.label_encoder.transform(y_test)
        else:
            y_test_encoded = y_test
        return self.model.score(X_test, y_test_encoded)

    """
    def save(self, filepath):
        #Save the trained MLP model to a file
        joblib.dump(self.model, filepath)

    @classmethod
    def load(cls, filepath):
        #Load a saved MLP model from file
        instance = cls()
        instance.model = joblib.load(filepath)
        return instance
    """
