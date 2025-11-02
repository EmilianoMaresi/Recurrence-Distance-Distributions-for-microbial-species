#import joblib
from sklearn.ensemble import RandomForestClassifier

class RandomForestModel:
    
    def __init__(self, **kwargs):
        self.model = RandomForestClassifier(**kwargs)

    def train(self, X_train, y_train):
        self.model.fit(X_train, y_train)

    def predict(self, X_test):
        return self.model.predict(X_test)

    def score(self, X_test, y_test):
        return self.model.score(X_test, y_test)

    """
    def save(self, filepath):
        #Save the trained model to a file
        joblib.dump(self.model, filepath)

    @classmethod
    def load(cls, filepath):
        #Load a model from a file and return an instance of RandomForestModel
        instance = cls()              # create a new instance
        instance.model = joblib.load(filepath)  # replace its model with loaded one
        return instance
    """
