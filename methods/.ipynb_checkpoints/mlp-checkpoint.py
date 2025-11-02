import joblib
from sklearn.neural_network import MLPClassifier

class MLPModel:
    """
    MLP classifier with increased max_iter, early stopping, and learning rate.
    """
    def __init__(self, hidden_layer_sizes=(100,), max_iter=2000, random_state=42):
        self.model = MLPClassifier(
            hidden_layer_sizes=hidden_layer_sizes,
            max_iter=max_iter,
            early_stopping=True,
            learning_rate_init=0.001,
            random_state=random_state
        )

    def train(self, X_train, y_train):
        self.model.fit(X_train, y_train)

    def predict(self, X_test):
        return self.model.predict(X_test)

    def score(self, X_test, y_test):
        return self.model.score(X_test, y_test)

    def save(self, filepath):
        """Save the trained MLP model to a file"""
        joblib.dump(self.model, filepath)

    @classmethod
    def load(cls, filepath):
        """Load a saved MLP model from file"""
        instance = cls()
        instance.model = joblib.load(filepath)
        return instance

