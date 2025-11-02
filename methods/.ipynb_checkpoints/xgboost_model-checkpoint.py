import joblib
from xgboost import XGBClassifier

class XGBoostModel:
    def __init__(self, n_estimators=100, random_state=42, eval_metric="mlogloss"):
        self.model = XGBClassifier(
            n_estimators=n_estimators,
            random_state=random_state,
            eval_metric=eval_metric
        )

    def train(self, X_train, y_train):
        self.model.fit(X_train, y_train)

    def predict(self, X_test):
        return self.model.predict(X_test)

    def score(self, X_test, y_test):
        return self.model.score(X_test, y_test)

    def save(self, filepath):
        """Save the trained XGBoost model to a file"""
        joblib.dump(self.model, filepath)

    @classmethod
    def load(cls, filepath):
        """Load a saved XGBoost model from file"""
        instance = cls()
        instance.model = joblib.load(filepath)
        return instance

