from datetime import datetime, timedelta
from pytz import timezone
from app.extensions import db

class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(120), unique=True, nullable=False)
    created_at = db.Column(db.DateTime, default=datetime.now(timezone('UTC')))
    results = db.relationship('Result', backref='user', lazy=True)
    devices = db.relationship('Device', backref='user', lazy=True)

    def __repr__(self):
        return f'<User {self.email}>'

    def is_subscription_active(self):
        expiration_date = self.created_at + timedelta(days=365)
        return datetime.now(timezone('UTC')) <= expiration_date.astimezone(timezone('UTC'))

