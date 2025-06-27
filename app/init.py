from flask import Flask
from app.extensions import db, migrate
from app.config import Config

def create_app(config_class=Config):
    app = Flask(__name__)
    app.config.from_object(config_class)
    app.config["JSON_AS_ASCII"] = False

    # Initialize extensions
    db.init_app(app)
    migrate.init_app(app, db)

    # Import models after db initialization
    from app.models import User, Result, Device

    # Create tables
    with app.app_context():
        db.create_all()
        print("Database tables created successfully!")

    return app

