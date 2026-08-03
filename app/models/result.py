"""Result model — stores NMR assignment results per user."""

import json
from datetime import datetime

from pytz import timezone

from app.extensions import db


class Result(db.Model):
    __tablename__ = "result"

    id            = db.Column(db.Integer, primary_key=True)
    user_id       = db.Column(db.Integer, db.ForeignKey("user.id"), nullable=False)
    smiles_string = db.Column(db.Text)
    weight        = db.Column(db.Float)
    MAE           = db.Column(db.Float)
    LAE           = db.Column(db.Float)
    # db.JSON() with SQLite fallback makes the column portable across MySQL,
    # Postgres, and SQLite without any additional configuration.
    json_result   = db.Column(
        db.JSON().with_variant(db.Text, "sqlite"), nullable=True
    )
    created_at    = db.Column(
        db.DateTime, default=datetime.now(timezone("UTC")), index=True
    )

    def __repr__(self) -> str:
        return f"<Result {self.id}>"

    @property
    def json_data(self) -> dict:
        """Return the stored JSON result as a Python dictionary."""
        if self.json_result:
            return json.loads(self.json_result)
        return {}
