"""Device model — represents a registered host machine."""

from datetime import datetime

from pytz import timezone

from app.extensions import db


class Device(db.Model):
    id             = db.Column(db.Integer, primary_key=True)
    hostid         = db.Column(db.String(100), unique=True, nullable=False)
    registered_at  = db.Column(db.DateTime, default=datetime.now(timezone("UTC")))
    user_id        = db.Column(db.Integer, db.ForeignKey("user.id"), nullable=False)
    usage_count    = db.Column(db.Integer, default=0)
    last_used      = db.Column(db.DateTime, default=datetime.now(timezone("UTC")))
    license_agreed = db.Column(db.Boolean, default=False, nullable=False)
    ml_consent     = db.Column(db.Boolean, default=False, nullable=False)

    def __repr__(self) -> str:
        return f"<Device {self.hostid}>"

    def increment_usage(self) -> None:
        """Increment usage count and update last_used timestamp."""
        self.usage_count += 1
        self.last_used = datetime.now(timezone("UTC"))
        db.session.commit()
