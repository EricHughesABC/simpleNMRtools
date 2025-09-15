# Using https://github.com/docker/awesome-compose/tree/master/flask
import sys
import os
import re
import platform
import socket
import copy
from functools import wraps
from flask import Flask, render_template, request, json, send_from_directory, url_for, jsonify, abort
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.exc import SQLAlchemyError
import time
from flask_migrate import Migrate
from flask_wtf.csrf import CSRFProtect
from sqlalchemy.dialects.mysql import JSON
from sqlalchemy import create_engine, Column, Integer, String

from datetime import datetime, timedelta
from pytz import timezone
from dotenv import load_dotenv
load_dotenv()

import numpy as np
import pandas as pd

import networkx as nx
from networkx.readwrite import json_graph


from pathlib import Path


# from flaskConfig import FlaskConfig

from html_from_assignments import NMRProblem
import jsonUtils

import nmrsolution
from globals import SVG_DIMENSIONS as svgDimensions

from simulatedAnnealing_v5 import SimulatedAnnealing2

global RUNNINGONPYTHONANYWHERE
global REGISTRATIONTIMEOUT

REGISTRATIONTIMEOUT = 100 * 60 * 60 * 24  # 100 days in seconds

# Path to the Sphinx-generated HTML documentation
DOCS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "docs/build/html")

# to test this using curl type the following in the terminal
# curl -X POST http://localhost:5000/simpleMNOVA -H "Content-Type: application/json" -d @/Users/vsmw51/Downloads/4Eric/EVB_330b_predicted/EVB_330b_predicted_assignments_mresnova.json
#
# move to the directory where the JSON file is located
# mac os
# curl -X POST "https://simplenmr.pythonanywhere.com/simpleMNOVA" -H "Content-Type: application/json" -d "@E72507_04164043_propyl_benzoate_assignments_mresnova.json" > test.html && open test.html
# 
# windows
# curl -X POST "https://simplenmr.pythonanywhere.com/simpleMNOVA" -H "Content-Type: application/json" -d "@E72507_04164043_propyl_benzoate_assignments_mresnova.json" > test.html && start test.html
#
# cross platform if python installed
# curl -X POST "https://simplenmr.pythonanywhere.com/simpleMNOVA" -H "Content-Type: application/json" -d "@E72507_04164043_propyl_benzoate_assignments_mresnova.json" > test.html && python -m webbrowser test.html

# Helper to convert numpy types to native types
def convert_numpy(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        raise TypeError(f"Type {type(obj)} not serializable")

# Determine environment (local or production)
def is_running_on_pythonanywhere():
    """Check if the application is running on PythonAnywhere"""
    # Check for PythonAnywhere-specific paths or environment variables
    return (os.path.exists('/var/www') and                # PythonAnywhere has this directory
            os.path.exists('/home/simpleNMR') or          # Check for your username's home directory
            'PYTHONANYWHERE_DOMAIN' in os.environ or      # PythonAnywhere sets this environment variable
            'PYTHONANYWHERE_SITE' in os.environ)

app = Flask(__name__)

# Database Configuration
RUNNINGONPYTHONANYWHERE = is_running_on_pythonanywhere()
if RUNNINGONPYTHONANYWHERE:
    # MySQL for PythonAnywhere (production)
    app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+pymysql://{username}:{password}@{hostname}/{database}'.format(
        username=os.environ.get('DB_USERNAME', 'simpleNMR'),
        password=os.environ.get('DB_PASSWORD', ''),
        hostname=os.environ.get('DB_HOST', 'simpleNMR.mysql.pythonanywhere-services.com'),
        database=os.environ.get('DB_NAME', 'simpleNMR$registrations')
    )
    
    # Configure connection pool for MySQL
    app.config['SQLALCHEMY_POOL_SIZE'] = int(os.environ.get('SQLALCHEMY_POOL_SIZE', 5))
    app.config['SQLALCHEMY_POOL_TIMEOUT'] = int(os.environ.get('SQLALCHEMY_POOL_TIMEOUT', 30))
    app.config['SQLALCHEMY_POOL_RECYCLE'] = int(os.environ.get('SQLALCHEMY_POOL_RECYCLE', 1800))
    app.config['SQLALCHEMY_MAX_OVERFLOW'] = int(os.environ.get('SQLALCHEMY_MAX_OVERFLOW', 2))
    
    print("Running in production mode with MySQL")
else:
    # SQLite for local development
    base_dir = os.path.abspath(os.path.dirname(__file__))
    app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///' + os.path.join(base_dir, 'registrations.db')
    print("Running in development mode with SQLite")

app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'default-secret-key')

# Print for debugging (you can remove this in production)
print("\nline Database URI:", app.config['SQLALCHEMY_DATABASE_URI'])

db = SQLAlchemy(app)
migrate = Migrate(app, db)

app.config["JSON_AS_ASCII"] = False  # contribution from Erdem

# This creates all tables based on your SQLAlchemy models
with app.app_context():
    db.create_all()
    print("Database tables created successfully!")

# Database Models
class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(120), unique=True, nullable=False)
    created_at = db.Column(db.DateTime, default=datetime.now(timezone('UTC')))
    results = db.relationship('Result', backref='user', lazy=True)

    # Relationship to devices (one-to-many)
    devices = db.relationship('Device', backref='user', lazy=True)

    def __repr__(self):
        return f'<User {self.email}>'

    def is_subscription_active(self):
        """Check if the user's subscription is still active"""
        # Example: 365-day subscription from registration date
        expiration_date = self.created_at + timedelta(days=365)
        return datetime.now(timezone('UTC')) <= expiration_date.astimezone(timezone('UTC'))

class Device(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    hostid = db.Column(db.String(100), unique=True, nullable=False)
    registered_at = db.Column(db.DateTime, default=datetime.now(timezone('UTC')))
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)
    usage_count = db.Column(db.Integer, default=0)
    last_used = db.Column(db.DateTime, default=datetime.now(timezone('UTC')))
    license_agreed = db.Column(db.Boolean, default=False, nullable=False)
    ml_consent = db.Column(db.Boolean, default=False, nullable=False)

    def __repr__(self):
        return f'<Device {self.hostid}>'

    def increment_usage(self):
        """Increment the usage count and update last_used timestamp"""
        self.usage_count += 1
        self.last_used = datetime.now(timezone('UTC'))
        db.session.commit()
    # Add any other properties or relationships here

# New Results table model
class Result(db.Model):
    __tablename__ = 'result'
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)
    smiles_string = db.Column(db.Text)
    weight = db.Column(db.Float)
    MAE = db.Column(db.Float)
    LAE = db.Column(db.Float)
    # Use JSON type for MySQL, fallback to Text for SQLite
    # json_result = db.Column(db.JSON().with_variant(db.Text, 'sqlite'), nullable=True)
    json_result = Column(JSON)
    created_at = db.Column(db.DateTime, default=datetime.now(timezone('UTC')), index=True)

    def __repr__(self):
        return f'<Result {self.id}>'

    @property
    def json_data(self):
        """Return the JSON result as a Python dictionary"""
        if self.json_result:
            return json.loads(self.json_result)
        return {}

def with_db_retries(max_retries=3, retry_delay=1):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            retries = 0

            while retries <= max_retries:
                try:
                    return func(*args, **kwargs)
                except SQLAlchemyError as e:
                    # Log the error if you have logging configured
                    if hasattr(app, 'logger'):
                        app.logger.warning(f"Database connection error in {func.__name__} (attempt {retries+1}/{max_retries+1}): {str(e)}")

                    retries += 1

                    # If this was the last retry, re-raise the exception
                    if retries > max_retries:
                        if hasattr(app, 'logger'):
                            app.logger.error(f"Failed after {max_retries+1} attempts in {func.__name__}: {str(e)}")
                        raise

                    # Wait before retrying
                    time.sleep(retry_delay)

            return None  # This should never be reached

        return wrapper

    return decorator


# Database Utility Functions
def create_database():
    """Create database tables if they don't exist"""
    with app.app_context():
        db.create_all()

create_database()

def register_device(email, hostid, license_agreed, ml_consent):
    """Register a new device or update an existing one with consent status"""

    try:
        with app.app_context():
            # Try to find existing user
            user = User.query.filter_by(email=email).first()

            # Create new user if doesn't exist
            if not user:
                user = User(email=email)
                db.session.add(user)
                db.session.flush()  # Get ID without committing

            # Check if device already exists
            device = Device.query.filter_by(hostid=hostid).first()

            if device:
                # If device exists but belongs to different user, reassign
                if device.user_id != user.id:
                    device.user_id = user.id
                device.registered_at = datetime.now(timezone('UTC'))  # Reset registration date
                # Always update consent flags when re-registering
                device.license_agreed = license_agreed
                device.ml_consent = ml_consent
            else:
                # Create new device with consent flags
                device = Device(
                    hostid=hostid,
                    user_id=user.id,
                    license_agreed=license_agreed,
                    ml_consent=ml_consent,
                    last_used=datetime.now(timezone('UTC'))  # Set last_used to now
                )
                db.session.add(device)

            db.session.commit()
            return True
    except Exception as e:
        print(f"Error during registration: {e}")
        return False
    
@with_db_retries(max_retries=3, retry_delay=1)
def is_device_registered(hostid):
    """Check if a device is registered and if the registration is still valid"""
    with app.app_context():
        device = Device.query.filter_by(hostid=hostid).first()

        print("is_device_registered")
        if not device:
            print("\tdevice is False")
            return False

        # Check if user's subscription is still active
        return device.user.is_subscription_active()

def record_usage(hostid):
    """Record usage of the software by this device"""
    with app.app_context():
        device = Device.query.filter_by(hostid=hostid).first()
        if device:
            device.increment_usage()
            return True
        return False

def get_device_stats(hostid):
    """Get statistics for a device"""
    with app.app_context():
        device = Device.query.filter_by(hostid=hostid).first()
        if not device:
            return None

        user = device.user
        days_remaining = 0

        if user:
            expiration_date = user.created_at + timedelta(days=365)
            days_remaining = (expiration_date.astimezone(timezone('UTC')) - datetime.now(timezone('UTC'))).days

        return {
            'email': user.email,
            'registered_at': device.registered_at,
            'usage_count': device.usage_count,
            'last_used': device.last_used,
            'days_remaining': max(0, days_remaining),
            'ml_consent': device.ml_consent,
        }

def get_user_id_hostid(hostid):
    """Get user ID for a given host ID"""
    with app.app_context():
        device = Device.query.filter_by(hostid=hostid).first()
        if device:
            return device.user_id
    return None

def has_device_expired(hostid):
    """Check if the device has expired"""
    with app.app_context():
        device = Device.query.filter_by(hostid=hostid).first()
        if device:
            user = device.user
            registered_at_UTC = device.registered_at.astimezone(timezone('UTC'))
            #  print all the attributes of the user
            print(f"User attributes: {user.__dict__}")
            if user:
                expiration_date_UTC = registered_at_UTC
                # convert expriation date to UTC
                now_UTC = (datetime.now()).astimezone(timezone('UTC'))

                # if (now_UTC-registered_at_UTC).seconds > 100*60*60*24:
                if (now_UTC-registered_at_UTC).seconds > REGISTRATIONTIMEOUT:
                    print(f"Device {hostid} has expired")
                    return True
                else:
                    print(f"Device {hostid} has not expired")

    return False

def update_ml_consent_for_all_user_devices(hostname, consent_value):
    # First, find the device with the provided hostname
    device = Device.query.filter_by(hostid=hostname).first()

    if device:
        # Get the user_id from this device
        user_id = device.user_id

        # Update ml_consent for all devices belonging to this user
        Device.query.filter_by(user_id=user_id).update({Device.ml_consent: consent_value})

        # Commit the changes to the database
        db.session.commit()

        # Return the number of devices updated
        return Device.query.filter_by(user_id=user_id).count()
    else:
        # Return 0 if no device was found with the provided hostname
        return 0

@app.shell_context_processor
def make_shell_context():
    """Provide objects for interactive Flask shell sessions.

    This function returns a dictionary of commonly used objects and models for use in the Flask shell.
    It allows for easier access to the database and model classes during development.

    Returns:
        dict: A dictionary mapping names to objects for the Flask shell context.
    """
    return {'sa': SQLAlchemy, 'db': db, 'User': User, 'Device': Device, 'Result': Result}

@app.route('/download')
def download_file():
    """Serve the MNOVAscripts.zip file for download.

    This endpoint allows users to download the MNOVAscripts.zip file from the downloads directory.
    Returns a 404 error if the file does not exist.

    Returns:
        Response: The file as an attachment if found, otherwise a 404 error.
    """
    # Define the path to the downloads directory using pathlib
    downloads_dir = Path(app.root_path) / 'downloads'
    filename = 'MNOVAscripts.zip'
    file_path = downloads_dir / filename

    # Debug print
    print(f"Looking for file at: {file_path}")

    # Check if the file exists
    if not file_path.exists():
        print("File not found!")
        abort(404)

    # Serve the file
    return send_from_directory(str(downloads_dir), filename, as_attachment=True)


@app.route("/documentation/")
def index():
    """Serve the main documentation index page.

    This endpoint returns the main HTML documentation page generated by Sphinx.
    
    Returns:
        Response: The index.html documentation file.
    """  
    return send_from_directory(DOCS_DIR, "index.html")


@app.route("/documentation/<path:filename>")
def docs(filename):
    """Serve a specific documentation file from the Sphinx-generated docs.

    This endpoint returns the requested documentation file from the documentation directory.

    Args:
        filename (str): The path to the documentation file to serve.

    Returns:
        Response: The requested documentation file.
    """
    return send_from_directory(DOCS_DIR, filename)


@app.route("/", methods=["GET"])
def display_front_page():
    """Render and serve the main front page of the application.

    This endpoint displays the mainPage.html template as the application's front page.

    Returns:
        Response: The rendered mainPage.html template.
    """
    return render_template("mainPage.html")


@app.route("/simpleMNOVAfinalHTML", methods=["POST"])
def simpleMNOVAfinalHTML():
    """Process and render the final HTML for a simpleMNOVA molecule visualization.

    This endpoint receives JSON data, processes molecular graph information, and renders an HTML template for visualization.
    It also handles device registration, machine learning consent, and optionally stores results in the database.

    Returns:
        Response: The rendered HTML template or a JSON response indicating registration status.
    """
    if request.method != "POST":
        return "Only POST requests are accepted", 400

    # Get JSON data from the request body (for curl POST requests)
    try:
        json_data = request.get_json()
        if json_data is None:
            print("No JSON data received")
            return "No JSON data received", 400

    except json.JSONDecodeError as e:
        print(f"Invalid JSON: {e}")
        return f"Invalid JSON: {e}", 400

    id_now_mapping_dict = {}

    # Create a dictionary for quick lookup of catoms_orig by atomNumber
    catoms_orig_lookup = {
        node["atomNumber"]: node["id"] for node in json_data["catoms_orig"]
    }

    # Iterate through nodes_now and create the mapping
    for node in json_data["nodes_now"]:
        atom_number = node["atomNumber"]
        if atom_number in catoms_orig_lookup:
            id_now_mapping_dict[node["id"]] = catoms_orig_lookup[atom_number]
        else:
            print(f"atomNumber {atom_number} not found in catoms_orig_lookup")

    catoms_orig_dict = {}

    for node in json_data["catoms_orig"]:
        catoms_orig_dict[node["id"]] = node

    links_moved = []

    for link in json_data["links"]:
        source_id = id_now_mapping_dict[link["source"]]
        target_id = id_now_mapping_dict[link["target"]]
        new_link = link.copy()
        new_link["source"] = source_id
        new_link["target"] = target_id
        links_moved.append(new_link)

    # update nodes
    nodes_moved_new = []

    for node in json_data["nodes_now"]:
        new_node = node.copy()
        new_node["id"] = id_now_mapping_dict[node["id"]]
        nodes_moved_new.append(new_node)

    rtn_html = render_template(
        "d3molplotmnova_template.html",
        graph_edges=links_moved,
        graph_nodes=nodes_moved_new,
        orig_nodes=nodes_moved_new,
        molgraph=json_data.get("molgraph", "'dummy'"),
        shortest_paths=json_data.get("shortest_paths", "'dummy'"),
        svg_container=json_data["svg"],
        title=json_data["title"],
        smilesString=json_data["smilesString"],
        molFile=json_data["molfile"],
        workingDirectory=json_data["workingDirectory"],
        workingFilename=json_data["workingFilename"],
        dataFrom=json_data["dataFrom"],
        catoms=json_data["catoms_orig"],
        oldjsondata=json_data["oldjsondata"],
        best_results=json_data["best_results"],
    )

    rtn_html = rtn_html.replace("True", "true")
    rtn_html = rtn_html.replace("False", "false")


    hostname = jsonUtils.extract_hostname(json_data["oldjsondata"]) # returns a string or None

        # Check if device is registered and registration is valid
    if not is_device_registered(hostname):
        return jsonify({
            'status': 'unregistered',
            'registration_url': url_for('registration_page', _external=True) + f"?hostid={hostname}"
        })

    if has_device_expired(hostname):
        email = get_device_stats(hostname)["email"]
        return jsonify({
            'status': 'registration_expired',
            'registration_url': url_for('registration_page', _external=True) + f"?hostid={hostname}" + f"&email={email}",
        })


    try:
        machine_learning_opt_in = json_data["oldjsondata"]["ml_consent"]["data"]["0"]
    except KeyError:
        # If the key doesn't exist, set machine_learning_opt_in to False
        machine_learning_opt_in = False

    updated_count = update_ml_consent_for_all_user_devices(hostname,  machine_learning_opt_in)

    if machine_learning_opt_in:
        # create a jinja2 template from the HTML

        jinja_template = {}
        jinja_template["svg_container"] = json_data["svg"]
        jinja_template["graph_edges"] = links_moved
        jinja_template["graph_nodes"] = nodes_moved_new
        jinja_template["orig_nodes"] = nodes_moved_new
        jinja_template["molgraph"] = json_data.get("molgraph", "'dummy'")
        jinja_template["shortest_paths"] = json_data.get("shortest_paths", "'dummy'")
        jinja_template["title"] = json_data["title"]
        jinja_template["smilesString"] = json_data["smilesString"]
        jinja_template["molFile"] = json_data["molfile"]
        jinja_template["workingDirectory"] = json_data["workingDirectory"]
        jinja_template["workingFilename"] = json_data["workingFilename"]
        jinja_template["dataFrom"] = json_data["dataFrom"]
        jinja_template["catoms"] = json_data["catoms_orig"]
        jinja_template["oldjsondata"] = json_data["oldjsondata"]
        jinja_template["best_results"] = json_data["best_results"]

        # check if the user has opted in for machine learning
        user_id = get_user_id_hostid(hostname)

        # Get the complete JSON result
        json_result = json.loads(json.dumps(jinja_template, default=convert_numpy))
        if not json_result:
            return rtn_html

        # Convert to string if it's already a dict
        if isinstance(json_result, dict):

            json_result_str = json.dumps(json_result, default=convert_numpy)
        else:
            json_result_str = json_result
            # Validate it's proper JSON
            try:
                json_data = json.loads(json_result_str)
            except json.JSONDecodeError:
                return rtn_html

        # check if the user id is valid
        if user_id is not None:

            # save the results to the database
            with app.app_context():
                if RUNNINGONPYTHONANYWHERE:
                    # MySQL for PythonAnywhere (production)
                    print("\n##########################################")
                    print("Running in production mode with MySQL")
                    print("##########################################\n")
                    new_result = Result(
                        user_id=user_id,
                        smiles_string=json_data["smilesString"],
                        weight=json_data["best_results"]["best_weight"],
                        MAE=json_data["best_results"]["best_mae"],
                        LAE=json_data["best_results"]["best_lae"],
                        json_result=json_result,
                        
                        # Use JSON type for MySQL, fallback to Text for SQLite
                    )
                else:

                    new_result = Result(
                        user_id=user_id,  # Replace with actual user ID
                        smiles_string=json_data["smilesString"],
                        weight=json_data["best_results"]["best_weight"],
                        MAE=json_data["best_results"]["best_mae"],
                        LAE=json_data["best_results"]["best_lae"],
                        json_result=json_result,
                        
                        # Use JSON type for MySQL, fallback to Text for SQLite
                        # json_result=json.dumps(json_result, default=convert_numpy),

                    )
                db.session.add(new_result)
                db.session.commit()

    return rtn_html


@app.route("/check_machine_learning", methods=["POST"])
def check_machine_learning():
    """Check the machine learning consent status for a registered device.

    This endpoint verifies device registration and expiration, then returns the machine learning consent status.
    It responds with registration status and a registration URL if the device is unregistered or expired.

    Returns:
        Response: A JSON object with registration status and machine learning consent information.
    """

    print("check_machine_learning")
    json_data = request.get_json()

    print("json_data\n", json_data)

    print("request\n", request)
    print("request.args\n", request.args)
    hostname = json_data.get("hostname")

    # Check if device is registered and registration is valid
    if not is_device_registered(hostname):
        print(f"\nDevice {hostname} is not registered\n")
        return jsonify({
            'status': 'unregistered',
            'registration_url': url_for('registration_page', _external=True) + f"?hostid={hostname}"
        })

    if has_device_expired(hostname):
        print(f"\nDevice {hostname} has expired\n")
        email = get_device_stats(hostname)["email"]
        return jsonify({
            'status': 'registration_expired',
            'registration_url': url_for('registration_page', _external=True) + f"?hostid={hostname}" + f"&email={email}",
        })

    # get machine learning consent status
    ml_consent = get_device_stats(hostname)["ml_consent"]
    return jsonify({"status": "registered", "ml_consent": ml_consent})


@app.route("/register", methods=["POST"])
def register_host():
    """Register a device and user with license and machine learning consent.

    This endpoint processes registration form data, validates required fields, and registers the device and user in the database.
    It returns a success or error page based on the registration outcome.

    Returns:
        Response: The rendered registration success or error template.
    """
    hostid = request.form.get("hostid")
    email = request.form.get("email")
    license_agreed = request.form.get("agree-license") == "on"
    ml_consent = request.form.get("ml-consent") == "on"

    print(f"Register host ID: {hostid}, email: {email}, license_agreed: {license_agreed}, ml_consent: {ml_consent}")

    if not hostid or not email:
        return "Both Host ID and email are required", 400

    if not license_agreed:
        return "You must agree to the license agreement", 400

    # Register the device in the database with consent status
    success = register_device(email, hostid, license_agreed, ml_consent)

    if success:
        return render_template("registration_success.html",
                              hostid=hostid,
                              email=email,
                              ml_consent=ml_consent)
    else:
        return render_template("registration_error.html")


@app.route("/registration", methods=["GET"])
def registration_page():
    """Render the registration and agreement page for a device.

    This endpoint serves the registration page, requiring a hostid query parameter to identify the device.
    It returns an error if the hostid is missing.

    Returns:
        Response: The rendered registration and agreement template or an error message.
    """

    hostid = request.args.get("hostid")

    # Check if hostid is provided
    if not hostid:
        return "No hostid provided. Please include a hostid query parameter.", 400

    # Render the registration page with the hostid
    return render_template("registration_and_agreement_page.html", hostid=hostid)


@app.route("/simpleMNOVA", methods=["POST"])
def simpleMNOVA_display_molecule():
    """Process and render the HTML visualization for a simpleMNOVA molecule.

    This endpoint receives JSON data, processes NMR molecular graph information, and renders an HTML template for visualization.
    It manages device registration, machine learning consent, and optionally stores results in the database.

    Returns:
        Response: The rendered HTML template or a JSON response indicating registration status or errors.
    """
    if request.method != "POST":
        return "Only POST requests are accepted", 400

    # Get JSON data from the request body (for curl POST requests)
    try:
        json_data = request.get_json()
        print("json_data.keys()\n", json_data.keys())
        if json_data is None:
            print("No JSON data received")
            return "No JSON data received", 400

    except json.JSONDecodeError as e:
        print(f"Invalid JSON: {e}")
        return f"Invalid JSON: {e}", 400

    # check if hostname is already known
    # print out the hostname

    hostname = jsonUtils.extract_hostname(json_data) # returns a string or None

    # Check if device is registered and registration is valid
    if not is_device_registered(hostname):
        return jsonify({
            'status': 'unregistered',
            'registration_url': url_for('registration_page', _external=True) + f"?hostid={hostname}"
        })

    if has_device_expired(hostname):
        email = get_device_stats(hostname)["email"]
        return jsonify({
            'status': 'registration_expired',
            'registration_url': url_for('registration_page', _external=True) + f"?hostid={hostname}" + f"&email={email}",
        })

    # Record this usage
    record_usage(hostname)

    # check if the user has agreed to the machine learning consent
    try:
        machine_learning_opt_in = json_data["ml_consent"]["data"]["0"]
    except KeyError:
        # If the key doesn't exist, set machine_learning_opt_in to False
        machine_learning_opt_in = False

    updated_count = update_ml_consent_for_all_user_devices(hostname,  machine_learning_opt_in)

    problemdata_json = NMRProblem.from_mnova_dict(json_data)

    print("problemdata_json\n", problemdata_json.dataframes.keys())


    # decide whether we are doing prediction or assignments
    if problemdata_json.is_prediction():

        solution = nmrsolution.NMRsolution(problemdata_json)


        print("solution.hsqc_df.shape ", solution.hsqc_df.shape[0])
        print(solution.hsqc_df.columns )


        if solution.nmrsolution_failed:
            return solution.nmrsolution_error_message, solution.nmrsolution_error_code

        if solution.expected_molecule.nmrshiftdb_failed:
            print("solution.expected_molecule.nmrshiftdb_failed_message\n")
            print(solution.expected_molecule.nmrshiftdb_failed_message)
            # return solution.expected_molecule.nmrshiftdb_failed_message, solution.solution_error_code
            return solution.solution_error_message, solution.solution_error_code

        ok, msg = solution.init_class_from_json()

        if not ok:
            return msg, 400

        # attempt to assign integrals to the HSQC dataframe

        # solution.assign_CH3_CH2_CH1_in_HSQC_using_Assignments()

        # solution.hsqc["numProtons"] = -1

        rtn_msg, rtn_num = solution.assign_CH3_CH2_CH1_overall()

        print("solution.hsqc.shape ", solution.hsqc.shape[0])
        print(solution.hsqc.columns )

        print("solution.hsqc.CH2 ", solution.hsqc[solution.hsqc.CH2].shape[0])
        print("solution.hsqc.CH3 ", solution.hsqc[solution.hsqc.CH3].shape[0])
        print("solution.hsqc.CH1 ", solution.hsqc[solution.hsqc.CH1].shape[0])
        print("solution.hsqc.CH3CH1 ", solution.hsqc[solution.hsqc.CH3CH1].shape[0])

        # print the mol_props_df

        df = solution.expected_molecule.molprops_df
        print("molprops_df.shape ", df.shape[0])

        print("CH2 count: ", df[df["CH2"]].shape[0])
        print("CH3 count: ", df[df["CH3"]].shape[0])
        print("CH1 count: ", df[df["CH1"]].shape[0])
        print("CH3CH1 count: ", df[df["CH3CH1"]].shape[0])

        # print information on the c13 dataframe
        print("solution.c13.shape ", solution.c13.shape[0])
        print(solution.c13.columns )    
        print("solution.c13.CH2 ", solution.c13[solution.c13.CH2].shape[0])
        print("solution.c13.CH3 ", solution.c13[solution.c13.CH3].shape[0])
        print("solution.c13.CH1 ", solution.c13[solution.c13.CH1].shape[0])
        print("solution.c13.CH3CH1 ", solution.c13[solution.c13.CH3CH1].shape[0])

        print("\nrtn_msg\n", rtn_msg)

        if rtn_msg != "ok":
            return rtn_msg, rtn_num

        solution.transfer_hsqc_info_to_c13()
        solution.transfer_hsqc_info_to_h1()

        print("solution.c13.shape ", solution.c13.shape[0])
        print(solution.c13.columns )    
        print("solution.c13.CH2 ", solution.c13[solution.c13.CH2].shape[0])
        print("solution.c13.CH3 ", solution.c13[solution.c13.CH3].shape[0])
        print("solution.c13.CH1 ", solution.c13[solution.c13.CH1].shape[0])
        print("solution.c13.CH3CH1 ", solution.c13[solution.c13.CH3CH1].shape[0])

        rtn_msg, rtn_num = solution.initialise_prior_to_carbon_assignment()

        if rtn_msg != "ok":
            return rtn_msg, rtn_num

        rtn_msg, rtn_num = solution.attempt_assignment_CH3_CH2_CH1_to_C13_table()

        if rtn_msg != "ok":
            return rtn_msg, rtn_num

        solution.update_assignments_expt_dataframes()

        G2 = nmrsolution.create_network_graph(solution.c13, solution.h1)
        G2 = nmrsolution.add_all_cosy_edges_to_graph(
            G2, solution.cosy, solution.hsqc_clipcosy, solution.h1
        )
        G2 = nmrsolution.add_all_hmbc_edges_to_graph(
            G2, solution.hmbc, solution.h1, solution.c13
        )

        jsonGraphData = json_graph.node_link_data(G2)
        jsonGraphData["moved_nodes"] = jsonGraphData["nodes"]

        # create a network graph of the expected molecule

        solution.initiate_molgraph( json_data, G2)

        #  convert the molgraph to a json object
        jsonGraphData_mol = json_graph.node_link_data(solution.molgraph)

        # calculate the shortest paths between all pairs of nodes in the molgraph
        shortest_paths = dict(nx.all_pairs_dijkstra_path_length(solution.molgraph))

        svgWidth = svgDimensions.svg_width
        svgHeight = svgDimensions.svg_height
        molWidth = svgDimensions.mol_width
        molHeight = svgDimensions.mol_height

        catoms_df = solution.expected_molecule.molprops_df[
            [
                "ppm",
                "atom_idx",
                "atomNumber",
                "numProtons",
                "x",
                "y",
            ]
        ].copy()

        # rename ppm to ppm_calculated
        catoms_df.rename(columns={"ppm": "ppm_calculated"}, inplace=True)
        # renam atom_idx to id
        catoms_df.rename(columns={"atom_idx": "id"}, inplace=True)

        # add columns ppm, jCouplingVals, jCouplingClass, visible, H1_ppm
        catoms_df["ppm"] = -1000000.0
        catoms_df["jCouplingVals"] = ""
        catoms_df["jCouplingClass"] = ""
        catoms_df["visible"] = "false"
        catoms_df["symbol"] = "C"
        catoms_df["iupacLabel"] = ""

        # catoms_str = json.dumps(catoms_df.to_dict(orient="records"), indent=4)
        json_data_str = json.dumps(json_data, indent=4)

        if problemdata_json.prediction_from_nmrshiftdb2():
            dataFrom = "nmrshiftdb2"
        else:
            dataFrom = "mnova"

        simAnneal = SimulatedAnnealing2.from_params(
            copy.deepcopy(jsonGraphData["nodes"]),
            copy.deepcopy(jsonGraphData["links"]),
            json_data["molfile"]["data"]["0"],
            json_data,
        )

        # setup the run using standaard parameters hardcoded in the class
        simAnneal.setup_run()

        if json_data["simulatedAnnealing"]["data"]["0"] and (simAnneal.predicted_weight > 0):  # True

            # simAnneal.setup_run(randomize_mapping=True)
            simAnneal.run_optimization(50)

            jsonGraphData, best_results = simAnneal.process_results(catoms_df, jsonGraphData)

        else:
            best_results = simAnneal.process_results_SA_skipped()


        #  start  to produce the html display

        id_atomNumber = []

        for node in jsonGraphData["moved_nodes"]:
            id_atomNumber.append([node["id"], node["atomNumber"]])

        # sort id_atomNumber
        id_atomNumber.sort(key=lambda x: x[0])
        id_x = [x[1] for x in id_atomNumber]

        # copy over the optimized nodes to the catoms_df


        for idx, row in catoms_df.iterrows():
            id = row["id"]
            for node in jsonGraphData["moved_nodes"]:
                if node["id"] == id:
                    catoms_df.loc[idx, "ppm"] = node["ppm"]
                    catoms_df.loc[idx, "x"] = node["x"]
                    catoms_df.loc[idx, "y"] = node["y"]
                    catoms_df.loc[idx, "ppm_calculated"] = node["ppm_calculated"]
                    catoms_df.loc[idx, "atomNumber"] = node["atomNumber"]


        catoms_str = json.dumps(catoms_df.to_dict(orient="records"), indent=4)
        # catoms_str = json.dumps(jsonGraphData["nodes"], indent=4)

        jinja_template = {
            "svg_container": solution.expected_molecule.svg_str,
            "graph_edges": jsonGraphData["links"],
            "graph_nodes": jsonGraphData["moved_nodes"],
            "orig_nodes": jsonGraphData["nodes"],
            "molgraph": jsonGraphData_mol,
            "shortest_paths": json.dumps(shortest_paths),
            "catoms": catoms_str,
            "title": "dummy_title",
            "smilesString": solution.smilesstr,
            "molFile": solution.molstr,
            "workingDirectory": problemdata_json.dataframes["workingDirectory"]
            .loc[0, "workingDirectory"]
            .strip(),
            "workingFilename": problemdata_json.dataframes["workingFilename"]
            .loc[0, "workingFilename"]
            .strip(),
            "dataFrom": dataFrom,
            "oldjsondata": json_data_str,
            "best_results": best_results,
        }
        msg = "ok"

    else:
        problemdata_json.prepare_network_graph()

        json_data_str = json.dumps(json_data, indent=4)

        jinja_template = problemdata_json.jinjadata
        jinja_template["smilesString"] = problemdata_json.dataframes["smiles"].loc[
            0, "smiles"
        ]
        jinja_template["molFile"] = problemdata_json.dataframes["molfile"].loc[
            0, "molfile"
        ]
        jinja_template["workingDirectory"] = (
            problemdata_json.dataframes["workingDirectory"]
            .loc[0, "workingDirectory"]
            .strip()
        )
        jinja_template["workingFilename"] = (
            problemdata_json.dataframes["workingFilename"]
            .loc[0, "workingFilename"]
            .strip()
        )
        jinja_template["oldjsondata"] = json_data_str

        if problemdata_json.prediction_from_nmrshiftdb2():
            dataFrom = "nmrshiftdb2"
        else:
            dataFrom = "mnova"

        jinja_template["dataFrom"] = dataFrom

        jinja_template["best_results"] = {
            "best_weight": 0,
            "best_mae": 0.0,
            "best_lae": 0.0,
        }

        # add molgaph to jinja_template

        msg = "ok"

    # self.problemdata_json.prediction_from_nmrshiftdb2()
    if msg == "ok":

        # check if machine learning opt in is true and if so save the data to the database
        # check if the user has opted in for machine learning

        rtn_html = render_template(
            "d3molplotmnova_template.html",
            graph_edges=jinja_template["graph_edges"],
            graph_nodes=jinja_template["graph_nodes"],
            orig_nodes=jinja_template["graph_nodes"],
            molgraph=jinja_template.get("molgraph", "'dummy'"),
            shortest_paths=jinja_template.get("shortest_paths", "'dummy'"),
            svg_container=jinja_template["svg_container"],
            title=jinja_template["title"],
            smilesString=jinja_template["smilesString"],
            molFile=jinja_template["molFile"],
            workingDirectory=jinja_template["workingDirectory"],
            workingFilename=jinja_template["workingFilename"],
            dataFrom=jinja_template["dataFrom"],
            catoms=jinja_template["catoms"],
            oldjsondata=jinja_template["oldjsondata"],
            best_results=jinja_template["best_results"],
        )

        rtn_html = rtn_html.replace("True", "true")
        rtn_html = rtn_html.replace("False", "false")

        # replace np.float64
        float_pattern = r'np\.float64\(([\d\.]+)\)'
        int_pattern = r'np\.int64\(([\d]+)\)'
        rtn_html =  re.sub(float_pattern, r'\1', rtn_html)
        rtn_html =  re.sub(int_pattern, r'\1', rtn_html)

        # add redults to database if the user has opted in for machine learning

        if machine_learning_opt_in:
            # check if the user has opted in for machine learning

            user_id = get_user_id_hostid(hostname)

            # Get the complete JSON result
            json_result = json.loads(json.dumps(jinja_template, default=convert_numpy))
            if not json_result:
                return rtn_html

            # Convert to string if it's already a dict
            if isinstance(json_result, dict):

                json_result_str = json.dumps(json_result, default=convert_numpy)
            else:
                json_result_str = json_result
                # Validate it's proper JSON
                try:
                    json_data = json.loads(json_result_str)
                except json.JSONDecodeError:
                    return rtn_html

            # # check if the user id is valid
            if user_id is not None:

            # save the results to the database

                print("json_data.keys()\n", json_data.keys())
                with app.app_context():
                    if RUNNINGONPYTHONANYWHERE:
                        # MySQL for PythonAnywhere (production)
                        print("\n##########################################")
                        print("Running in production mode with MySQL 2")
                        print("##########################################\n")
                        new_result = Result(
                            user_id=user_id,
                            smiles_string=solution.smilesstr,
                            weight=best_results["best_weight"],
                            MAE=best_results["best_mae"],
                            LAE=best_results["best_lae"],
                            json_result=json_result,
                            
                            # Use JSON type for MySQL, fallback to Text for SQLite
                        )
                    else:

                        new_result = Result(
                            user_id=user_id,  # Replace with actual user ID
                            smiles_string=solution.smilesstr,
                            weight=best_results["best_weight"],
                            MAE=best_results["best_mae"],
                            LAE=best_results["best_lae"],
                            
                            # Use JSON type for MySQL, fallback to Text for SQLite
                            # json_result=json.dumps(json_result, default=convert_numpy),
                            json_result=json_result,

                        )
                    db.session.add(new_result)
                    db.session.commit()
        return rtn_html
    else:
        print("Error in processing the data\n", msg)
        return msg, 400


if __name__ == "__main__":
    # reload(sys)
    # sys.setdefaultencoding("utf-8")
    create_database()  # Create tables if they don't exist
    app.run(debug=True)
