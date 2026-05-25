"""
Flask routes for simpleNMR.

All route handlers and database helper functions extracted from
simpleNMRtest_app.py as a Flask Blueprint.
"""

import os
import re
import copy
import json
import time
from datetime import datetime, timedelta
from functools import wraps
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph
from sqlalchemy.exc import SQLAlchemyError
from pytz import timezone
from loguru import logger

from flask import (
    Blueprint,
    abort,
    current_app,
    jsonify,
    render_template,
    request,
    send_from_directory,
    url_for,
)

from app.extensions import db
from app.models import User, Device, Result

from core.html_from_assignments import NMRProblem
import utils.json_utils as jsonUtils
import core.expectedmolecule as expectedmolecule
import core.nmrsolution as nmrsolution
from config.globals import SVG_DIMENSIONS as svgDimensions
from core.simulated_annealing import SimulatedAnnealing2

bp = Blueprint("main", __name__)

# ── Module-level constants ────────────────────────────────────────────────────

REGISTRATIONTIMEOUT = 1000 * 60 * 60 * 24  # 100 days in seconds

# Path to Sphinx-generated HTML docs — project root is one level above app/
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
DOCS_DIR = str(_PROJECT_ROOT / "docs" / "build" / "html")

# ── Utilities ────────────────────────────────────────────────────────────────
def _node_link_data(G) -> dict:
    """Version-safe wrapper for json_graph.node_link_data.

    NX < 3.3  — no edges kwarg; default produces {"links": [...]}
    NX 3.3–3.5 — edges kwarg added; explicit "links" silences FutureWarning
    NX >= 3.6  — default changed to "edges"; explicit "links" preserves key

    Always produces {"links": [...]} so downstream code is stable across
    all versions. Revisit when upgrading server beyond NX 3.6.
    """
    try:
        return json_graph.node_link_data(G, edges="links")
    except TypeError:  # NX < 3.3 — edges kwarg not yet supported
        return json_graph.node_link_data(G)
from networkx.readwrite import json_graph

def convert_numpy(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        raise TypeError(f"Type {type(obj)} not serializable")


# MODIFIED: replaced is_running_on_pythonanywhere() with get_database_uri().

# ── Database retry decorator ─────────────────────────────────────────────────
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
                    if hasattr(app, "logger"):
                        logger.warning(
                            f"Database connection error in {func.__name__} (attempt {retries+1}/{max_retries+1}): {str(e)}"
                        )

                    retries += 1

                    # If this was the last retry, re-raise the exception
                    if retries > max_retries:
                        if hasattr(app, "logger"):
                            logger.error(
                                f"Failed after {max_retries+1} attempts in {func.__name__}: {str(e)}"
                            )
                        raise

                    # Wait before retrying
                    time.sleep(retry_delay)

            return None  # This should never be reached

        return wrapper

    return decorator


# ── Database helper functions ────────────────────────────────────────────────
def create_database():
    """Create database tables if they don't exist"""
    with current_app.app_context():
        db.create_all()


# create_database() # removed as it is called in __main__


def register_device(email, hostid, license_agreed, ml_consent):
    """Register a new device or update an existing one with consent status"""

    try:
        with current_app.app_context():
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
                device.registered_at = datetime.now(
                    timezone("UTC")
                )  # Reset registration date
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
                    last_used=datetime.now(timezone("UTC")),  # Set last_used to now
                )
                db.session.add(device)

            db.session.commit()
            return True
    except Exception as e:
        logger.error(f"Registration error: {e}")
        return False


@with_db_retries(max_retries=3, retry_delay=1)
def is_device_registered(hostid):
    """Check if a device is registered and if the registration is still valid"""
    with current_app.app_context():
        device = Device.query.filter_by(hostid=hostid).first()

        logger.debug("is_device_registered called")
        if not device:
            logger.debug("Device not found")
            return False
        
        else:
            logger.debug("Device found")
            return True

        # Check if user's subscription is still active
        # return device.user.is_subscription_active()


def record_usage(hostid):
    """Record usage of the software by this device"""
    with current_app.app_context():
        device = Device.query.filter_by(hostid=hostid).first()
        if device:
            device.increment_usage()
            return True
        return False


def get_device_stats(hostid):
    """Get statistics for a device"""
    with current_app.app_context():
        device = Device.query.filter_by(hostid=hostid).first()
        if not device:
            return None

        user = device.user
        days_remaining = 0

        if user:
            expiration_date = user.created_at + timedelta(days=365)
            days_remaining = (
                expiration_date.astimezone(timezone("UTC"))
                - datetime.now(timezone("UTC"))
            ).days

        return {
            "email": user.email,
            "registered_at": device.registered_at,
            "usage_count": device.usage_count,
            "last_used": device.last_used,
            "days_remaining": max(0, days_remaining),
            "ml_consent": device.ml_consent,
        }


def get_user_id_hostid(hostid):
    """Get user ID for a given host ID"""
    with current_app.app_context():
        device = Device.query.filter_by(hostid=hostid).first()
        if device:
            return device.user_id
    return None


def has_device_expired(hostid):
    """Check if the device has expired"""
    with current_app.app_context():
        device = Device.query.filter_by(hostid=hostid).first()
        if device:
            user = device.user
            registered_at_UTC = device.registered_at.astimezone(timezone("UTC"))
            #  print all the attributes of the user
            logger.debug(f"User attributes: {user.__dict__}")
            if user:
                expiration_date_UTC = registered_at_UTC
                # convert expriation date to UTC
                now_UTC = (datetime.now()).astimezone(timezone("UTC"))

                # if (now_UTC-registered_at_UTC).seconds > 100*60*60*24:
                if (now_UTC - registered_at_UTC).seconds > REGISTRATIONTIMEOUT:
                    logger.info(f"Device {hostid} has expired")
                    return True
                else:
                    logger.debug(f"Device {hostid} is active")

    return False


def update_ml_consent_for_all_user_devices(hostname, consent_value):
    # First, find the device with the provided hostname
    device = Device.query.filter_by(hostid=hostname).first()

    if device:
        # Get the user_id from this device
        user_id = device.user_id

        # Update ml_consent for all devices belonging to this user
        Device.query.filter_by(user_id=user_id).update(
            {Device.ml_consent: consent_value}
        )

        # Commit the changes to the database
        db.session.commit()

        # Return the number of devices updated
        return Device.query.filter_by(user_id=user_id).count()
    else:
        # Return 0 if no device was found with the provided hostname
        return 0


# ── Routes ───────────────────────────────────────────────────────────────────
@bp.route("/download")
def download_file():
    """Serve the MNOVAscripts.zip file for download.

    This endpoint allows users to download the MNOVAscripts.zip file from the downloads directory.
    Returns a 404 error if the file does not exist.

    Returns:
        Response: The file as an attachment if found, otherwise a 404 error.
    """
    # Define the path to the downloads directory using pathlib
    downloads_dir = _PROJECT_ROOT / "downloads"
    filename = "MNOVAscripts.zip"
    file_path = downloads_dir / filename

    # Debug print
    logger.debug(f"Looking for file at: {file_path}")

    # Check if the file exists
    if not file_path.exists():
        logger.warning(f"File not found: {file_path}")
        abort(404)

    # Serve the file
    return send_from_directory(str(downloads_dir), filename, as_attachment=True)


@bp.route("/documentation/")
def index():
    """Serve the main documentation index page.

    This endpoint returns the main HTML documentation page generated by Sphinx.

    Returns:
        Response: The index.html documentation file.
    """
    return send_from_directory(DOCS_DIR, "index.html")


@bp.route("/documentation/<path:filename>")
def docs(filename):
    """Serve a specific documentation file from the Sphinx-generated docs.

    This endpoint returns the requested documentation file from the documentation directory.

    Args:
        filename (str): The path to the documentation file to serve.

    Returns:
        Response: The requested documentation file.
    """
    return send_from_directory(DOCS_DIR, filename)


# create a set of example routes for displaying a number of final examples for showing some results
@bp.route("/example1/")
def example1():
    """Render the first example page.

    This endpoint serves the example1.html template, which can be used to demonstrate specific features or results.

    Returns:
        Response: The rendered example1.html template.
    """
    return render_template("example1.html")


@bp.route("/example2/")
def example2():
    """Render the second example page.

    This endpoint serves the example2.html template, which can be used to demonstrate specific features or results.

    Returns:
        Response: The rendered example2.html template.
    """
    return render_template("example2.html")


@bp.route("/simulatedAnnealingDemo/")
def example3():
    """Render the simulated annealing demo page.

    This endpoint serves the run_minimization_v4.html template, which can be used to demonstrate the simulated annealing process.

    Returns:
        Response: The rendered run_minimization_v4.html template.
    """
    return render_template("run_minimization_v4.html")


# create a route that accepts a molfile string via post and calculates the C13 ppm predictions using nmrshiftDB and returns the predictions as a JSON object
@bp.route("/predict_c13_shifts/", methods=["POST"])
def predict_c13_shifts():
    """Predict C13 NMR chemical shifts from a given molfile.

    This endpoint accepts a molfile string via POST request, processes it using nmrshiftDB,
    and returns the predicted C13 chemical shifts as a JSON object.

    Returns:
        Response: A JSON object containing the predicted C13 chemical shifts or an error message.
    """

    if request.method != "POST":
        return "Only POST requests are accepted", 400

    json_data = request.get_json()
    if json_data is None or "molstring" not in json_data:
        logger.warning("No molfile data received")
        return "No molfile data received", 400

    molfile = json_data["molstring"]

    # Use nmrshiftDB to predict C13 shifts
    predicted_shifts = expectedmolecule.calc_c13_chemical_shifts_using_nmrshift2D(
        molfile
    )

    # convert pandas dataframe to dictionary include AtomIndex in the output which is the index of the dataframe
    predicted_shifts_dict = predicted_shifts.to_dict(orient="index")

    return {"predicted_c13_shifts": predicted_shifts_dict}


@bp.route("/", methods=["GET"])
def display_front_page():
    """Render and serve the main front page of the application.

    This endpoint displays the mainPage.html template as the application's front page.

    Returns:
        Response: The rendered mainPage.html template.
    """
    return render_template("mainPage.html")


@bp.route("/simpleMNOVAfinalHTML", methods=["POST"])
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
            logger.warning("No JSON data received")
            return "No JSON data received", 400

    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON: {e}")
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
            logger.warning(f"atomNumber {atom_number} not found in catoms_orig_lookup")

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

    hostname = jsonUtils.extract_hostname(
        json_data["oldjsondata"]
    )  # returns a string or None

    # Check if device is registered and registration is valid
    if not is_device_registered(hostname):
        return jsonify(
            {
                "status": "unregistered",
                "registration_url": url_for("registration_page", _external=True)
                + f"?hostid={hostname}",
            }
        )

    if has_device_expired(hostname):
        email = get_device_stats(hostname)["email"]
        return jsonify(
            {
                "status": "registration_expired",
                "registration_url": url_for("registration_page", _external=True)
                + f"?hostid={hostname}"
                + f"&email={email}",
            }
        )

    try:
        machine_learning_opt_in = json_data["oldjsondata"]["ml_consent"]["data"]["0"]
    except KeyError:
        # If the key doesn't exist, set machine_learning_opt_in to False
        machine_learning_opt_in = False

    updated_count = update_ml_consent_for_all_user_devices(
        hostname, machine_learning_opt_in
    )

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
            with current_app.app_context():
                # MODIFIED: removed if RUNNINGONPYTHONANYWHERE branch — both the
                # MySQL and SQLite branches created an identical Result object, so
                # the check was doing nothing. A single branch is used now.
                new_result = Result(
                    user_id=user_id,
                    smiles_string=json_data["smilesString"],
                    weight=json_data["best_results"]["best_weight"],
                    MAE=json_data["best_results"]["best_mae"],
                    LAE=json_data["best_results"]["best_lae"],
                    json_result=json_result,
                )
                db.session.add(new_result)
                db.session.commit()

    return rtn_html


@bp.route("/check_machine_learning", methods=["POST"])
def check_machine_learning():
    """Check the machine learning consent status for a registered device.

    This endpoint verifies device registration and expiration, then returns the machine learning consent status.
    It responds with registration status and a registration URL if the device is unregistered or expired.

    Returns:
        Response: A JSON object with registration status and machine learning consent information.
    """

    json_data = request.get_json()

    hostname = json_data.get("hostname")

    # Check if device is registered and registration is valid
    if not is_device_registered(hostname):
        return jsonify(
            {
                "status": "unregistered",
                "registration_url": url_for("registration_page", _external=True)
                + f"?hostid={hostname}",
            }
        )

    if has_device_expired(hostname):
        email = get_device_stats(hostname)["email"]
        return jsonify(
            {
                "status": "registration_expired",
                "registration_url": url_for("registration_page", _external=True)
                + f"?hostid={hostname}"
                + f"&email={email}",
            }
        )

    # get machine learning consent status
    ml_consent = get_device_stats(hostname)["ml_consent"]
    return jsonify({"status": "registered", "ml_consent": ml_consent})


@bp.route("/register", methods=["POST"])
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

    if not hostid or not email:
        return "Both Host ID and email are required", 400

    if not license_agreed:
        return "You must agree to the license agreement", 400

    # Register the device in the database with consent status
    success = register_device(email, hostid, license_agreed, ml_consent)

    if success:
        return render_template(
            "registration_success.html",
            hostid=hostid,
            email=email,
            ml_consent=ml_consent,
        )
    else:
        return render_template("registration_error.html")


@bp.route("/registration", methods=["GET"])
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


@bp.route("/simpleMNOVA", methods=["POST"])
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
        if json_data is None:
            return "No JSON data received", 400

    except json.JSONDecodeError as e:
        return f"Invalid JSON: {e}", 400

    # check if hostname is already known
    # print out the hostname

    hostname = jsonUtils.extract_hostname(json_data)  # returns a string or None

    # Check if device is registered and registration is valid
    if not is_device_registered(hostname):
        return jsonify(
            {
                "status": "unregistered",
                "registration_url": url_for("registration_page", _external=True)
                + f"?hostid={hostname}",
            }
        )

    if has_device_expired(hostname):
        email = get_device_stats(hostname)["email"]
        return jsonify(
            {
                "status": "registration_expired",
                "registration_url": url_for("registration_page", _external=True)
                + f"?hostid={hostname}"
                + f"&email={email}",
            }
        )

    # Record this usage
    record_usage(hostname)

    # check if the user has agreed to the machine learning consent
    try:
        machine_learning_opt_in = json_data["ml_consent"]["data"]["0"]
    except KeyError:
        # If the key doesn't exist, set machine_learning_opt_in to False
        machine_learning_opt_in = False

    updated_count = update_ml_consent_for_all_user_devices(
        hostname, machine_learning_opt_in
    )

    problemdata_json = NMRProblem.from_mnova_dict(json_data)

    # decide whether we are doing prediction or assignments
    if problemdata_json.is_prediction():

        solution = nmrsolution.NMRsolution(problemdata_json)

        if solution.nmrsolution_failed:
            return solution.nmrsolution_error_message, solution.nmrsolution_error_code

        if solution.expected_molecule.nmrshiftdb_failed:
            return solution.solution_error_message, solution.solution_error_code

        ok, msg = solution.init_class_from_json()

        if not ok:
            return msg, 400

        rtn_msg, rtn_num = solution.assign_CH3_CH2_CH1_overall()

        df = solution.expected_molecule.molprops_df

        if rtn_msg != "ok":
            return rtn_msg, rtn_num

        solution.transfer_hsqc_info_to_c13()
        solution.transfer_hsqc_info_to_h1()

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

        jsonGraphData = _node_link_data(G2)
        jsonGraphData["moved_nodes"] = jsonGraphData["nodes"]

        # create a network graph of the expected molecule

        solution.initiate_molgraph(json_data, G2)

        #  convert the molgraph to a json object
        jsonGraphData_mol = _node_link_data(solution.molgraph)

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
                "sym_atom_idx",
                "sym_atomNumber",
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

        possible_symmetry = (
            solution.c13_df.shape[0] < solution.expected_molecule.num_carbon_atoms
        )

        simAnneal = SimulatedAnnealing2.from_params(
            copy.deepcopy(jsonGraphData["nodes"]),
            copy.deepcopy(jsonGraphData["links"]),
            json_data["molfile"]["data"]["0"],
            json_data,
            possible_symmetry,
        )

        # setup the run using standard parameters hardcoded in the class
        simAnneal.setup_run(cooling_rate=0.999)

        if json_data["simulatedAnnealing"]["data"]["0"] and (
            simAnneal.predicted_weight > 0
        ):  # True

            # simAnneal.setup_run(randomize_mapping=True)
            simAnneal.run_optimization(100)

            jsonGraphData, best_results = simAnneal.process_results(
                catoms_df, jsonGraphData
            )

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
        float_pattern = r"np\.float64\(([\d\.]+)\)"
        int_pattern = r"np\.int64\(([\d]+)\)"
        rtn_html = re.sub(float_pattern, r"\1", rtn_html)
        rtn_html = re.sub(int_pattern, r"\1", rtn_html)

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
                with current_app.app_context():
                    # MODIFIED: removed if RUNNINGONPYTHONANYWHERE branch — both the
                    # MySQL and SQLite branches created an identical Result object, so
                    # the check was doing nothing. A single branch is used now.
                    new_result = Result(
                        user_id=user_id,
                        smiles_string=solution.smilesstr,
                        weight=best_results["best_weight"],
                        MAE=best_results["best_mae"],
                        LAE=best_results["best_lae"],
                        json_result=json_result,
                    )
                    db.session.add(new_result)
                    db.session.commit()
        return rtn_html
    else:
        return msg, 400
