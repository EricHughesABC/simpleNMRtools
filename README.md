# simpleNMRtools

simpleNMRtools is a Python/Flask web application for analysing small molecule NMR data
to aid structure validation and elucidation. It uses RDKit for cheminformatics,
NetworkX for graph processing, JPype to bridge NMR shift predictions from a Java library
(NMRShiftDB), and Flask with SQLAlchemy for its web server and database backend.

Getting it to run locally is quite involved. My recommendation would be to try the simpleNMR version hosted at [simpleNMR](https://simplenmr.pythonanywhere.com/). At the home page one can download the MNOVA QtScripts and there are instructions to install them in the documentation that can be found at the home page.

If you still prefer to run the simpleNMR server locally, I have written a step-by-step guide below, but it is difficult to cover all computer and python configurations, so the instructions maynot be complete. Feel free to get in touch if you run into trouble.  

[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/EricHughesABC/simpleNMRtools)

---

## Local Installation

### Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or
  [Anaconda](https://www.anaconda.com/download) installed on your system

> **Note:** Git and Java are both required but neither needs to be installed system-wide.
> Git will be installed via conda in Step 1, and Java will be installed locally inside
> the conda environment in Step 3.

---

### Step 1 — Install Git and clone the repository

Git is not bundled with Anaconda or Miniconda, so the first step is to install it into
the conda base environment. Open a terminal (or Anaconda Prompt on Windows) and run:

```bash
conda install -c conda-forge git
```

Then clone the repository:

```bash
git clone https://github.com/EricHughesABC/simpleNMRtools.git
```

```bash
cd simpleNMRtools
```

---

### Step 2 — Create a conda environment with Python

Create a dedicated conda environment. Python 3.10 is recommended for compatibility
with all dependencies (RDKit, JPype, Flask):

```bash
conda create -n nmrtools python=3.10
```

```bash
conda activate nmrtools
```

---

### Step 3 — Install conda-managed packages

Install packages that are best obtained through conda. This includes RDKit (distributed
via `conda-forge`) and OpenJDK, which installs a Java Runtime Environment locally inside
the conda environment — no system-wide Java installation is needed:

```bash
conda install conda-forge::openjdk
```

```bash
conda install -c conda-forge rdkit numpy scipy pandas networkx pillow
```

The OpenJDK installed by conda is self-contained within the `nmrtools` environment and
will be automatically found by JPype when the application runs.

---

### Step 4 — Install remaining dependencies via pip

Because pyenv or other Python version managers can intercept the `pip` command via
PATH shims — causing packages to install into the wrong Python — always invoke pip
through the active environment's Python directly:

```bash
python -m pip install -r requirements.txt
```

`python -m pip` bypasses most shims and is guaranteed to use the Python belonging to the
active `nmrtools` conda environment.

To confirm you are using the right Python before installing:

```bash
which python        # macOS / Linux
where python        # Windows (Anaconda Prompt)
```

The path should contain `nmrtools` (e.g. `.../envs/nmrtools/bin/python`).

**If the path still shows a pyenv shim** (e.g. `~/.pyenv/shims/python`) even with
`(nmrtools)` shown in the prompt, pyenv is overriding conda in your shell's `PATH`.
See the Troubleshooting section for how to fix this. As an immediate workaround, find
the full path to the conda environment's Python and use it directly:

```bash
conda env list      # note the full path shown next to nmrtools
```

Then run all Python and pip commands using that path explicitly, for example:

```bash
/Users/you/miniconda3/envs/nmrtools/bin/python -m pip install -r requirements.txt
```

```bash
/Users/you/miniconda3/envs/nmrtools/bin/python simpleNMRtest_app.py
```

This installs the following key packages:

| Package | Purpose |
|---|---|
| Flask | Web server framework |
| Flask-SQLAlchemy | Database ORM integration |
| Flask-Migrate | Database schema migrations |
| Jinja2 | HTML templating |
| jpype1 | Java bridge for NMRShiftDB predictions |
| python-dotenv | Environment variable management |
| PyMySQL | MySQL driver for production deployment |
| scipy | Scientific calculations (stats, optimisation) |
| networkx | Molecular graph processing |
| sphinx / sphinx-rtd-theme | Documentation generation |

---

### Step 5 — Configure the environment

#### Local development

No `.env` file is needed for local development:

- Uses **SQLite** as the database. The file `registrations.db` is created automatically
  in the project directory on first run — no configuration required.
- Falls back to a built-in default for `SECRET_KEY`. This is fine for local use; do
  not use the default in production.

#### Production deployment on PythonAnywhere

PythonAnywhere does **not** use a `.env` file. It injects environment variables
directly into the application through its web dashboard. There is no `.env` file to
create or edit.

To configure the database, go to the **Web** tab in the PythonAnywhere dashboard,
scroll to the **Environment variables** section, and set the following variables:

| Variable | Value |
|---|---|
| `DB_USERNAME` | `simpleNMR` |
| `DB_PASSWORD` | your MySQL password |
| `DB_HOST` | `simpleNMR.mysql.pythonanywhere-services.com` |
| `DB_NAME` | `simpleNMR$registrations` |
| `SECRET_KEY` | a long random string (recommended for production) |

After saving and reloading the web app, the startup log (error log tab) will confirm
the database is active:

```
[db] Using DB_* environment variables — host: simpleNMR.mysql.pythonanywhere-services.com/simpleNMR$registrations
```

#### Production deployment on other servers

For any server that is not PythonAnywhere, create a `.env` file in the same folder as
`simpleNMRtest_app.py` and set a single `DATABASE_URL` variable:

```
SECRET_KEY=your-strong-secret-key-here
DATABASE_URL=mysql+pymysql://your_db_username:your_db_password@your_db_host/your_db_name
```

The application reads this file automatically via `python-dotenv` at startup.

---

### Step 6 — Configure the MestReNova client

The client side of simpleNMRtools runs inside **MestReNova** (Mnova) as a set of
ECMAScript (`.qs`) scripts. Before using them, you must tell the scripts which server
to connect to by editing the file `server_address.qs`.

Open `server_address.qs` in a text editor. The file contains a `server_address()`
function with several return statements, all but one commented out:

```javascript
function server_address() {
    // return "http://localhost:8000/";
    // return "http://localhost:5000/";
    return "http://simplenmr.pythonanywhere.com/";
    // return "http://test-simplenmr.pythonanywhere.com/";
    // return "http://simplenmr.awh.durham.ac.uk/";
}
```

Only one line must be active (uncommented) at a time. Comment out the current active
line and uncomment the one that matches your setup:

| Scenario | Line to uncomment |
|---|---|
| Local development (Flask default port) | `return "http://localhost:5000/";` |
| Local development (alternative port) | `return "http://localhost:8000/";` |
| Production (PythonAnywhere) | `return "http://simplenmr.pythonanywhere.com/";` |
| Test server (PythonAnywhere) | `return "http://test-simplenmr.pythonanywhere.com/";` |

For example, to use the local Flask server, the file should look like this:

```javascript
function server_address() {
    // return "http://localhost:8000/";
    return "http://localhost:5000/";
    // return "http://simplenmr.pythonanywhere.com/";
    // return "http://test-simplenmr.pythonanywhere.com/";
    // return "http://simplenmr.awh.durham.ac.uk/";
}
```

Save the file and reload the scripts in MestReNova for the change to take effect.

---

### Step 7 — Run the application

Start the Flask development server:

```bash
python simpleNMRtest_app.py
```

The application will be available at `http://localhost:5000`.

---

### Notes on Java and NMRShiftDB

The NMR shift prediction module (`javaUtils.py`) uses JPype to call a Java-based
NMRShiftDB prediction tool. The OpenJDK installed via conda in Step 3 is scoped to the
`nmrtools` environment, so no system-wide Java is required. Ensure the required JAR
files are present in the `lib/` subdirectory of the project.

To verify that Java is available inside the active environment:

```bash
java -version
```

If this fails, confirm the `nmrtools` environment is active (`conda activate nmrtools`)
and that the OpenJDK install completed without errors.

---

### Troubleshooting

- **RDKit import errors**: Ensure you installed RDKit via conda (`conda install -c conda-forge rdkit`),
  not pip, as the pip builds can be incomplete on some platforms.
- **JPype / JVM errors**: Confirm the `nmrtools` conda environment is active and that `conda install conda-forge::openjdk` completed successfully. Since Java is installed locally by conda, no system `JAVA_HOME` configuration should be needed.
- **MySQL errors on local development**: These can be safely ignored — if no `DATABASE_URL` is set in the environment, the app falls back to SQLite automatically.
- **Package conflicts**: If any package fails to install, try running
  `python -m pip install <package> --no-deps` and then resolve dependencies individually.
- **Activate environment**: Always ensure the `nmrtools` conda environment is active
  before running any commands (`conda activate nmrtools`).
- **pyenv overriding conda** (macOS): If `which python` still shows a pyenv shim even
  with `(nmrtools)` active, the pyenv initialisation block in `~/.zshrc` is placed after
  the conda block, so pyenv shims take precedence. Fix this by moving the pyenv block
  **above** the conda initialisation block in `~/.zshrc`, then run `source ~/.zshrc`.
  The pyenv block looks like this and must come first:
  ```bash
  export PYENV_ROOT="$HOME/.pyenv"
  export PATH="$PYENV_ROOT/bin:$PATH"
  eval "$(pyenv init -)"
  ```
  Until this is fixed, use the full conda environment path instead of `python` directly
  (see Step 4).

For issues with the software, raise an issue on the
[GitHub repository](https://github.com/EricHughesABC/simpleNMRtools).
