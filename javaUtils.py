import jpype
import jpype.imports
from jpype.types import *
import pandas as pd
from io import StringIO
from functools import lru_cache
import os
import platform
import threading

# Change to script directory
os.chdir(os.path.dirname(os.path.realpath(__file__)))


class NMRShiftDBBridge:
    """Singleton class to manage JPype JVM and MolFileParser"""

    _instance = None
    _lock = threading.Lock()

    def __new__(cls):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        if self._initialized:
            return

        self.parser = None
        self.is_ready = False
        self._initialize_jvm()

    def _initialize_jvm(self):
        """Initialize JVM and load MolFileParser class"""
        if jpype.isJVMStarted():
            print("JVM already started")
            self._load_parser()
            return

        try:
            # Build classpath based on platform
            current_dir = os.path.dirname(os.path.abspath(__file__))
            lib_dir = os.path.join(current_dir, "lib")
            build_dir = os.path.join(current_dir, "build")

            # Collect all paths that might contain classes
            classpath_parts = []

            # Add lib directory if it exists
            if os.path.exists(lib_dir):
                # Add all JAR files in lib directory
                for file in os.listdir(lib_dir):
                    if file.endswith('.jar'):
                        classpath_parts.append(os.path.join(lib_dir, file))
                print(f"Found {len([f for f in os.listdir(lib_dir) if f.endswith('.jar')])} JAR files in lib/")

            # Add build directory if it exists (for compiled .class files)
            if os.path.exists(build_dir):
                classpath_parts.append(build_dir)
                print(f"Added build directory: {build_dir}")

            # Add current directory (for .class files in root)
            classpath_parts.append(current_dir)

            # Also check for specific JAR files in current directory (fallback)
            for jar_name in ["predictorc.jar", "cdk-2.7.1.jar"]:
                jar_path = os.path.join(current_dir, jar_name)
                if os.path.exists(jar_path) and jar_path not in classpath_parts:
                    classpath_parts.append(jar_path)

            if not classpath_parts:
                raise Exception("No JAR files or class directories found. Please check lib/ directory.")

            classpath = os.pathsep.join(classpath_parts)

            print(f"Platform: {platform.system()}")
            print(f"Classpath: {classpath}")

            # Start JVM
            jpype.startJVM(
                jpype.getDefaultJVMPath(),
                f"-Djava.class.path={classpath}",
                "-Xmx2048m",  # 2GB heap for chemistry calculations
                convertStrings=False
            )

            print("JVM started successfully")
            self._load_parser()

        except Exception as e:
            print(f"Failed to initialize JVM: {e}")
            import traceback
            traceback.print_exc()
            self.is_ready = False

    def _load_parser(self):
        """Load the MolFileParser class"""
        try:
            # Load the Java class
            print("Attempting to load MolFileParser class...")
            MolFileParser = jpype.JClass('MolFileParser')
            print("MolFileParser class loaded")

            self.parser = MolFileParser()
            print("MolFileParser instance created")

            # Check if initialization was successful
            if hasattr(self.parser, 'isInitialized') and callable(self.parser.isInitialized):
                self.is_ready = self.parser.isInitialized()
                print(f"Parser isInitialized: {self.is_ready}")
            else:
                # Older version without isInitialized method
                self.is_ready = True
                print("Parser ready (no isInitialized method)")

            if self.is_ready:
                print("✓ MolFileParser initialized successfully")
            else:
                print("✗ MolFileParser failed to initialize (PredictionTool error)")

            self._initialized = True

        except Exception as e:
            print(f"Failed to load MolFileParser: {e}")
            import traceback
            traceback.print_exc()
            self.is_ready = False
            self._initialized = True

    def predict(self, molfile_string):
        """
        Predict NMR shifts for a molfile string
        Returns pandas DataFrame with columns: min, mean, max
        """
        if not self.is_ready:
            print("Parser not ready")
            return pd.DataFrame(columns=["min", "mean", "max"])

        try:
            # Call Java method to get CSV string
            csv_result = self.parser.predictFromMolfileAsCSV(molfile_string)

            # Parse CSV into DataFrame
            mol_df = pd.read_csv(StringIO(str(csv_result)), index_col=0)

            # Adjust index to be 0-based (Java returns 1-based)
            mol_df.index = mol_df.index - 1

            # Rename columns to lowercase to match original function
            mol_df.columns = ["min", "mean", "max"]

            return mol_df

        except Exception as e:
            print(f"Error during prediction: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame(columns=["min", "mean", "max"])

    def get_carbon_count(self, molfile_string):
        """Quick check for number of carbon atoms"""
        if not self.is_ready:
            return 0

        try:
            return self.parser.getCarbonCount(molfile_string)
        except Exception as e:
            print(f"Error getting carbon count: {e}")
            return 0

    def shutdown(self):
        """Shutdown JVM (call this when application exits)"""
        if jpype.isJVMStarted():
            jpype.shutdownJVM()
            print("JVM shut down")


# Global bridge instance
_bridge = None
_bridge_lock = threading.Lock()


def _get_bridge():
    """Get or create the global NMRShiftDB bridge"""
    global _bridge
    if _bridge is None:
        with _bridge_lock:
            if _bridge is None:
                _bridge = NMRShiftDBBridge()
    return _bridge


@lru_cache(maxsize=None)
def calculate_nmrpredictions_nmrshiftdb(molstr: str) -> pd.DataFrame:
    """Returns carbon chemical shift prediction using nmrshiftdb

    Args:
        molstr (str): mol string of molecule

    Returns:
        pd.DataFrame: carbon 13 chemical shift prediction with columns min, mean, max
                     or empty dataframe if an error occurred
    """
    bridge = _get_bridge()
    return bridge.predict(molstr)


def get_carbon_count(molstr: str) -> int:
    """Get the number of carbon atoms in a molecule

    Args:
        molstr (str): mol string of molecule

    Returns:
        int: number of carbon atoms
    """
    bridge = _get_bridge()
    return bridge.get_carbon_count(molstr)


def shutdown_nmr_bridge():
    """Shutdown the JVM (call when application exits)"""
    global _bridge
    if _bridge is not None:
        _bridge.shutdown()


# Register shutdown hook
import atexit
atexit.register(shutdown_nmr_bridge)


if __name__ == "__main__":

    import time

    molstr1 = """
RDKit          2D

 20 22  0  0  0  0  0  0  0  0999 V2000
    2.6401   -3.1891    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2135   -2.7256    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    2.5981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    3.8971    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    5.1962    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    6.4952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    3.8971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2135   -2.7256    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4271   -1.8439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4271   -3.6073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9635   -5.0339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4635   -5.0339    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4182   -6.2474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -3.6073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  1  0
  8  9  1  0
  9 10  1  0
  8 11  2  0
  6 12  1  0
 12 13  2  0
 13 14  1  0
 14 15  1  0
 14 16  1  0
 16 17  1  0
 17 18  1  0
 18 19  1  0
 18 20  1  0
 20  2  1  0
 13  3  1  0
 20 14  1  0
M  END
"""

    molstr5 = "\nMnova   09232416552D\n\n 10 11  0  0  0  0  0  0  0  0999 V2000\n   33.5340   16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   48.0542    8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   48.0542   -8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   62.5761  -16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   77.0963   -8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   77.0963    8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   62.5761   16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   91.6165    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  100.0000   14.5202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  100.0000  -14.5202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  3  4  1  0  0  0  0\n  5  4  1  6  0  0  0\n  5  6  1  0  0  0  0\n  6  7  1  6  0  0  0\n  7  2  1  0  0  0  0\n  6  8  1  0  0  0  0\n  8  5  1  0  0  0  0\n  8  9  1  0  0  0  0\n  8 10  1  0  0  0  0\nM  ZZC   1 10\nM  ZZC   2 6\nM  ZZC   3 5\nM  ZZC   4 4\nM  ZZC   5 3\nM  ZZC   6 2\nM  ZZC   7 7\nM  ZZC   8 1\nM  ZZC   9 8\nM  ZZC  10 9\nM  END\n"

    print("=" * 60)
    print("Testing NMR Predictions via JPype")
    print("=" * 60)

    print("\nTest 1: molstr1")
    print("-" * 60)
    start = time.time()
    result1 = calculate_nmrpredictions_nmrshiftdb(molstr1)
    end = time.time()
    print(result1)
    print(f"Cached call took: {(end-start)*1000:.2f} ms")

    print("\nTest 2: molstr5")
    print("-" * 60)
    start = time.time()
    result3 = calculate_nmrpredictions_nmrshiftdb(molstr5)
    end = time.time()
    print(result3)
    print(f"Cached call took: {(end-start)*1000:.2f} ms")

    print("\nTest 3: Carbon count for molstr5")
    print("-" * 60)
    carbon_count = get_carbon_count(molstr5)
    print(f"Carbon atoms: {carbon_count}")

    print("\nTest 4: Cache test (calling molstr1 again - should be instant)")
    print("-" * 60)

    start = time.time()
    result1_cached = calculate_nmrpredictions_nmrshiftdb(molstr1)
    end = time.time()
    print(f"Cached call took: {(end-start)*1000:.2f} ms")
    print(result1_cached)

    print("\n" + "=" * 60)
    print("All tests completed")
    print("=" * 60)


