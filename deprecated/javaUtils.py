import subprocess
import os
import platform
import pandas as pd
from io import StringIO
from functools import lru_cache

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# JAVAPARAMS = (
#         "java",
#         "-classpath",
#         "predictorc.jar:cdk-2.7.1.jar:.",
#         "MolFileParser",
#         ""
# )


# @lru_cache(maxsize=None)
def calculate_nmrpredictions_nmrshiftdb(molstr: str) -> pd.DataFrame:
    """returns carbon chemical shift prediction using nmrshiftdb

    Args:
        javaparams (List): parameters to call Java nmrshittdb c13prectictor using subprocess
        molstr (str): mol string of molecule

    Returns:
        pd.DataFrame: carbon 13 chemical shift prediction with columns min, mean, max or empty dataframe if an error occurred
    """
    #             r"jre\\javawindows\\bin\\java",

    print("platform.system()", platform.system())
    if platform.system() == "Windows":
        print("Windows")
        javaparams = [
            r"jre\javawindows\bin\java.exe",
            "-classpath",
            "predictorc.jar;cdk-2.7.1.jar;.",
            "MolFileParser",
            "",
        ]
    else:
        print("Not Windows")
        javaparams = [
            "java",
            "-classpath",
            "predictorc.jar:cdk-2.7.1.jar:.",
            "MolFileParser",
            "",
        ]

    # print("\njavaparams\n")
    # print(javaparams)

    javaparams[-1] = molstr

    result = subprocess.run(
        javaparams, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
    )
    # print("result\n", result)
    # print("calculate_nmrpredictions_nmrshiftdb(  molstr: str )-> pd.DataFrame: CALLED")
    # print("result.returncode", result.returncode)
    # print("result.stdout", result.stdout)
    # print("result.stderr", result.stderr)
    # Check if the command was successful
    if result.returncode != 0:
        return pd.DataFrame(columns=["min", "mean", "max"])
    else:
        print("return code == 0")
        output = result.stdout.decode()
        # print("output", output)
        try:
            mol_df = pd.read_csv(StringIO(output), index_col=0)
            mol_df.index = mol_df.index - 1
            return mol_df
        except:
            return pd.DataFrame(columns=["min", "mean", "max"])


if __name__ == "__main__":

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

    molstr2 = """
JME 2024-04-29 Mon Sep 16 20:41:40 GMT+100 2024

 18 20  0  0  0  0  0  0  0  0999 V2000
    1.9124    7.5124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    6.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5124    7.5124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    5.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    4.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    3.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.4000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4249    1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4249    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6373    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8498    1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0622    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0622    3.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8498    4.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6373    3.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4249    4.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4249    5.6000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  3  2  1  0  0  0  0
  4  2  1  0  0  0  0
  5  4  1  0  0  0  0
  6  5  1  0  0  0  0
  7  6  1  0  0  0  0
  8  7  2  0  0  0  0
  9  7  1  0  0  0  0
 10  9  2  0  0  0  0
 11  9  1  0  0  0  0
 12 11  2  0  0  0  0
 13 12  1  0  0  0  0
 14 13  2  0  0  0  0
 15 14  1  0  0  0  0
 16 15  2  0  0  0  0
 16 11  1  0  0  0  0
 17 16  1  0  0  0  0
 17  6  2  0  0  0  0
 18 17  1  0  0  0  0
 18  2  1  0  0  0  0
M  END

"""

    molstr3 = """
  Mrv2011 09232414502D

 10 11  0  0  0  0  0  0  0  0999 V2000
    3.0480   -0.8048   -0.3219 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0882   -1.5427   -0.1718 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9764   -0.0278   -0.1718 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4378    1.0900   -0.2919 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1365    1.7665    0.0292 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7438   -1.0832    0.7376 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7233   -0.1862    0.0190 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9324    0.9055    0.6303 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6265   -0.5239    0.9877 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5391   -0.2170   -0.1718 C   0  0  0  0  0  0  0  0  0  0  0  0
  5  4  1  0  0  0  0
  6  7  1  0  0  0  0
  7  1  1  0  0  0  0
  7  4  2  0  0  0  0
  8  5  1  0  0  0  0
  9  6  1  0  0  0  0
  9  8  1  0  0  0  0
 10  2  1  0  0  0  0
 10  3  1  0  0  0  0
 10  8  1  0  0  0  0
 10  9  1  0  0  0  0
M  END

"""

    molstr4 = """
 Mnova   09232416552D

 10 11  0  0  0  0  0  0  0  0999 V2000
   33.5340   16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   48.0542    8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   48.0542   -8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   62.5761  -16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   77.0963   -8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   77.0963    8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   62.5761   16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   91.6165    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  100.0000   14.5202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  100.0000  -14.5202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  5  4  1  6  0  0  0
  5  6  1  0  0  0  0
  6  7  1  6  0  0  0
  7  2  1  0  0  0  0
  6  8  1  0  0  0  0
  8  5  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
M  END
"""
    molstr6 = """
JME 2024-04-29 Mon Sep 23 21:25:06 GMT+100 2024

 10 11  0  0  0  0  0  0  0  0999 V2000
    0.0000    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4249    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6373    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6373    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4249    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8497    1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5498    2.6124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5498    0.1876    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  5  4  1  6  0  0  0
  5  6  1  0  0  0  0
  6  7  1  6  0  0  0
  7  2  1  0  0  0  0
  6  8  1  0  0  0  0
  8  5  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
M  END
"""
    # M  ZZC   1 10
    # M  ZZC   2 6
    # M  ZZC   3 5
    # M  ZZC   4 4
    # M  ZZC   5 3
    # M  ZZC   6 2
    # M  ZZC   7 7
    # M  ZZC   8 1
    # M  ZZC   9 8
    # M  ZZC  10 9

    molstr7 = "\r\nMnova   09232416552D\r\n\r\n 10 11  0  0  0  0  0  0  0  0999 V2000\r\n   33.5340   16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   48.0542    8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   48.0542   -8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   62.5761  -16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   77.0963   -8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   77.0963    8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   62.5761   16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n   91.6165    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n  100.0000   14.5202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n  100.0000  -14.5202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n  1  2  1  0  0  0  0\r\n  2  3  2  0  0  0  0\r\n  3  4  1  0  0  0  0\r\n  5  4  1  6  0  0  0\r\n  5  6  1  0  0  0  0\r\n  6  7  1  6  0  0  0\r\n  7  2  1  0  0  0  0\r\n  6  8  1  0  0  0  0\r\n  8  5  1  0  0  0  0\r\n  8  9  1  0  0  0  0\r\n  8 10  1  0  0  0  0\r\nM  ZZC   1 10\r\nM  ZZC   2 6\r\nM  ZZC   3 5\r\nM  ZZC   4 4\r\nM  ZZC   5 3\r\nM  ZZC   6 2\r\nM  ZZC   7 7\r\nM  ZZC   8 1\r\nM  ZZC   9 8\r\nM  ZZC  10 9\r\nM  END\r\n"
    molstr5 = "\nMnova   09232416552D\n\n 10 11  0  0  0  0  0  0  0  0999 V2000\n   33.5340   16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   48.0542    8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   48.0542   -8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   62.5761  -16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   77.0963   -8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   77.0963    8.3835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   62.5761   16.7670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   91.6165    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  100.0000   14.5202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  100.0000  -14.5202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  3  4  1  0  0  0  0\n  5  4  1  6  0  0  0\n  5  6  1  0  0  0  0\n  6  7  1  6  0  0  0\n  7  2  1  0  0  0  0\n  6  8  1  0  0  0  0\n  8  5  1  0  0  0  0\n  8  9  1  0  0  0  0\n  8 10  1  0  0  0  0\nM  ZZC   1 10\nM  ZZC   2 6\nM  ZZC   3 5\nM  ZZC   4 4\nM  ZZC   5 3\nM  ZZC   6 2\nM  ZZC   7 7\nM  ZZC   8 1\nM  ZZC   9 8\nM  ZZC  10 9\nM  END\n"

    print(calculate_nmrpredictions_nmrshiftdb(molstr5))
    # print( calculate_nmrpredictions_nmrshiftdb( JAVAPARAMS, "" ))
