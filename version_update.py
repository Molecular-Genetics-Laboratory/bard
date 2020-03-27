# bumpy_id.py
# Usage:
# ~$ python bump_id.py script.py

# Update version string with current timestamp


import os
import sys
import time
from datetime import datetime as dt

def main():
    idstr = str(dt.now().strftime("%H-%M-%S-%d-%m-%Y"))
    fi    = sys.argv[1]
    fo    = "tmp.xyz"

    fhin  = open(fi, "r")
    fhout = open(fo, "w")

    for line in fhin:
        if "versionstr =" in line:
            fhout.write(
            "versionstr = \"bard v1.0 ID={}\"\n".format(idstr)
            )
        else:
            fhout.write(line)

    fhin.close()
    fhout.close()
    os.system("rm {}".format(fi))
    os.system("mv tmp.xyz {}".format(fi))


if __name__ == "__main__":
    main()
