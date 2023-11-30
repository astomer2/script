import argparse
import json
import logging
import math
import os
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from icecream import ic
from modeller import energy_data

list  = []
for root, dirs, files in os.walk("/mnt/sdc/lanwei/software/cPEPmatch/5w59/colab/cycol/"):
    for file in files:
        if file.endswith('.pdb'):
            list.append(file.split('.')[0])    

with open("log.txt", "a+") as f:
    for pdb in list:
        f.write(pdb + "\n")