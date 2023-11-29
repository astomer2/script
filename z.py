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

sys.path.append(os.path.abspath('.'))

from openMM.caculate_potential_energy import batch_calc_and_rank_by_raw_diff_energy

list  = []
for root, dirs, files in os.walk("/mnt/nas1/lanwei-125/PRLR/MD/"):
    for file in files:
        if file.endswith('.pdb'):
            pdb_path = os.path.join(root, file)
            list.append(pdb_path)    


ic(batch_calc_and_rank_by_raw_diff_energy(list))