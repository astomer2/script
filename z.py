import argparse
from pathlib import Path
import os

parser = argparse.ArgumentParser()
parser.add_argument( "-g","--gpu_id", type=str, default="0", help="gpu ids")
parser.add_argument( "-t","--task", type=Path, default="test", help="task path")
parser.parse_args()
gpu_id = parser.parse_args().gpu_id
task = parser.parse_args().task
print("testing: %s" % gpu_id, "testing: %s" % task)
