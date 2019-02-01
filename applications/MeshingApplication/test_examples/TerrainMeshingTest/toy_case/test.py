from scipy import spatial
import math
import numpy as np
import os
import sys

sys.path.append('../python_scripts')
from terrain_processor import terrain_processor


tool = terrain_processor()

tool.Import( "toy_terrain.xyz")

tool.MeshConcaveHull()
