import os

os.chdir("adaptive_mesher2d.gid")	

import remesh_build_reference.py

os.chdir("..")
os.chdir("adaptive_mesher3d.gid")	
import remesh_build_reference.py
# Add other examples here
