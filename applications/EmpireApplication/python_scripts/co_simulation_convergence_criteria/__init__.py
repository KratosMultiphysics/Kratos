# Adding the current directory to the path such that the modules can be imported with __import__ in the factory
import sys, os
sys.path.append(os.path.dirname(__file__))