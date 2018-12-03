import sys
import os.path
kratos_libs=os.path.abspath(os.path.join(os.path.dirname(__file__),'../libs'))
kratos_applications=os.path.abspath(os.path.join(os.path.dirname(__file__),'../applications'))
kratos_scripts=os.path.abspath(os.path.join(os.path.dirname(__file__),'../kratos/python_scripts'))
kratos_tests=os.path.abspath(os.path.join(os.path.dirname(__file__),'../kratos/tests'))
sys.path.append(kratos_libs)
sys.path.append(kratos_applications)
sys.path.append(kratos_scripts)
sys.path.append(kratos_tests)
