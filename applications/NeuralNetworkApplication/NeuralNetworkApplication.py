# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosNeuralNetworkApplication import *
application = KratosNeuralNetworkApplication()
application_name = "KratosNeuralNetworkApplication"

_ImportApplication(application, application_name)
