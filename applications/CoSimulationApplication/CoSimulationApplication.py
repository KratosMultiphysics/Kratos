from __future__ import print_function, absolute_import, division

# Application dependent names and paths
import KratosMultiphysics as KM
application = 'KratosCoSimulationApplication'
application_name = 'KratosCoSimulationApplication'
application_folder = 'CoSimulationApplication'

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)