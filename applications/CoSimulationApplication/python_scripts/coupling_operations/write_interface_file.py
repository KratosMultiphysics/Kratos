# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

def Create(*args):
    return WriteInterfaceFile(*args)

class WriteInterfaceFile(CoSimulationCouplingOperation):
    """This operation is used to write an interface file."""
    def __init__(self, settings, solver_wrappers, process_info=None):
        super().__init__(settings)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.execution_point = self.settings["execution_point"].GetString()
        self.file_name = self.settings["file_name"].GetString()
        self.model_part_name = self.settings["output_parameters"]["model_part_name"].GetString()
        # self.base_output_file_name = "{}_{}_{}_".format(self.settings["solver"].GetString(), self.model_part_name, self.execution_point)
        # number_of_points = self.model["Structure"].GetSubModelPart(self.model_part_name).NumberOfNodes()
        # print('NUMBER OF NODES = ', number_of_points)

        available_execution_points = [
            "finalize_coupling_iteration"
        ]

        if self.execution_point not in available_execution_points:
            err_msg  = 'Execution point "{}" is not available, only the following options are available:\n    '.format(self.execution_point)
            err_msg += "\n    ".join(available_execution_points)
            raise Exception(err_msg)

        self.step = 0 # this should come from self.process_info
        # TODO check if restarted. If not delete the folder => check self.process_info
        #self.output = KM.VtkOutput(self.model[model_part_name], self.settings["output_parameters"]) # currently hardcoded to vtk

    def FinalizeCouplingIteration(self):
        # if self.execution_point == "finalize_coupling_iteration":
        #     output_file_name = self.base_output_file_name + "{}_{}".format(self.step, self.coupling_iteration)
        #     self.output.PrintOutput(output_file_name)

        # upper_membrane = self._GetSolver().GetComputingModelPart().GetSubModelPart('PointLoad3D_UpperMembrane')
        '''
        number_of_points = self.model[self.model_part_name].NumberOfNodes()
        print('NUMBER OF NODES = ', number_of_points)

        print('Writing interface into: ', self.file_name)
        from scipy.io import netcdf
        ncf = netcdf.netcdf_file(self.file_name, 'w')
        nops = 'no_of_points'

        ncf.createDimension(nops, number_of_points)

        # define variables
        global_node_ids = ncf.createVariable('global_id', 'i', (nops,))
        nodal_x_coordinates = ncf.createVariable('x', 'd', (nops,))
        nodal_y_coordinates = ncf.createVariable('y', 'd', (nops,))
        nodal_z_coordinates = ncf.createVariable('z', 'd', (nops,))
        nodal_x_displacements = ncf.createVariable('dx', 'd', (nops,))
        nodal_y_displacements = ncf.createVariable('dy', 'd', (nops,))
        nodel_z_displacements = ncf.createVariable('dz', 'd', (nops,))

        i = 0
        for node in self.model[self.model_part_name].Nodes:
            global_node_ids[i] = node.Id
            nodal_x_coordinates[i] = node.X
            nodal_y_coordinates[i] = node.Y
            nodal_z_coordinates[i] = node.Z
            nodal_x_displacements[i] = node.GetSolutionStepValue(KM.DISPLACEMENT_X)
            nodal_y_displacements[i] = node.GetSolutionStepValue(KM.DISPLACEMENT_Y)
            nodel_z_displacements[i] = node.GetSolutionStepValue(KM.DISPLACEMENT_Z)
            i += 1
        ncf.close()
        #'''

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"            : "UNSPECIFIED",
            "execution_point"   : "UNSPECIFIED",
            "output_format"     : "vtk",
            "file_name"         : "UNSPECIFIED",
            "output_parameters" : { }
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultSettings())
        return this_defaults