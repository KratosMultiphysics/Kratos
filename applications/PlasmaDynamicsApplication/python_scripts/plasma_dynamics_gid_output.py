import KratosMultiphysics as Kratos
from KratosMultiphysics import MultiFileFlag
from KratosMultiphysics import GiDPostMode

from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_gid_output import SwimmingDEMGiDOutput as SDEMGiDOutput


class PlasmaDynamicsGiDOutput(SDEMGiDOutput):

    def __init__(self,
                 file_name,
                 vol_output=True,
                 post_mode="Binary",
                 multifile="Single",
                 deformed_mesh=False,
                 write_conditions=True
                 ):
                 

        super().__init__(file_name,
                            vol_output,
                            post_mode,
                            multifile,
                            deformed_mesh,
                            write_conditions
                            )
	
	
    def initialize_plasma_dynamics_results(self, DEM_model_part, clusters_model_part, rigid_faces_model_part, mixed_model_part):

        if self.multi_file == MultiFileFlag.SingleFile:
            Kratos.Logger.PrintWarning("Singlefile option is not available for the plasma dynamics application!")
            mesh_name = 0.0
            self.io.InitializeMesh(mesh_name)
            self.io.WriteSphereMesh(DEM_model_part.GetMesh())
            self.io.WriteMesh(clusters_model_part.GetMesh())
            self.io.WriteMesh(rigid_faces_model_part.GetMesh())
            self.io.WriteMesh(mixed_model_part.GetMesh())
            self.io.FinalizeMesh()
            self.io.InitializeResults(mesh_name, mixed_model_part.GetMesh())

        # Initialize list file
        with open(self.listfilename, "w") as listfile:

            if self.multi_file == MultiFileFlag.MultipleFiles:
                listfile.write("Multiple\n")

            elif self.multi_file == MultiFileFlag.SingleFile:
                listfile.write("Single\n")

        if self.multi_file == MultiFileFlag.SingleFile:

            if self.post_mode == GiDPostMode.GiD_PostBinary:
                self.write_step_to_list(0)

            else:
                self.write_step_to_list(0)


    def write_plasma_dynamics_results(self, label,
                            fluid_model_part,
		                    DEM_model_part,
		                    clusters_model_part,
		                    rigid_faces_model_part,
		                    mixed_model_part,
		                    fluid_nodal_variables,
                            DEM_nodal_variables,
                            cluster_variables,
	                        rigid_faces_variables,
	                        mixed_nodal_variables,
	                        fluid_gp_variables):

        out_model_part = self.get_out_model_part(fluid_model_part)

        # update cut data if necessary
        if not self.volume_output:
            self.cut_app.UpdateCutData(out_model_part, fluid_model_part)

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self.io.InitializeMesh(label)
            self.io.WriteSphereMesh(DEM_model_part.GetMesh())
            self.io.WriteMesh(mixed_model_part.GetMesh())
            self.io.WriteMesh(rigid_faces_model_part.GetMesh())
            self.io.FinalizeMesh()
            self.io.InitializeResults(label, mixed_model_part.GetMesh())

        for var in fluid_nodal_variables:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, fluid_model_part, kratos_variable)

        for var in DEM_nodal_variables:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, DEM_model_part, kratos_variable)

        for var in cluster_variables:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, clusters_model_part, kratos_variable)

        for var in rigid_faces_variables:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, rigid_faces_model_part, kratos_variable)

        for var in mixed_nodal_variables:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, mixed_model_part, kratos_variable)

        for var in fluid_gp_variables:
            kratos_variable = Kratos.KratosGlobals.GetVariable(var)
            self._write_gp_results(label, fluid_model_part, kratos_variable)

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self._finalize_results()

            self.write_step_to_list(label)

            self.write_step_to_outer_list(label)

