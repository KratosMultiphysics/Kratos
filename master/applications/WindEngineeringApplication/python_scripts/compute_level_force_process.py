# Kratos imports
import KratosMultiphysics
import KratosMultiphysics.WindEngineeringApplication as WindEngineering
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

# STL imports
import pathlib


def Factory(parameters, Model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("Expecting a Parameters object containing a JSON string, but got {} instead".format(type(parameters)))
    return ComputeLevelForceProcess(Model, parameters["Parameters"])


class ComputeLevelForceProcess(KratosMultiphysics.Process):
    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters):
        """Reduce nodal reaction forces and torques on stacked slab domains.
        A region of space between 'bottom_point' and 'top_point' is subdivided into 'number_of_slabs'
        parallel slabs. Then, nodes from the specified model part are sorted into sub model parts
        based on which slab they are located in. Finally, for each sub model part, the reaction forces
        are summed up, and their torque (plus MOMENT if applicable) is reduced to 'moment_reference_point'.
        The reduced values are written to output files for each sub model part.
        Default parameters:
        {
            "model_part_name"           : "",
            "moment_reference_point"    : [0.0, 0.0, 0.0],
            "bottom_point"              : [0.0, 0.0, 0.0],
            "top_point"                 : [0.0, 0.0, 0.0],
            "number_of_slabs"           : 1,
            "open_domain"               : false,
            "time_domain"               : [0.0, 1e100],
            "output_name_stub"          : "slab_"
        }"""

        KratosMultiphysics.Process.__init__(self)
        self.model_part = model[parameters["model_part_name"].GetString()]

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.moment_reference_point = parameters["moment_reference_point"].GetVector()
        self.bottom_point = parameters["bottom_point"].GetVector()
        self.top_point = parameters["top_point"].GetVector()
        self.time_domain = parameters["time_domain"].GetVector()
        self.number_of_slabs = parameters["number_of_slabs"].GetInt()
        self.is_open_domain = parameters["open_domain"].GetBool()
        self.output_name_stub = pathlib.Path(parameters["output_name_stub"].GetString())


    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.time_domain[0] <= time and time <= self.time_domain[1]:
            sub_model_parts = WindEngineering.ModelSubdivisionUtilities.SortNodesBySlabs(
                self.model_part,
                self.bottom_point,
                self.top_point,
                self.number_of_slabs,
                self.is_open_domain)

            for model_part_index, model_part in enumerate(sub_model_parts):
                force, torque = KratosMultiphysics.ForceAndTorqueUtils.ComputeEquivalentForceAndTorque(
                    model_part,
                    self.moment_reference_point,
                    KratosMultiphysics.REACTION,
                    KratosMultiphysics.MOMENT)

                output_params = KratosMultiphysics.Parameters("""
                {
                    "output_path" : "",
                    "file_name" : ""
                }
                """)

                output_params["output_path"].SetString(str(self.output_name_stub.parent))
                output_params["file_name"].SetString("{}{}".format(self.output_name_stub.name, str(model_part_index)))

                # TODO: this is a barebone output, implement something more user friendly
                output_header = ""

                output_file = TimeBasedAsciiFileWriterUtility(
                    model_part,
                    output_params,
                    output_header
                ).file

                output_file.write("force:\t{}\ntorque:\t{}".format(force, torque))
                output_file.flush()
                output_file.close()


    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""
        {
            "model_part_name"           : "",
            "moment_reference_point"    : [0.0, 0.0, 0.0],
            "bottom_point"              : [0.0, 0.0, 0.0],
            "top_point"                 : [0.0, 0.0, 0.0],
            "number_of_slabs"           : 1,
            "open_domain"               : false,
            "time_domain"               : [0.0, 1e100],
            "output_name_stub"          : "slab_"
        }
        """)


    @staticmethod
    def ParseOutput(fileName: pathlib.Path):
        """Get output values from a file written by this process (temporary implementation)."""

        def ParseList(line: str):
            begin = line.index('(') + 1
            end = line.index(')')
            values = [float(component.strip()) for component in line[begin:end].split(',') if component ]
            return values

        force = []
        torque = []
        with open(fileName, 'r') as file:
            for line in file:
                if "force" in line:
                    force = ParseList(line)
                elif "torque" in line:
                    torque = ParseList(line)

        return force, torque 