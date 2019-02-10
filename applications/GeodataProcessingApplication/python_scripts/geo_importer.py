import KratosMultiphysics
from geo_processor import GeoProcessor

class GeoImporter( GeoProcessor ):

    ### --- overwriting the import function --- ###

    def SetGeoModelPart( self, modelPartIn ):

        KratosMultiphysics.Logger.PrintWarning("GeoImporter: SetGeoModelPart", "This function cannot be used. The model is created be reading a file.")


    ### --- function imports the nodes from a *.stl file --- ###

    def stl_import( self, stl_file_name_input ):

        self._initialize_model_part()

        print( "I will keep my hands off this one - you are the expert, Nicola." )

        self.HasModelPart = True


    ### --- function imports the nodes from a *.xyz file --- ###

    def xyz_import( self, xyz_file_name_input ):

        self._initialize_model_part()

        with open (xyz_file_name_input) as read_file:
            node_id = 1
            for row in read_file.readlines():
                row = row.split()
                if ( all( self._isfloat(n) for n in row) ):
				    # row[0] = x coordinate; row[1] = y coordinate; row[2] = z coordinate
                    X_coord, Y_coord, Z_coord = [float(coord) for coord in row[0:3]]
                    n = self.ModelPart.CreateNewNode(node_id, X_coord, Y_coord, Z_coord)
                    node_id += 1

        self.HasModelPart = True


########################################################################
# auxiliary functions
#########################################################################

    def _isfloat( self, value):
        # detecting if a conversion in float is possible
        try:
            float( value )
            return True
        except ValueError:
            return False


    def _initialize_model_part( self ):
        # this function adds all variables that will be needed during the pipeline (please add)

        current_model = KratosMultiphysics.Model()
        self.ModelPart = current_model.CreateModelPart("ModelPart")
        self.ModelPart.SetBufferSize(1)

        self.ModelPart.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self.ModelPart.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        self.ModelPart.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)

        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FLUID)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FREE_SURFACE)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.IS_BOUNDARY)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)

        self.ModelPart.Nodes.clear()
        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()
