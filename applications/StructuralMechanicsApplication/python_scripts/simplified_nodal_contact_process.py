# Importing the Kratos Library
import KratosMultiphysics

import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return SimplifiedNodalContactProcess(Model, settings["Parameters"])

class SimplifiedNodalContactProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        self.settings = settings;
        self.Model = Model


    def ExecuteInitialize(self):
        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "help"                       :"This process computes a simplified contact problem nodally using the distance to define the penalty value",
                "background_domain"          : "name_of_background_model_part",
                "background_contact_surface" : "name_of_surface_of_contact",
                "background_contact_volume"  : "inside_of_the_contact_volume",
                "active_contact_surface"     : "active_contact_surface",
                "active_contact_body"        : "active_contact_body",
                "contact_property_id"        : 1
            }
            """
        );

        self.settings.ValidateAndAssignDefaults(default_settings)

        KratosMultiphysics.Process.__init__(self)

        #modelparts in the background (to be used in detecting contact and computing distances to the wall)
        self.background_domain                = self.Model[self.settings["background_domain"].GetString()]
        self.background_contact_surface       = self.Model[self.settings["background_contact_surface"].GetString()]
        self.background_contact_volume        = self.Model[self.settings["background_contact_volume"].GetString()]
        self.background_all = self.background_domain

        #modelparts on the structure
        self.active_contact_body            = self.Model[self.settings["active_contact_body"].GetString()]
        self.active_contact_surface       = self.Model[self.settings["active_contact_surface"].GetString()]

        self.domain_size = self.background_domain.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        zero = KratosMultiphysics.Vector(3)
        zero[0] = 0.0
        zero[1] = 0.0
        zero[2] = 0.0

        ################## STRUCTURAL SIDE
        #compute the normals for the structure
        KratosMultiphysics.BodyNormalCalculationUtils().CalculateBodyNormals(self.active_contact_body,self.domain_size)

        max_cond_id = 1
        for cond in self.active_contact_body.Conditions:
            if(cond.Id > max_cond_id):
                max_cond_id = cond.Id

        prop_id = self.settings["contact_property_id"].GetInt()
        prop = self.active_contact_body.Properties[prop_id]

        if not prop.Has(KratosMultiphysics.YOUNG_MODULUS):
            raise Exception("property with Id ",prop_id," does not define YOUNG_MODULUS and cannot be used in contact")


        if(self.domain_size == 2):
            self.locate_on_background = KratosMultiphysics.BinBasedFastPointLocator2D(self.background_all)
            contact_condition_type = "PointContactCondition2D1N"
        else:
            self.locate_on_background = KratosMultiphysics.BinBasedFastPointLocator3D(self.background_all)
            contact_condition_type = "PointContactCondition3D1N"
        self.locate_on_background.UpdateSearchDatabase()

        i = max_cond_id + 1
        for node in self.active_contact_surface.Nodes: #nodes on the skin of the rotor
            self.active_contact_surface.GetRootModelPart().GetSubModelPart("computing_domain").CreateNewCondition(contact_condition_type,i,[node.Id], prop)
            i+=1
            node.Fix(KratosMultiphysics.DISTANCE)
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,0.0)
            node.SetValue(KratosMultiphysics.DISTANCE,0.0)
            node.SetValue(KratosMultiphysics.DISPLACEMENT,zero)


        for node in self.active_contact_body.Nodes:
            if(not node.IsFixed(KratosMultiphysics.DISTANCE)):
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,-1.0)
                print(node.Id)



        ################## BACKGROUND SIDE

        ##assigning the distances
        for node in self.background_contact_surface.Nodes: #nodes on the contact surface
            node.Fix(KratosMultiphysics.DISTANCE)
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,0.0)

        for node in self.background_domain.Nodes:
            if(node.IsFixed(KratosMultiphysics.DISTANCE) == False):
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,-1.0)

        for node in self.background_contact_volume.Nodes: #nodes inside the "air"
            if(node.IsFixed(KratosMultiphysics.DISTANCE) == False):
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,1.0)




        #computing distance from the contact surface on the background mesh
        distance_linear_solver_settings = KratosMultiphysics.Parameters( """{
                                       "solver_type" : "amgcl"
                                   } """)
        distance_linear_solver = linear_solver_factory.ConstructSolver(distance_linear_solver_settings)

        max_iterations=30
        if(self.domain_size == 2):
            self.distance_calculator = KratosMultiphysics.VariationalDistanceCalculationProcess2D(self.background_all, distance_linear_solver, max_iterations)
        else:
            self.distance_calculator = KratosMultiphysics.VariationalDistanceCalculationProcess3D(self.background_all, distance_linear_solver, max_iterations)
        self.distance_calculator.Execute()

        for node in self.background_all.Nodes:
            node.SetValue(KratosMultiphysics.DISTANCE, node.GetSolutionStepValue(KratosMultiphysics.DISTANCE))

        if(self.domain_size == 2):
            KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess2D(self.background_all, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA).Execute()
        else:
            KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess3D(self.background_all, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA).Execute()
        print("finished initialize")

    def ExecuteInitializeSolutionStep(self):
        KratosMultiphysics.BodyNormalCalculationUtils().CalculateBodyNormals(self.active_contact_body,self.domain_size)

        zero = KratosMultiphysics.Vector(3)
        zero[0] = 0.0
        zero[1] = 0.0
        zero[2] = 0.0

        N = KratosMultiphysics.Vector(self.domain_size+1)
        coords =  KratosMultiphysics.Array3()
        pelem =  KratosMultiphysics.Element(-1) #UGLY! here i create an empty pointer
        grad =  KratosMultiphysics.Vector(3)

        for node in self.active_contact_surface.Nodes: #nodes on the skin of the rotor

            #save the displacement
            disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            node.SetValue(KratosMultiphysics.DISPLACEMENT,disp)

            #now find if inside
            coords[0] = node.X
            coords[1] = node.Y
            coords[2] = node.Z
            found = self.locate_on_background.FindPointOnMesh(coords, N, pelem, 1000, 1e-9)

            if(found):
                d = 0.0
                grad[0] = 0.0
                grad[1] = 0.0
                grad[2] = 0.0
                k = 0
                for p in pelem.GetNodes():
                    d += N[k]*p.GetSolutionStepValue(KratosMultiphysics.DISTANCE)

                    g = p.GetValue(KratosMultiphysics.DISTANCE_GRADIENT)
                    grad[0] += N[k]*g[0]
                    grad[1] += N[k]*g[1]
                    grad[2] += N[k]*g[2]

                    k+=1
                node.SetValue(KratosMultiphysics.DISTANCE,d)
                node.SetValue(KratosMultiphysics.DISTANCE_GRADIENT,grad)
            else:
                node.SetValue(KratosMultiphysics.DISTANCE,0.0)
                node.SetValue(KratosMultiphysics.DISTANCE_GRADIENT,zero)





