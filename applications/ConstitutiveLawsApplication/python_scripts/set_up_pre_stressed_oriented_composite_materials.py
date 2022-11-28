import KratosMultiphysics

class SetUpPreStressedOrientedCompositeMaterials(KratosMultiphysics.Process):

    """This process sets a proper orientation of the local axes of the elements intersected by line elements (steel tendons). Besides it also computes and sets a volumetric participation of steel within the concrete FE as well as an indicated pre-stressing strain. It also creates a submodelpart for each steel tendon intersected FE.

    Format of the file:
    --> Intersection points:
    Begin Tendon - Hexahedra intersection: Tendon_Inf
	    4807	        0.01000	        0.05000	        0.01500		        0.02000	        0.05000	        0.01500
        ...
    End Tendon - Hexahedra intersection

    --> FE inside each tendon:
    Begin Tendon - Hexahedra: Tendon_Inf
        4807    4812    4817    4822    4827    4832    4837    4842    4847    4852   ...
    End Tendon - Hexahedra


    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """