from Kratos import *
from KratosIncompressibleFluidApplication import *
from KratosMeshingApplication import *

class AdaptiveRemeshing:
    """ This class allows refining the problem mesh at run time.

        It will split all elements where some error estimate surpasses a given
        threshold. In 2D, some edge swapping will be performed to improve the
        quality of the refined mesh. The current version of this class is
        intended to work with VMS2D and VMS3D elements (by calling
        RefineOnSubscaleError), but it can work with other fluid elements if
        an error estimate is calculated and stored an elemental
        variable (in that case, call RefineOnErrorRatio instead).
    """

    def __init__(self,model_part,domain_size,solver):
        """ Constructor for AdaptiveRemeshing.
            model_part ModelPart containing the mesh to be refined
            domain_size Spatial dimension (2 or 3)
            solver The fluid solver, to ensure that the Dofs are regenerated
            after remeshing
        """
        self.model_part = model_part
        self.domain_size = domain_size
        self.fluid_solver = solver

        # Error estimation tools
        self.refinement_utilities = RefinementUtilities()

        # Refiner, swapping and time estimation classes
        if (domain_size == 2):
            self.refinement_process = LocalRefineTriangleMesh(model_part)
            self.swapping_process = EdgeSwapping2DModeler()
            self.time_estimator = EstimateDt2D(model_part)
        else:
            self.refinement_process = LocalRefineTetrahedraMesh(model_part)
            self.time_estimator = EstimateDt3D(model_part)

        # Neigbour search tool instance
        AvgElemNum = 10
        AvgNodeNum = 10
        self.nodal_neighbour_search = FindNodalNeighboursProcess(model_part,\
                                                                 AvgElemNum,\
                                                                 AvgNodeNum)

        # Find neighbours
        self.nodal_neighbour_search.Execute()

    def EstimateTimeStep(self,CFL,max_dt):
        """ Return the maximum time step for which the Courant number is
            lesser or equal to CFL.

            The max_dt parameter is an upper bound for situations with low
            velocities.
        """
        NewDt = self.time_estimator.EstimateDt(CFL,max_dt)
        return NewDt

    def RefineOnErrorRatio(self,refine_var,refine_tol,min_area,max_refinements):
        """ Refine all elements where the refine_var variable value is greater
            than refine_tol.

            refine_var Variable containing an error estimation
            refine_tol Maximum allowed value for the error estimate
            min_area Minimum allowed size (area or volume) for the new elements
            max_refinements Maximum refinement steps allowed from the same
            original element (to preserve mesh quality
        """

        # Clear the Dof set stored by the solver to force its regeneration
        # in the next solution step, using the new mesh
        self.fluid_solver.solver.Clear()

        # Identify refinement candidates
        self.refinement_utilities.MarkForRefinement(refine_var,\
                                                    self.model_part,\
                                                    self.domain_size,\
                                                    refine_tol,\
                                                    min_area,\
                                                    max_refinements)

        # Refine
        refine_on_reference = True
        interpolate_internal_variables = False
        self.refinement_process.LocalRefineMesh(refine_on_reference,\
                                                interpolate_internal_variables)

        # In 2D, swap edges to improve mesh quality
        if (self.domain_size == 2):
            self.swapping_process.ReGenerateMesh(self.model_part)

        # Update neigbours
        self.nodal_neighbour_search.Execute()

    def RefineOnSubscaleError(self,refine_tol,min_area,max_refinements):
        """ Refine elements using an error estimate based on the subscales.

            Equivalent to RefineErrorRatio, but specific for VMS2D and VMS3D
            elements and including an internal call to the error estimation
            process.

            refine_tol Maximum allowed value for the error estimate
            min_area Minimum allowed size (area or volume) for the new elements
            max_refinements Maximum refinement steps allowed from the same
            original element (to preserve mesh quality
        """

        self.refinement_utilities.SubscaleErrorEstimate(self.model_part)
        self.RefineOnErrorRatio(ERROR_RATIO,refine_tol,min_area,max_refinements)
        
