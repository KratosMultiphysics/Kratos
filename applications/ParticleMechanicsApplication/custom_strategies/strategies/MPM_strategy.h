//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


#if !defined(KRATOS_MPM_STRATEGY )
#define  KRATOS_MPM_STRATEGY

/* System includes */
#include <set>

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/kratos_flags.h"
#include "geometries/geometry.h"
#include "includes/element.h"
#include "boost/numeric/ublas/matrix.hpp"

// Application includes
#include "particle_mechanics_application_variables.h"

// Geometry utilities
#include "utilities/geometry_utilities.h"

// Custom includes
#include "custom_strategies/schemes/MPM_residual_based_bossak_scheme.hpp"
#include "custom_strategies/strategies/MPM_residual_based_newton_raphson_strategy.hpp"

// Core includes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/mpm_search_element_utility.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */
/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.
Detail class definition.

 */
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
         unsigned int TDim>
class MPMStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */
    //		typedef std::set<Dof::Pointer,ComparePDof> DofSetType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef typename TSparseSpace::DataType TDataType;
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TSparseSpace::MatrixPointerType TSystemMatrixPointerType;
    typedef typename TSparseSpace::VectorPointerType TSystemVectorPointerType;

    typedef Node < 3 > NodeType;
    typedef Geometry<NodeType> GeometryType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    typedef Scheme<TSparseSpace, TDenseSpace> TSchemeType;
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> TBuilderAndSolverType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > SolvingStrategyType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MPMStrategy);

    typedef typename ModelPart::DofType TDofType;
    typedef typename ModelPart::DofsArrayType DofsArrayType;

    typedef typename PointerVectorSet<TDofType, IndexedObject>::iterator DofIterator;
    typedef typename PointerVectorSet<TDofType, IndexedObject>::const_iterator DofConstantIterator;


    /**
     * @brief Default constructor of MPM Strategy.
     * @details
     * The grid model part contains all the information about ID, geometry and initial/boundary conditions
     * of the computational mesh.
     * In the costructor of time scheme the model part of material points is defined as:
     *
     * STEP 1:
     * The nodes, properties and process info of grid_model_part are assigned to the material points' mdpa.
     *
     * STEP 2:
     * loop over grid elements to evaluate:
     * - the MP integration weight and MP mass
     * - rGeo : connectivity of the grid element
     * - shape function values of the integration points of the grid element
     *
     * STEP 3:
     * loop over the integration points of a grid element to
     * - create a new MP element
     * - evaluate xg which is the coordinate of the integration point
     * - to assign all MP variables
     *
    */

    /*@{ */

    MPMStrategy(ModelPart& grid_model_part, ModelPart& initial_model_part, ModelPart& mpm_model_part, typename TLinearSolver::Pointer plinear_solver,
        Element const& rNewElement, std::string SolutionType = "static", int MaxIteration = 10, bool ComputeReaction = false, bool BlockBuilder = false,
        bool IsMixedFormulation = false, bool MoveMeshFlag = false)
        : SolvingStrategyType(grid_model_part, MoveMeshFlag), mr_grid_model_part(grid_model_part), mr_initial_model_part(initial_model_part),
        mr_mpm_model_part(mpm_model_part)
    {

        // Assigning the nodes to the new model part
        mr_mpm_model_part.Nodes() = mr_grid_model_part.Nodes();

        mr_mpm_model_part.SetProcessInfo(mr_grid_model_part.pGetProcessInfo());
        mr_mpm_model_part.SetBufferSize(mr_grid_model_part.GetBufferSize());
        mr_mpm_model_part.SetProperties(mr_initial_model_part.pProperties());

        // Prepare Dimension and Block Size
        unsigned int TBlock = TDim;
        if (IsMixedFormulation) TBlock ++;

        KRATOS_INFO("MPM_Strategy") << "Dimension Size = " << TDim << " and Block Size = " << TBlock << std::endl;

        // Create Material Point Element
        this->CreateMaterialPointElement(rNewElement, IsMixedFormulation);

        // Create Material Point Condition
        this->CreateMaterialPointCondition();

        // Define a standard static strategy to be used in the calculation
        if(SolutionType == "static" || SolutionType == "Static")
        {
            typename TSchemeType::Pointer pscheme = typename TSchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace >() );

            typename TBuilderAndSolverType::Pointer pBuilderAndSolver;
            if(BlockBuilder == true){
                KRATOS_INFO("MPM_Strategy") << "Block Builder is used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }
            else{
                KRATOS_INFO("MPM_Strategy") << "Block Builder is not used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }

            const double ratio_tolerance = 1e-04;
            const double always_converged_norm = 1e-09;
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));

            bool reform_DOF_at_each_iteration = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new MPMResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,MaxIteration,ComputeReaction,reform_DOF_at_each_iteration,MoveMeshFlag) );
        }

        // Define a quasi-static strategy to be used in the calculation
        else if(SolutionType == "quasi_static" || SolutionType == "Quasi-static")
        {
            double alpha_M;
            double dynamic;
            typename TSchemeType::Pointer pscheme = typename TSchemeType::Pointer( new MPMResidualBasedBossakScheme< TSparseSpace,TDenseSpace >(mr_grid_model_part, TDim, TBlock, alpha_M = 0.00, dynamic=0) );

            typename TBuilderAndSolverType::Pointer pBuilderAndSolver;
            if(BlockBuilder == true){
                KRATOS_INFO("MPM_Strategy") << "Block Builder is used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }
            else{
                KRATOS_INFO("MPM_Strategy") << "Block Builder is not used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }

            const double ratio_tolerance = 0.0001;
            const double always_converged_norm = 1e-09;
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));

            bool reform_DOF_at_each_iteration = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new MPMResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,MaxIteration,ComputeReaction,reform_DOF_at_each_iteration,MoveMeshFlag) );
        }

        // Define a dynamic strategy to be used in the calculation
        else if(SolutionType == "dynamic" || SolutionType == "Dynamic")
        {
            double alpha_M;
            double dynamic;
            typename TSchemeType::Pointer pscheme = typename TSchemeType::Pointer( new MPMResidualBasedBossakScheme< TSparseSpace,TDenseSpace >(mr_grid_model_part, TDim, TBlock, alpha_M = 0.0, dynamic=1) );

            typename TBuilderAndSolverType::Pointer pBuilderAndSolver;
            if(BlockBuilder == true){
                KRATOS_INFO("MPM_Strategy") << "Block Builder is used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }
            else{
                KRATOS_INFO("MPM_Strategy") << "Block Builder is not used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }

            const double ratio_tolerance = 0.00005;
            const double always_converged_norm = 1e-09;

            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));
            bool reform_DOF_at_each_iteration = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new MPMResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,MaxIteration,ComputeReaction,reform_DOF_at_each_iteration,MoveMeshFlag) );
        }

    }

    /*@} */

    /** Destructor.
     */

    /*@{ */
    virtual ~MPMStrategy()
    {
    }
    /*@} */


    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/

    /**
     * @brief This sets the level of echo for the solution strategy
     * @param Level of echo for the solution strategy
     * @details
     * {
     * 0->Mute... no echo at all
     * 1->Printing time and basic informations
     * 2->Printing linear solver data
     * 3->Print of debug informations: Echo of stiffness matrix, Dx, b...
     * }
     */
    void SetEchoLevel(const int Level) override
    {
        BaseType::mEchoLevel = Level;
        mp_solving_strategy->SetEchoLevel(Level);
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        // Initialize solving strategy: only to be done at the beginning of time step
        mp_solving_strategy->Initialize();
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     */
    void InitializeSolutionStep() override
    {
        // The nodal initial conditions are computed
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - InitializeSolutionStep" <<std::endl;
        mp_solving_strategy->InitializeSolutionStep();
    }

    /**
     * @brief Operation to predict the solution
     */
    void Predict() override
    {
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - Predict" <<std::endl;
        mp_solving_strategy->Predict();
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // Do solution iterations
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - SolveSolutionStep" <<std::endl;
        return mp_solving_strategy->SolveSolutionStep();
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     */
    void FinalizeSolutionStep() override
    {
        // The nodal solution are mapped from mesh to MP
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - FinalizeSolutionStep" <<std::endl;
        mp_solving_strategy->FinalizeSolutionStep();
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - Clear" <<std::endl;
        mp_solving_strategy->Clear();
    }

    /**
     * @brief Search element connectivity for each particle
     * @details A search is performed to know in which grid element the material point falls.
     * If one or more material points fall in the grid element, the grid element is
     * set to be active and its connectivity is associated to the material point
     * element.
     *
    */
    virtual void SearchElement(
        const std::size_t MaxNumberOfResults = 1000,
        const double Tolerance = 1.0e-5)
    {
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - Search Element" <<std::endl;
        MPMSearchElementUtility::SearchElement<TDim>(mr_grid_model_part, mr_mpm_model_part, MaxNumberOfResults, Tolerance);
    }


    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        const auto& mpm_process_info  = mr_mpm_model_part.GetProcessInfo();
        for (auto& r_element : mr_mpm_model_part.Elements())
            r_element.Check(mpm_process_info);

        const auto& grid_process_info = mr_grid_model_part.GetProcessInfo();
        for (auto& r_condition : mr_grid_model_part.Conditions())
            r_condition.Check(grid_process_info);

        return 0;
        KRATOS_CATCH("")
    }

    /**
     * @brief Function to Initiate material point element.
     * @details It is designed to be called ONCE by the class constructor.
     */
    virtual void CreateMaterialPointElement(Element const& rNewElement, bool IsMixedFormulation = false)
    {
        // Initialize zero the variables needed
        array_1d<double,3> xg = ZeroVector(3);
        array_1d<double,3> MP_Displacement = ZeroVector(3);
        array_1d<double,3> MP_Velocity = ZeroVector(3);
        array_1d<double,3> MP_Acceleration = ZeroVector(3);
        array_1d<double,3> MP_Volume_Acceleration = ZeroVector(3);

        Vector MP_Cauchy_Stress_Vector = ZeroVector(6);
        Vector MP_Almansi_Strain_Vector = ZeroVector(6);
        double MP_Pressure = 0.0;

        double MP_Mass;
        double MP_Volume;

        // Determine element index: This convention is done in order for the purpose of visualization in GiD
        const unsigned int number_elements = mr_grid_model_part.NumberOfElements() + mr_initial_model_part.NumberOfElements();
        const unsigned int number_nodes = mr_grid_model_part.NumberOfNodes();
        unsigned int last_element_id = (number_nodes > number_elements) ? (number_nodes + 1) : (number_elements+1);

        // Loop over the submodelpart of mr_initial_model_part
        for (ModelPart::SubModelPartIterator submodelpart_it = mr_initial_model_part.SubModelPartsBegin();
                submodelpart_it != mr_initial_model_part.SubModelPartsEnd(); submodelpart_it++)
        {
            ModelPart&  submodelpart = *submodelpart_it;
            std::string submodelpart_name = submodelpart.Name();

            mr_mpm_model_part.CreateSubModelPart(submodelpart_name);

            // Loop over the element of submodelpart and generate mpm element to be appended to the mr_mpm_model_part
            for (ModelPart::ElementIterator i = submodelpart.ElementsBegin();
                    i != submodelpart.ElementsEnd(); i++)
            {
                if(i->IsDefined(ACTIVE))
                {
                    Properties::Pointer properties = i->pGetProperties();
                    const int material_id = i->GetProperties().Id();
                    const double density  = i->GetProperties()[DENSITY];

                    // Check number of particles per element to be created
                    unsigned int particles_per_element;
                    if (i->GetProperties().Has( PARTICLES_PER_ELEMENT )){
                        particles_per_element = i->GetProperties()[PARTICLES_PER_ELEMENT];
                    }
                    else{
                        std::string warning_msg = "PARTICLES_PER_ELEMENT is not specified in Properties, ";
                        warning_msg += "1 Particle per element is assumed.";
                        KRATOS_WARNING("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                        particles_per_element = 1;
                    }

                    const Geometry< Node < 3 > >& rGeom = i->GetGeometry(); // current element's geometry
                    const GeometryData::KratosGeometryType rGeoType = rGeom.GetGeometryType();
                    Matrix shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                    if (rGeoType == GeometryData::Kratos_Tetrahedra3D4  || rGeoType == GeometryData::Kratos_Triangle2D3)
                    {
                        switch (particles_per_element)
                        {
                            case 1:
                                shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                                break;
                            case 3:
                                shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                                break;
                            case 6:
                                shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                                break;
                            case 12:
                                shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_5);
                                break;
                            case 16:
                                if (TDim==2){
                                    shape_functions_values = this->MP16ShapeFunctions();
                                    break;
                                }
                            case 33:
                                if (TDim==2) {
                                    shape_functions_values = this->MP33ShapeFunctions();
                                    break;
                                }
                            default:
                                std::string warning_msg = "The input number of PARTICLES_PER_ELEMENT: " + std::to_string(particles_per_element);
                                warning_msg += " is not available for Triangular" + std::to_string(TDim) + "D.\n";
                                warning_msg += "Available options are: 1, 3, 6, 12, 16 (only 2D), and 33 (only 2D).\n";
                                warning_msg += "The default number of particle: 3 is currently assumed.";
                                KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                break;
                        }
                    }
                    else if(rGeoType == GeometryData::Kratos_Hexahedra3D8  || rGeoType == GeometryData::Kratos_Quadrilateral2D4)
                    {
                        switch (particles_per_element)
                        {
                            case 1:
                                shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                                break;
                            case 4:
                                shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                                break;
                            case 9:
                                shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_3);
                                break;
                            case 16:
                                shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                                break;
                            default:
                                std::string warning_msg = "The input number of PARTICLES_PER_ELEMENT: " + std::to_string(particles_per_element);
                                warning_msg += " is not available for Quadrilateral" + std::to_string(TDim) + "D.\n";
                                warning_msg += "Available options are: 1, 4, 9, 16.\n";
                                warning_msg += "The default number of particle: 4 is currently assumed.";
                                KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                break;
                        }
                    }
                    else{
                        std::string error_msg = "The Geometry type of the Element given is invalid or currently not available. ";
                        error_msg += "Please remesh the problem domain to Triangle2D3N or Quadrilateral2D4N for 2D or ";
                        error_msg += "Tetrahedral3D4N or Hexahedral3D8N for 3D.";
                        KRATOS_ERROR << error_msg << std::endl;
                    }

                    // Number of MP per elements
                    const unsigned int integration_point_per_elements = shape_functions_values.size1();

                    // Evaluation of element area/volume
                    const double area = rGeom.Area();
                    if(TDim == 2 && i->GetProperties().Has( THICKNESS )){
						const double thickness = i->GetProperties()[THICKNESS];
						MP_Mass = area * thickness * density / integration_point_per_elements;
					}
					else {
                        MP_Mass = area * density / integration_point_per_elements;
                    }
                    MP_Volume = area / integration_point_per_elements;

                    // Loop over the material points that fall in each grid element
                    unsigned int new_element_id = 0;
                    for ( unsigned int PointNumber = 0; PointNumber < integration_point_per_elements; PointNumber++ )
                    {
                        // Create new material point element
                        new_element_id = last_element_id + PointNumber;
                        Element::Pointer p_element = rNewElement.Create(new_element_id, mr_grid_model_part.ElementsBegin()->GetGeometry(), properties);

                        const double MP_Density  = density;
                        const int MP_Material_Id = material_id;

                        xg.clear();

                        // Loop over the nodes of the grid element
                        for (unsigned int dim = 0; dim < rGeom.WorkingSpaceDimension(); dim++)
                        {
                            for ( unsigned int j = 0; j < rGeom.size(); j ++)
                            {
                                xg[dim] = xg[dim] + shape_functions_values(PointNumber, j) * rGeom[j].Coordinates()[dim];
                            }
                        }

                        // Setting particle element's initial condition
                        p_element->SetValue(MP_MATERIAL_ID, MP_Material_Id);
                        p_element->SetValue(MP_DENSITY, MP_Density);
                        p_element->SetValue(MP_MASS, MP_Mass);
                        p_element->SetValue(MP_VOLUME, MP_Volume);
                        p_element->SetValue(MP_COORD, xg);
                        p_element->SetValue(MP_DISPLACEMENT, MP_Displacement);
                        p_element->SetValue(MP_VELOCITY, MP_Velocity);
                        p_element->SetValue(MP_ACCELERATION, MP_Acceleration);
                        p_element->SetValue(MP_VOLUME_ACCELERATION, MP_Volume_Acceleration);
                        p_element->SetValue(MP_CAUCHY_STRESS_VECTOR, MP_Cauchy_Stress_Vector);
                        p_element->SetValue(MP_ALMANSI_STRAIN_VECTOR, MP_Almansi_Strain_Vector);

                        if(IsMixedFormulation)
                        {
                            p_element->SetValue(MP_PRESSURE, MP_Pressure);
                        }

                        // Add the MP Element to the model part
                        mr_mpm_model_part.GetSubModelPart(submodelpart_name).AddElement(p_element);
                    }

                    last_element_id += integration_point_per_elements;

                }
            }
        }
    }

    /**
     * @brief Function to Initiate material point condition.
     * @details It is designed to be called ONCE by the class constructor.
     */
    virtual void CreateMaterialPointCondition()
    {
        // Initialize zero the variables needed
        array_1d<double,3> mpc_xg = ZeroVector(3);
        array_1d<double,3> MPC_Normal = ZeroVector(3);
        array_1d<double,3> MPC_Displacement = ZeroVector(3);
        array_1d<double,3> MPC_Velocity = ZeroVector(3);
        array_1d<double,3> MPC_Acceleration = ZeroVector(3);

        double MPC_Area = 0.0;
        double MPC_Penalty_Factor = 0.0;

        // Determine condition index: This convention is done in order for the purpose of visualization in GiD
        const unsigned int number_conditions = mr_grid_model_part.NumberOfConditions();
        const unsigned int number_elements = mr_grid_model_part.NumberOfElements() + mr_initial_model_part.NumberOfElements() + mr_mpm_model_part.NumberOfElements();
        const unsigned int number_nodes = mr_grid_model_part.NumberOfNodes();
        unsigned int last_condition_id;
        if (number_elements > number_nodes && number_elements > number_conditions)
            last_condition_id = number_elements + 1;
        else if (number_nodes > number_elements && number_nodes > number_conditions)
            last_condition_id = number_nodes + 1;
        else
            last_condition_id = number_conditions + 1;

        // Loop over the submodelpart of mr_grid_model_part
        for (ModelPart::SubModelPartIterator submodelpart_it = mr_grid_model_part.SubModelPartsBegin();
                submodelpart_it != mr_grid_model_part.SubModelPartsEnd(); submodelpart_it++)
        {
            ModelPart& submodelpart = *submodelpart_it;

            // For submodelpart without condition, exit
            if (submodelpart.NumberOfConditions() != 0){

                std::string submodelpart_name = submodelpart.Name();
                mr_mpm_model_part.CreateSubModelPart(submodelpart_name);

                // For regular conditions: straight copy all conditions
                if (!submodelpart.ConditionsBegin()->Is(BOUNDARY)){
                    mr_mpm_model_part.SetConditions(submodelpart.pConditions());
                    mr_mpm_model_part.GetSubModelPart(submodelpart_name).SetConditions(submodelpart.pConditions());
                }
                // For boundary conditions: create particle conditions for all the necessary conditions
                else{

                    // NOTE: To create Particle Condition, we consider both the nodal position as well as the position of integration point
                    // Loop over the conditions of submodelpart and generate mpm condition to be appended to the mr_mpm_model_part
                    for (ModelPart::ConditionIterator i = submodelpart.ConditionsBegin();
                            i != submodelpart.ConditionsEnd(); i++)
                    {
                        Properties::Pointer properties = i->pGetProperties();
                        const int MPC_Condition_Id = i->GetProperties().Id();

                        // Flag whether condition is Neumann or Dirichlet
                        const bool is_neumann_condition = i->GetValue(MPC_IS_NEUMANN);

                        // Check number of particles per condition to be created
                        unsigned int particles_per_condition = 0; // Default zero
                        if (i->Has( PARTICLES_PER_CONDITION )){
                            particles_per_condition = i->GetValue(PARTICLES_PER_CONDITION);
                        }
                        else{
                            std::string warning_msg = "PARTICLES_PER_CONDITION is not specified, ";
                            warning_msg += "Only using nodal position is assumed: 1 (Point), 2 (Line), 3 (Triangular), 4 (Quadrilateral)";
                            KRATOS_WARNING("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                        }

                        // Get condition variables:
                        // Normal vector (normalized)
                        if (i->Has(NORMAL)) MPC_Normal = i->GetValue(NORMAL);
                        const double denominator = std::sqrt(MPC_Normal[0]*MPC_Normal[0] + MPC_Normal[1]*MPC_Normal[1] + MPC_Normal[2]*MPC_Normal[2]);
                        if (std::abs(denominator) > std::numeric_limits<double>::epsilon() ) MPC_Normal *= 1.0 / denominator;

                        // Get shape_function_values from defined particle_per_condition
                        auto& rGeom = i->GetGeometry(); // current condition's geometry
                        const GeometryData::KratosGeometryType rGeoType = rGeom.GetGeometryType();
                        Matrix shape_functions_values;
                        if (rGeoType == GeometryData::Kratos_Point2D  || rGeoType == GeometryData::Kratos_Point3D)
                        {
                            switch (particles_per_condition)
                            {
                                case 0: // Default case
                                    break;
                                case 1: // Only nodal
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Point" + std::to_string(TDim) + "D.\n";
                                    warning_msg += "Available option is: 1 (default).\n";
                                    warning_msg += "The default number of particle: 1 is currently assumed.";
                                    KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                    break;
                            }

                            if(is_neumann_condition)
                                KRATOS_ERROR << "Particle point load condition is not yet implemented." << std::endl;

                        }
                        else if (rGeoType == GeometryData::Kratos_Line2D2  || rGeoType == GeometryData::Kratos_Line3D2)
                        {
                            switch (particles_per_condition)
                            {
                                case 0: // Default case
                                    break;
                                case 2: // Only nodal
                                    break;
                                case 3:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                                    break;
                                case 4:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                                    break;
                                case 5:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_3);
                                    break;
                                case 6:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                                    break;
                                case 7:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_5);
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Line" + std::to_string(TDim) + "D.\n";
                                    warning_msg += "Available options are: 2 (default), 3, 4, 5, 6, 7.\n";
                                    warning_msg += "The default number of particle: 2 is currently assumed.";
                                    KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                    break;
                            }

                            if(is_neumann_condition)
                                KRATOS_ERROR << "Particle line load condition is not yet implemented." << std::endl;

                        }
                        else if(rGeoType == GeometryData::Kratos_Triangle3D3)
                        {
                            switch (particles_per_condition)
                            {
                                case 0: // Default case
                                    break;
                                case 3: // Only nodal
                                    break;
                                case 4:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                                    break;
                                case 6:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                                    break;
                                case 9:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                                    break;
                                case 15:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_5);
                                    break;
                                case 19:
                                    shape_functions_values = this->MP16ShapeFunctions();
                                    break;
                                case 36:
                                    shape_functions_values = this->MP33ShapeFunctions();
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Triangular" + std::to_string(TDim) + "D.\n";
                                    warning_msg += "Available options are: 3 (default), 4, 6, 9, 15, 19 and 36.\n";
                                    warning_msg += "The default number of particle: 3 is currently assumed.";
                                    KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                    break;
                            }

                            if(is_neumann_condition)
                                KRATOS_ERROR << "Particle surface load condition is not yet implemented." << std::endl;

                        }
                        else if(rGeoType == GeometryData::Kratos_Quadrilateral3D4)
                        {
                            switch (particles_per_condition)
                            {
                                case 0: // Default case
                                    break;
                                case 4: // Only nodal
                                    break;
                                case 5:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                                    break;
                                case 8:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                                    break;
                                case 13:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_3);
                                    break;
                                case 20:
                                    shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                                    break;
                                default:
                                    std::string warning_msg = "The input number of PARTICLES_PER_CONDITION: " + std::to_string(particles_per_condition);
                                    warning_msg += " is not available for Quadrilateral" + std::to_string(TDim) + "D.\n";
                                    warning_msg += "Available options are: 4 (default), 5, 8, 13, and 20.\n";
                                    warning_msg += "The default number of particle: 4 is currently assumed.";
                                    KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                    break;
                            }

                            if(is_neumann_condition)
                                KRATOS_ERROR << "Particle surface load condition is not yet implemented." << std::endl;

                        }
                        else{
                            std::string error_msg = "The Geometry type of the Condition given is invalid or currently not available. ";
                            error_msg += "Please remesh the problem domain to Point2D1N or Line2D2N for 2D or ";
                            error_msg += "Point3D1N, Line3D2N, Triangle3D3N or Quadrilateral3D4N for 3D.";
                            KRATOS_ERROR << error_msg << std::endl;
                        }

                        // Number of integration point per condition
                        const unsigned int integration_point_per_conditions = shape_functions_values.size1();

                        // Evaluation of geometric length/area
                        const double area = rGeom.Area();
                        MPC_Area = area / (rGeom.size() + integration_point_per_conditions);

                        // Check condition variables
                        if (i->Has(DISPLACEMENT))
                            MPC_Displacement = i->GetValue(DISPLACEMENT);
                        if (i->Has(VELOCITY))
                            MPC_Velocity = i->GetValue(VELOCITY);
                        if (i->Has(ACCELERATION))
                            MPC_Acceleration = i->GetValue(ACCELERATION);
                        if (i->Has(PENALTY_FACTOR))
                            MPC_Penalty_Factor = i->GetValue(PENALTY_FACTOR);
                        const bool is_slip = i->Is(SLIP);
                        const bool is_contact = i->Is(CONTACT);
                        const bool flip_normal_direction = i->Is(MODIFIED);

                        // If dirichlet boundary
                        std::string condition_type_name;
                        const GeometryData::KratosGeometryType rBackgroundGeoType = mr_grid_model_part.ElementsBegin()->GetGeometry().GetGeometryType();
                        if (!is_neumann_condition){
                            if (TDim==2){
                                if (rBackgroundGeoType == GeometryData::Kratos_Triangle2D3)
                                    condition_type_name = "MPMParticlePenaltyDirichletCondition2D3N";
                                else if (rBackgroundGeoType == GeometryData::Kratos_Quadrilateral2D4)
                                    condition_type_name = "MPMParticlePenaltyDirichletCondition2D4N";
                            }
                            else if (TDim==3){
                                if (rBackgroundGeoType == GeometryData::Kratos_Tetrahedra3D4)
                                    condition_type_name = "MPMParticlePenaltyDirichletCondition3D4N";
                                else if (rBackgroundGeoType == GeometryData::Kratos_Hexahedra3D8)
                                    condition_type_name = "MPMParticlePenaltyDirichletCondition3D8N";
                            }
                        }

                        // Get new condition
                        const Condition& new_condition = KratosComponents<Condition>::Get(condition_type_name);

                        // 1. Loop over the conditions to create inner particle condition
                        unsigned int new_condition_id = 0;
                        for ( unsigned int PointNumber = 0; PointNumber < integration_point_per_conditions; PointNumber++ )
                        {
                            // Create new material point condition
                            new_condition_id = last_condition_id + PointNumber;
                            Condition::Pointer p_condition = new_condition.Create(new_condition_id, mr_grid_model_part.ElementsBegin()->GetGeometry(), properties);

                            mpc_xg.clear();

                            // Loop over the nodes of the grid condition
                            for (unsigned int dim = 0; dim < rGeom.WorkingSpaceDimension(); dim++){
                                for ( unsigned int j = 0; j < rGeom.size(); j ++){
                                    mpc_xg[dim] = mpc_xg[dim] + shape_functions_values(PointNumber, j) * rGeom[j].Coordinates()[dim];
                                }
                            }

                            // Check Normal direction
                            if (flip_normal_direction) MPC_Normal *= -1.0;

                            // Setting particle condition's initial condition
                            // TODO: If any variable is added or remove here, please add and remove also at the second loop below
                            p_condition->SetValue(MPC_CONDITION_ID, MPC_Condition_Id);
                            p_condition->SetValue(MPC_COORD, mpc_xg);
                            p_condition->SetValue(MPC_AREA, MPC_Area);
                            p_condition->SetValue(MPC_NORMAL, MPC_Normal);
                            p_condition->SetValue(MPC_DISPLACEMENT, MPC_Displacement);
                            p_condition->SetValue(MPC_VELOCITY, MPC_Velocity);
                            p_condition->SetValue(MPC_ACCELERATION, MPC_Acceleration);
                            p_condition->SetValue(PENALTY_FACTOR, MPC_Penalty_Factor);
                            if (is_slip)
                                p_condition->Set(SLIP);
                            if (is_contact)
                                p_condition->Set(CONTACT);

                            // Add the MP Condition to the model part
                            mr_mpm_model_part.GetSubModelPart(submodelpart_name).AddCondition(p_condition);
                        }

                        last_condition_id += integration_point_per_conditions;

                        // 2. Loop over the nodes associated to each condition to create nodal particle condition
                        for ( unsigned int j = 0; j < rGeom.size(); j ++)
                        {
                            // Nodal normal vector is used
                            if (rGeom[j].Has(NORMAL)) MPC_Normal = rGeom[j].FastGetSolutionStepValue(NORMAL);
                            const double denominator = std::sqrt(MPC_Normal[0]*MPC_Normal[0] + MPC_Normal[1]*MPC_Normal[1] + MPC_Normal[2]*MPC_Normal[2]);
                            if (std::abs(denominator) > std::numeric_limits<double>::epsilon() ) MPC_Normal *= 1.0 / denominator;

                            // Check Normal direction
                            if (flip_normal_direction) MPC_Normal *= -1.0;

                            // Create new material point condition
                            new_condition_id = last_condition_id + j;
                            Condition::Pointer p_condition = new_condition.Create(new_condition_id, mr_grid_model_part.ElementsBegin()->GetGeometry(), properties);

                            mpc_xg.clear();
                            for (unsigned int dim = 0; dim < rGeom.WorkingSpaceDimension(); dim++){
                                mpc_xg[dim] = rGeom[j].Coordinates()[dim];
                            }

                            // Setting particle condition's initial condition
                            // TODO: If any variable is added or remove here, please add and remove also at the first loop above
                            p_condition->SetValue(MPC_CONDITION_ID, MPC_Condition_Id);
                            p_condition->SetValue(MPC_COORD, mpc_xg);
                            p_condition->SetValue(MPC_AREA, MPC_Area);
                            p_condition->SetValue(MPC_NORMAL, MPC_Normal);
                            p_condition->SetValue(MPC_DISPLACEMENT, MPC_Displacement);
                            p_condition->SetValue(MPC_VELOCITY, MPC_Velocity);
                            p_condition->SetValue(MPC_ACCELERATION, MPC_Acceleration);
                            p_condition->SetValue(PENALTY_FACTOR, MPC_Penalty_Factor);
                            if (is_slip)
                                p_condition->Set(SLIP);
                            if (is_contact)
                                p_condition->Set(CONTACT);

                            // Add the MP Condition to the model part
                            mr_mpm_model_part.GetSubModelPart(submodelpart_name).AddCondition(p_condition);
                        }

                        last_condition_id += rGeom.size();
                    }
                }
            }
        }
    }

    /**
     * @brief Function that return matrix of shape function value for 16 particles.
     * @details It is only possible to be used in 2D Triangular.
     */
    virtual Matrix MP16ShapeFunctions()
    {
        const double Na1 = 0.33333333333333;
        const double Nb1 = 0.45929258829272;
        const double Nb2 = 0.08141482341455;
        const double Nc1 = 0.17056930775176;
        const double Nc2 = 0.65886138449648;

        const double Nd1 = 0.05054722831703;
        const double Nd2 = 0.89890554336594;

        const double Ne1 = 0.26311282963464;
        const double Ne2 = 0.72849239295540;
        const double Ne3 = 0.00839477740996;

        BoundedMatrix<double,16,3> MP_ShapeFunctions;

        MP_ShapeFunctions(0,0) = Na1;
        MP_ShapeFunctions(0,1) = Na1;
        MP_ShapeFunctions(0,2) = Na1;

        MP_ShapeFunctions(1,0) = Nb1;
        MP_ShapeFunctions(1,1) = Nb1;
        MP_ShapeFunctions(1,2) = Nb2;

        MP_ShapeFunctions(2,0) = Nb1;
        MP_ShapeFunctions(2,1) = Nb2;
        MP_ShapeFunctions(2,2) = Nb1;

        MP_ShapeFunctions(3,0) = Nb2;
        MP_ShapeFunctions(3,1) = Nb1;
        MP_ShapeFunctions(3,2) = Nb1;

        MP_ShapeFunctions(4,0) = Nc1;
        MP_ShapeFunctions(4,1) = Nc1;
        MP_ShapeFunctions(4,2) = Nc2;

        MP_ShapeFunctions(5,0) = Nc1;
        MP_ShapeFunctions(5,1) = Nc2;
        MP_ShapeFunctions(5,2) = Nc1;

        MP_ShapeFunctions(6,0) = Nc2;
        MP_ShapeFunctions(6,1) = Nc1;
        MP_ShapeFunctions(6,2) = Nc1;

        MP_ShapeFunctions(7,0) = Nd1;
        MP_ShapeFunctions(7,1) = Nd1;
        MP_ShapeFunctions(7,2) = Nd2;

        MP_ShapeFunctions(8,0) = Nd1;
        MP_ShapeFunctions(8,1) = Nd2;
        MP_ShapeFunctions(8,2) = Nd1;

        MP_ShapeFunctions(9,0) = Nd2;
        MP_ShapeFunctions(9,1) = Nd1;
        MP_ShapeFunctions(9,2) = Nd1;

        MP_ShapeFunctions(10,0) = Ne1;
        MP_ShapeFunctions(10,1) = Ne2;
        MP_ShapeFunctions(10,2) = Ne3;

        MP_ShapeFunctions(11,0) = Ne2;
        MP_ShapeFunctions(11,1) = Ne3;
        MP_ShapeFunctions(11,2) = Ne1;

        MP_ShapeFunctions(12,0) = Ne3;
        MP_ShapeFunctions(12,1) = Ne1;
        MP_ShapeFunctions(12,2) = Ne2;

        MP_ShapeFunctions(13,0) = Ne2;
        MP_ShapeFunctions(13,1) = Ne1;
        MP_ShapeFunctions(13,2) = Ne3;

        MP_ShapeFunctions(14,0) = Ne1;
        MP_ShapeFunctions(14,1) = Ne3;
        MP_ShapeFunctions(14,2) = Ne2;

        MP_ShapeFunctions(15,0) = Ne3;
        MP_ShapeFunctions(15,1) = Ne2;
        MP_ShapeFunctions(15,2) = Ne1;

        //MP_ShapeFunctions = [(Na1, Na1, Na1),(Nb1, Nb1, Nb2),(Nb1, Nb2, Nb1),(Nb2, Nb1, Nb1),
        //                    (Nc1, Nc1, Nc2),(Nc1, Nc2, Nc1),(Nc2, Nc1, Nc1),(Nd1, Nd1, Nd2),
        //                    (Nd1, Nd2, Nd1),(Nd2, Nd1, Nd1),(Ne1, Ne2, Ne3),(Ne2, Ne3, Ne1),
        //                    (Ne3, Ne1, Ne2),(Ne2, Ne1, Ne3),(Ne1, Ne3, Ne2),(Ne3, Ne2, Ne1)];

        return MP_ShapeFunctions;
    }

    /**
     * @brief Function that return matrix of shape function value for 33 particles.
     * @details It is only possible to be used in 2D Triangular.
     */
    virtual Matrix MP33ShapeFunctions()
    {
        const double Na2 = 0.02356522045239;
        const double Na1 = 0.488217389773805;

        const double Nb2 = 0.120551215411079;
        const double Nb1 = 0.43972439229446;

        const double Nc2 = 0.457579229975768;
        const double Nc1 = 0.271210385012116;

        const double Nd2 = 0.744847708916828;
        const double Nd1 = 0.127576145541586;

        const double Ne2 = 0.957365299093579;
        const double Ne1 = 0.021317350453210;

        const double Nf1 = 0.115343494534698;
        const double Nf2 = 0.275713269685514;
        const double Nf3 = 0.608943235779788;

        const double Ng1 = 0.022838332222257;
        const double Ng2 = 0.281325580989940;
        const double Ng3 = 0.695836086787803;

        const double Nh1 = 0.025734050548330;
        const double Nh2 = 0.116251915907597;
        const double Nh3 = 0.858014033544073;

        BoundedMatrix<double,33,3> MP_ShapeFunctions;

        MP_ShapeFunctions(0,0) = Na1;
        MP_ShapeFunctions(0,1) = Na1;
        MP_ShapeFunctions(0,2) = Na2;

        MP_ShapeFunctions(1,0) = Na1;
        MP_ShapeFunctions(1,1) = Na2;
        MP_ShapeFunctions(1,2) = Na1;

        MP_ShapeFunctions(2,0) = Na2;
        MP_ShapeFunctions(2,1) = Na1;
        MP_ShapeFunctions(2,2) = Na1;


        MP_ShapeFunctions(3,0) = Nb1;
        MP_ShapeFunctions(3,1) = Nb1;
        MP_ShapeFunctions(3,2) = Nb2;

        MP_ShapeFunctions(4,0) = Nb1;
        MP_ShapeFunctions(4,1) = Nb2;
        MP_ShapeFunctions(4,2) = Nb1;

        MP_ShapeFunctions(5,0) = Nb2;
        MP_ShapeFunctions(5,1) = Nb1;
        MP_ShapeFunctions(5,2) = Nb1;

        MP_ShapeFunctions(6,0) = Nc1;
        MP_ShapeFunctions(6,1) = Nc1;
        MP_ShapeFunctions(6,2) = Nc2;

        MP_ShapeFunctions(7,0) = Nc1;
        MP_ShapeFunctions(7,1) = Nc2;
        MP_ShapeFunctions(7,2) = Nc1;

        MP_ShapeFunctions(8,0) = Nc2;
        MP_ShapeFunctions(8,1) = Nc1;
        MP_ShapeFunctions(8,2) = Nc1;

        MP_ShapeFunctions(9,0) = Nd1;
        MP_ShapeFunctions(9,1) = Nd1;
        MP_ShapeFunctions(9,2) = Nd2;

        MP_ShapeFunctions(10,0) = Nd1;
        MP_ShapeFunctions(10,1) = Nd2;
        MP_ShapeFunctions(10,2) = Nd1;

        MP_ShapeFunctions(11,0) = Nd2;
        MP_ShapeFunctions(11,1) = Nd1;
        MP_ShapeFunctions(11,2) = Nd1;

        MP_ShapeFunctions(12,0) = Ne1;
        MP_ShapeFunctions(12,1) = Ne1;
        MP_ShapeFunctions(12,2) = Ne2;

        MP_ShapeFunctions(13,0) = Ne1;
        MP_ShapeFunctions(13,1) = Ne2;
        MP_ShapeFunctions(13,2) = Ne1;

        MP_ShapeFunctions(14,0) = Ne2;
        MP_ShapeFunctions(14,1) = Ne1;
        MP_ShapeFunctions(14,2) = Ne1;

        MP_ShapeFunctions(15,0) = Nf1;
        MP_ShapeFunctions(15,1) = Nf2;
        MP_ShapeFunctions(15,2) = Nf3;

        MP_ShapeFunctions(16,0) = Nf2;
        MP_ShapeFunctions(16,1) = Nf3;
        MP_ShapeFunctions(16,2) = Nf1;

        MP_ShapeFunctions(17,0) = Nf3;
        MP_ShapeFunctions(17,1) = Nf1;
        MP_ShapeFunctions(17,2) = Nf2;

        MP_ShapeFunctions(18,0) = Nf2;
        MP_ShapeFunctions(18,1) = Nf1;
        MP_ShapeFunctions(18,2) = Nf3;

        MP_ShapeFunctions(19,0) = Nf1;
        MP_ShapeFunctions(19,1) = Nf3;
        MP_ShapeFunctions(19,2) = Nf2;

        MP_ShapeFunctions(20,0) = Nf3;
        MP_ShapeFunctions(20,1) = Nf2;
        MP_ShapeFunctions(20,2) = Nf1;

        MP_ShapeFunctions(21,0) = Ng1;
        MP_ShapeFunctions(21,1) = Ng2;
        MP_ShapeFunctions(21,2) = Ng3;

        MP_ShapeFunctions(22,0) = Ng2;
        MP_ShapeFunctions(22,1) = Ng3;
        MP_ShapeFunctions(22,2) = Ng1;

        MP_ShapeFunctions(23,0) = Ng3;
        MP_ShapeFunctions(23,1) = Ng1;
        MP_ShapeFunctions(23,2) = Ng2;

        MP_ShapeFunctions(24,0) = Ng2;
        MP_ShapeFunctions(24,1) = Ng1;
        MP_ShapeFunctions(24,2) = Ng3;

        MP_ShapeFunctions(25,0) = Ng1;
        MP_ShapeFunctions(25,1) = Ng3;
        MP_ShapeFunctions(25,2) = Ng2;

        MP_ShapeFunctions(26,0) = Ng3;
        MP_ShapeFunctions(26,1) = Ng2;
        MP_ShapeFunctions(26,2) = Ng1;

        MP_ShapeFunctions(27,0) = Nh1;
        MP_ShapeFunctions(27,1) = Nh2;
        MP_ShapeFunctions(27,2) = Nh3;

        MP_ShapeFunctions(28,0) = Nh2;
        MP_ShapeFunctions(28,1) = Nh3;
        MP_ShapeFunctions(28,2) = Nh1;

        MP_ShapeFunctions(29,0) = Nh3;
        MP_ShapeFunctions(29,1) = Nh1;
        MP_ShapeFunctions(29,2) = Nh2;

        MP_ShapeFunctions(30,0) = Nh2;
        MP_ShapeFunctions(30,1) = Nh1;
        MP_ShapeFunctions(30,2) = Nh3;

        MP_ShapeFunctions(31,0) = Nh1;
        MP_ShapeFunctions(31,1) = Nh3;
        MP_ShapeFunctions(31,2) = Nh2;

        MP_ShapeFunctions(32,0) = Nh3;
        MP_ShapeFunctions(32,1) = Nh2;
        MP_ShapeFunctions(32,2) = Nh1;

        return MP_ShapeFunctions;
    }

protected:
    /**@name Protected static Member Variables */
    /*@{ */

    ModelPart& mr_grid_model_part;
    ModelPart& mr_initial_model_part;
    ModelPart& mr_mpm_model_part;

    SolvingStrategyType::Pointer mp_solving_strategy;

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */


private:

    /*@} */
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /** Copy constructor.
     */

    /*@} */

}; /* Class NewSolvingStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_MPM_STRATEGY  defined */

