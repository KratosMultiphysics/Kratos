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

// Application includes
#include "particle_mechanics_application.h"

// Geometry utilities
#include "utilities/geometry_utilities.h"

// Custom includes
#include "custom_strategies/schemes/MPM_residual_based_bossak_scheme.hpp"
#include "custom_strategies/strategies/MPM_residual_based_newton_raphson_strategy.hpp"
#include "custom_elements/updated_lagrangian.hpp"

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
#include "utilities/binbased_fast_point_locator.h"

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

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
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


    /** Constructor.
     * The grid model part contains all the information about ID, geometry and initial/boundary conditions
     * of the computational mesh.
     * In the costructor of time scheme the model part of material points is defined as:
     * STEP 1:
     * The nodes, properties and process info of grid_model_part are assigned to the material points' mdpa.
     *
     * STEP 2:
     *loop over grid elements to evaluate:
     * - the MP integration weight and MP mass
     * - rGeo : connectivity of the grid element
     * - shape function values of the integration points of the grid element
     *
     * STEP 3:
     *loop over the integration points of a grid element to
     * - create a new MP element
     * - evaluate xg which is the coordinate of the integration point
     * - to assign all MP variables
     *
    */

    /*@{ */

    MPMStrategy(ModelPart& grid_model_part, ModelPart& initial_model_part, ModelPart& mpm_model_part, typename TLinearSolver::Pointer plinear_solver,
        Element const& NewElement, bool MoveMeshFlag = false, std::string SolutionType = "StaticType", std::string GeometryElement = "Triangle",
        int NumPar = 3, bool BlockBuilder = false, bool isMixedFormulation = false)
        : SolvingStrategyType(grid_model_part, MoveMeshFlag), mr_grid_model_part(grid_model_part), mr_initial_model_part(initial_model_part),
        mr_mpm_model_part(mpm_model_part), m_GeometryElement(GeometryElement), m_NumPar(NumPar)
    {

        // Assigning the nodes to the new model part
        mpm_model_part.Nodes() = grid_model_part.Nodes();

        mpm_model_part.SetProcessInfo(grid_model_part.pGetProcessInfo());
        mpm_model_part.SetBufferSize(grid_model_part.GetBufferSize());
        mpm_model_part.SetProperties(initial_model_part.pProperties());
        mpm_model_part.SetConditions(grid_model_part.pConditions());

        array_1d<double,3> xg = ZeroVector(3);
        array_1d<double,3> MP_Displacement = ZeroVector(3);
        array_1d<double,3> MP_Velocity = ZeroVector(3);
        array_1d<double,3> MP_Acceleration = ZeroVector(3);
        array_1d<double,3> Aux_MP_Velocity = ZeroVector(3);
        array_1d<double,3> Aux_MP_Acceleration = ZeroVector(3);
        array_1d<double,3> MP_Volume_Acceleration = ZeroVector(3);

        Vector MP_Cauchy_Stress_Vector = ZeroVector(6);
        Vector MP_Almansi_Strain_Vector = ZeroVector(6);
        double MP_Pressure = 0.0;
        double Aux_MP_Pressure = 0.0;

        double MP_Mass;
        double MP_Volume;

        // Prepare Dimension and Block Size
        unsigned int TBlock = TDim;
        if (isMixedFormulation) TBlock ++;

        KRATOS_INFO("MPM_Strategy") << "Dimension Size = " << TDim << " and Block Size = " << TBlock << std::endl;

        unsigned int k = 0;
        const unsigned int number_elements = grid_model_part.NumberOfElements();
        const unsigned int number_nodes = grid_model_part.NumberOfNodes();
        int new_element_id = 0;

        // Loop over the submodelpart of initial_model_part
        for (ModelPart::SubModelPartIterator submodelpart_it = initial_model_part.SubModelPartsBegin();
                submodelpart_it != initial_model_part.SubModelPartsEnd(); submodelpart_it++)
        {
            ModelPart& submodelpart = *submodelpart_it;
            std::string submodelpart_name = submodelpart.Name();

            mpm_model_part.CreateSubModelPart(submodelpart_name);

            // Loop over the element of submodelpart's submodelpart and generate mpm element to be appended to the mpm_model_part
            for (ModelPart::ElementIterator i = submodelpart.ElementsBegin();
                    i != submodelpart.ElementsEnd(); i++)
            {
                if(i->IsDefined(ACTIVE))
                {
                    Properties::Pointer properties = i->pGetProperties();
                    int material_id = i->GetProperties().Id();
                    double density  = i->GetProperties()[DENSITY];

                    Geometry< Node < 3 > >& rGeom = i->GetGeometry(); // current element's connectivity
                    Matrix shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                    if (m_GeometryElement == "Triangle")
                    {
                        switch (m_NumPar)
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
                                std::string warning_msg = "The input number of particle: " + std::to_string(m_NumPar);
                                warning_msg += " is not available for Triangular" + std::to_string(TDim) + "D.\n";
                                warning_msg += "Available options are: 1, 3, 6, 12, 16 (only 2D), and 33 (only 2D).\n";
                                warning_msg += "The default number of particle: 3 is currently assumed.";
                                KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                break;
                        }
                    }
                    else if(m_GeometryElement == "Quadrilateral")
                    {
                        switch (m_NumPar)
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
                                std::string warning_msg = "The input number of particle: " + std::to_string(m_NumPar);
                                warning_msg += " is not available for Quadrilateral" + std::to_string(TDim) + "D.\n";
                                warning_msg += "Available options are: 1, 4, 9, 16.\n";
                                warning_msg += "The default number of particle: 4 is currently assumed.";
                                KRATOS_INFO("MPM_Strategy") << "WARNING: " << warning_msg << std::endl;
                                break;
                        }
                    }

                    // Number of MP per elements
                    const unsigned int integration_point_per_elements = shape_functions_values.size1();

                    // Evaluation of element area/volume
                    const double area = rGeom.Area();

                    MP_Mass   = area * density / integration_point_per_elements;
                    MP_Volume = area / integration_point_per_elements;

                    // Loop over the material points that fall in each grid element
                    for ( unsigned int PointNumber = 0; PointNumber < integration_point_per_elements; PointNumber++ )
                    {
                        if(number_elements > number_nodes)
                        {
                            new_element_id = (1+PointNumber+number_elements)+(integration_point_per_elements*k);
                        }
                        else
                        {
                            new_element_id = (1+PointNumber+number_nodes)+(integration_point_per_elements*k);
                        }
                        Element::Pointer p_element = NewElement.Create(new_element_id, rGeom, properties);
                        double MP_Density  = density;
                        int MP_Material_Id = material_id;

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
                        p_element->SetValue(MP_NUMBER, integration_point_per_elements);
                        p_element->SetValue(MP_MATERIAL_ID, MP_Material_Id);
                        p_element->SetValue(MP_DENSITY, MP_Density);
                        p_element->SetValue(MP_MASS, MP_Mass);
                        p_element->SetValue(MP_VOLUME, MP_Volume);
                        p_element->SetValue(GAUSS_COORD, xg);
                        p_element->SetValue(MP_DISPLACEMENT, MP_Displacement);
                        p_element->SetValue(MP_VELOCITY, MP_Velocity);
                        p_element->SetValue(MP_ACCELERATION, MP_Acceleration);
                        p_element->SetValue(AUX_MP_VELOCITY, Aux_MP_Velocity);
                        p_element->SetValue(AUX_MP_ACCELERATION, Aux_MP_Acceleration);
                        p_element->SetValue(MP_VOLUME_ACCELERATION, MP_Volume_Acceleration);
                        p_element->SetValue(MP_CAUCHY_STRESS_VECTOR, MP_Cauchy_Stress_Vector);
                        p_element->SetValue(MP_ALMANSI_STRAIN_VECTOR, MP_Almansi_Strain_Vector);

                        if(isMixedFormulation)
                        {
                            p_element->SetValue(MP_PRESSURE, MP_Pressure);
                            p_element->SetValue(AUX_MP_PRESSURE, Aux_MP_Pressure);
                        }

                        // Add the MP Element to the model part
                        mpm_model_part.GetSubModelPart(submodelpart_name).AddElement(p_element);
                    }

                    k +=1;

                }


            }

        }

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

            int max_iteration = 20;
            bool calculate_reaction = false;
            bool reform_DOF_at_each_iteration = false;
            bool move_mesh_flag = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new MPMResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,max_iteration,calculate_reaction,reform_DOF_at_each_iteration,move_mesh_flag) );
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

            int max_iteration = 100;
            bool calculate_reaction = false;
            bool reform_DOF_at_each_iteration = false;
            bool move_mesh_flag = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new MPMResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,max_iteration,calculate_reaction,reform_DOF_at_each_iteration,move_mesh_flag) );
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
            int max_iteration = 20;
            bool calculate_reaction = false;
            bool reform_DOF_at_each_iteration = false;
            bool move_mesh_flag = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new MPMResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,max_iteration,calculate_reaction,reform_DOF_at_each_iteration,move_mesh_flag) );
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

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/
    /*@{ */

    /**
    operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
    void Predict() override
    {
    }

    /**
    Initialization of member variables and prior operations
     */
    void Initialize() override
    {
    }

    /**
    the problem of interest is solved
     */
    double Solve() override
    {
        // Check which nodes and elements are ACTIVE and populate the MPM model part
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - Start" <<std::endl;
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Search Element - Start" <<std::endl;
        this->SearchElement(mr_grid_model_part, mr_mpm_model_part);

        // Only perform this once
        mp_solving_strategy->Initialize();

        // The nodal initial conditions are computed
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - InitializeSolutionStep" <<std::endl;
        mp_solving_strategy->InitializeSolutionStep();

        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - Predict" <<std::endl;
        mp_solving_strategy->Predict();

        // Do solution iterations
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - SolveSolutionStep" <<std::endl;
        mp_solving_strategy->SolveSolutionStep();

        // The nodal solution are mapped from mesh to MP
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - FinalizeSolutionStep" <<std::endl;
        mp_solving_strategy->FinalizeSolutionStep();

        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - Clear" <<std::endl;
        mp_solving_strategy->Clear();

        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - End" <<std::endl;

		return 0.00;
    }

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

    /** SearchElement.
     * A search is performed to know in which grid element the material point falls.
     *
     * If one or more material points fall in the grid element, the grid element is
     * set to be active and its connectivity is associated to the material point
     * element.
     *
     * STEPS:
     * 1) All the elements are set to be INACTIVE
     * 2) A searching is performed and the grid elements which contain at least a MP are set to be ACTIVE
     *
    */
    virtual void SearchElement(
        ModelPart& grid_model_part,
        ModelPart& mpm_model_part,
        const std::size_t MaxNumberOfResults = 1000,
        const double Tolerance = 1.0e-5)
    {
        // Reset elements to inactive
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(grid_model_part.Elements().size()); ++i){

			auto element_itr = grid_model_part.Elements().begin() + i;
			element_itr->Reset(ACTIVE);
            if (m_GeometryElement == "Triangle"){
                element_itr->GetGeometry()[0].Reset(ACTIVE);
                element_itr->GetGeometry()[1].Reset(ACTIVE);
                element_itr->GetGeometry()[2].Reset(ACTIVE);

                if (TDim ==3)
                {
                    element_itr->GetGeometry()[3].Reset(ACTIVE);
                }
            }
            else if (m_GeometryElement == "Quadrilateral"){
                element_itr->GetGeometry()[0].Reset(ACTIVE);
                element_itr->GetGeometry()[1].Reset(ACTIVE);
                element_itr->GetGeometry()[2].Reset(ACTIVE);
                element_itr->GetGeometry()[3].Reset(ACTIVE);

                if (TDim ==3)
                {
                    element_itr->GetGeometry()[4].Reset(ACTIVE);
                    element_itr->GetGeometry()[5].Reset(ACTIVE);
                    element_itr->GetGeometry()[6].Reset(ACTIVE);
                    element_itr->GetGeometry()[7].Reset(ACTIVE);
                }
            }
		}

        // Search background grid and make element active
        Vector N;
        const int max_result = 1000;

        #pragma omp parallel
        {
            BinBasedFastPointLocator<TDim> SearchStructure(grid_model_part);
            SearchStructure.UpdateSearchDatabase();

            typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_result);

            #pragma omp for
            for(int i = 0; i < static_cast<int>(mpm_model_part.Elements().size()); ++i){

                auto element_itr = mpm_model_part.Elements().begin() + i;

                array_1d<double,3> xg = element_itr->GetValue(GAUSS_COORD);
                typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                Element::Pointer pelem;

                // FindPointOnMesh find the element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg, N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                if (is_found == true)
                {
                    pelem->Set(ACTIVE);
                    element_itr->GetGeometry() = pelem->GetGeometry();
                    if (m_GeometryElement == "Triangle")
                    {
                        element_itr->GetGeometry()[0].Set(ACTIVE);
                        element_itr->GetGeometry()[1].Set(ACTIVE);
                        element_itr->GetGeometry()[2].Set(ACTIVE);
                        if (TDim ==3)
                        {
                            element_itr->GetGeometry()[3].Set(ACTIVE);
                        }
                    }
                    else if(m_GeometryElement == "Quadrilateral")
                    {
                        element_itr->GetGeometry()[0].Set(ACTIVE);
                        element_itr->GetGeometry()[1].Set(ACTIVE);
                        element_itr->GetGeometry()[2].Set(ACTIVE);
                        element_itr->GetGeometry()[3].Set(ACTIVE);
                        if (TDim ==3)
                        {
                            element_itr->GetGeometry()[4].Set(ACTIVE);
                            element_itr->GetGeometry()[5].Set(ACTIVE);
                            element_itr->GetGeometry()[6].Set(ACTIVE);
                            element_itr->GetGeometry()[7].Set(ACTIVE);
                        }
                    }
                }
                else{
                    KRATOS_INFO("MPM_Strategy.SearchElement") << "WARNING: Search Element for Particle " << element_itr->Id()
                        << " is failed. Geometry is cleared." << std::endl;

                    element_itr->GetGeometry().clear();
                    element_itr->Reset(ACTIVE);
                    element_itr->Set(TO_ERASE);
                }
            }
        }
    }


    /**
     * Function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        for (ModelPart::ElementsContainerType::iterator it = mr_mpm_model_part.ElementsBegin();
                it != mr_mpm_model_part.ElementsEnd(); it++)
        {
            it->Check(mr_mpm_model_part.GetProcessInfo());
        }

        for (ModelPart::ConditionsContainerType::iterator it = mr_grid_model_part.ConditionsBegin();
                it != mr_grid_model_part.ConditionsEnd(); it++)
        {
            it->Check(mr_grid_model_part.GetProcessInfo());
        }
        return 0;
        KRATOS_CATCH("")
    }

    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */

    //level of echo for the solving strategy
    //int mEchoLevel;

    //settings for the rebuilding of the stiffness matrix
    //int mRebuildLevel;
    //bool mStiffnessMatrixIsBuilt;

    ModelPart& mr_grid_model_part;
    ModelPart& mr_initial_model_part;
    ModelPart& mr_mpm_model_part;
    std::string m_GeometryElement;
    int m_NumPar;

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

