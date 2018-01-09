/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//   Project Name:        KratosParticleMechanicsApplication $
//   Last modified by:    $Author:                    ilaria $
//   Date:                $Date:                   July 2015 $
//   Revision:            $Revision:                     0.0 $


#if !defined(KRATOS_MPM_STRATEGY )
#define  KRATOS_MPM_STRATEGY


/* System includes */
#include <set>
//#include <chrono>
/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "spaces/ublas_space.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/kratos_flags.h"
#include "geometries/geometry.h"
#include "includes/element.h"
#include "solid_mechanics_application.h"
//geometry utilities
#include "utilities/geometry_utilities.h"

#include "particle_mechanics_application.h"

#include "custom_elements/updated_lagrangian.hpp"


#include "custom_strategies/schemes/MPM_residual_based_bossak_scheme.hpp"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "custom_strategies/convergence_criteria/displacement_convergence_criterion.hpp"

#include "custom_strategies/strategies/MPM_residual_based_newton_raphson_strategy.hpp"
//#include "custom_strategies/strategies/MPM_strategy.h"

#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "linear_solvers/linear_solver.h"

#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

#include "utilities/binbased_fast_point_locator.h"
#include "custom_utilities/quad_binbased_fast_point_locator.h"


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

  \URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

        \URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

          \URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

                \URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


                        \URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

                          \URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

                                \URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

                                  \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


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

    /** Counted pointer of ClassName */
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

    MPMStrategy(ModelPart& grid_model_part, ModelPart& initial_model_part, ModelPart& mpm_model_part, typename TLinearSolver::Pointer plinear_solver,Element const& NewElement, bool MoveMeshFlag = false, std::string SolutionType = "StaticType", std::string GeometryElement = "Triangle", int NumPar = 3)
        : SolvingStrategyType(grid_model_part, MoveMeshFlag), mr_grid_model_part(grid_model_part), mr_initial_model_part(initial_model_part), mr_mpm_model_part(mpm_model_part), m_GeometryElement(GeometryElement), m_NumPar(NumPar)
    {

        //populate for the first time the mpm_model_part

        //assigning the nodes to the new model part
        mpm_model_part.Nodes() = grid_model_part.Nodes();

        mpm_model_part.SetProcessInfo(grid_model_part.pGetProcessInfo());
        mpm_model_part.SetBufferSize(grid_model_part.GetBufferSize());
        mpm_model_part.SetProperties(initial_model_part.pProperties());
        mpm_model_part.SetConditions(grid_model_part.pConditions());




        array_1d<double,3> xg = ZeroVector(3);
        array_1d<double,3> MP_Displacement = ZeroVector(3);
        array_1d<double,3> MP_Velocity = ZeroVector(3);
        //double MP_KineticEnergy = 0.0;
        //double MP_StrainEnergy = 0.0;
        //Vector MP_CauchyVector = ZeroVector(3);
        //Vector MP_AlmansiVector = ZeroVector(3);
        //Matrix MP_ConstitutiveMatrix = ZeroMatrix(6,6);
        double MP_Mass;
        double MP_Volume;

        unsigned int k = 0;
        const unsigned int number_elements = grid_model_part.NumberOfElements();
        const unsigned int number_nodes = grid_model_part.NumberOfNodes();
        int new_element_id = 0;

        for (ModelPart::ElementIterator i = initial_model_part.ElementsBegin();
                i != initial_model_part.ElementsEnd(); i++)
        {
            if(i->IsDefined(ACTIVE))
            {


                Properties::Pointer properties = i->pGetProperties();
                double Density = i->GetProperties()[DENSITY];
                //std::cout<< "Density "<< Density<<std::endl;
                Geometry< Node < 3 > >& rGeom = i->GetGeometry(); // current element's connectivity
                Matrix shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                if (m_GeometryElement == "Triangle")
                {
                    if(m_NumPar == 1)
                    {
                        shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                    }
                    else if(m_NumPar == 3)
                    {
                        shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                    }
                    else if(m_NumPar == 6)
                    {
                        shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                    }
                    else if(m_NumPar == 16)
                    {
                        shape_functions_values = this->MP16ShapeFunctions();
                    }
                }
                else if(m_GeometryElement == "Quadrilateral")
                {
                    if(m_NumPar == 1)
                    {
                        shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                    }
                    else if(m_NumPar == 4)
                    {
                        shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2);
                    }
                    else if(m_NumPar == 9)
                    {
                        shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_3);
                    }
                    else if(m_NumPar == 16)
                    {
                        shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_4);
                    }
                }
                //Matrix shape_functions_values = this->MP16ShapeFunctions();
                //std::cout<<"shape_functions_values "<< shape_functions_values<<std::endl;
                //const GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints( GeometryData::GI_GAUSS_4);

                //INITIAL NUMBER OF MATERIAL POINTS PER ELEMENT
                unsigned int integration_point_per_elements = shape_functions_values.size1();


                //evaluation of element area/volume
                //double area = GeometryUtils::CalculateVolume2D(rGeom);

                double area = rGeom.Area();
                //int integration_point_per_element = 3;

                MP_Mass = area * Density / integration_point_per_elements;
                //std::cout<<"MP_Mass "<<MP_Mass<<std::endl;
                MP_Volume = area / integration_point_per_elements;


                //loop over the material points that fall in each grid element
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
                    //Element::Pointer p_element = NewElement.Create((1+PointNumber+number_nodes)+(integration_point_per_element*k), rGeom, properties);
                    double MP_Density = Density;
                    //std::cout<<"MPM ID "<<(1+PointNumber+number_elements)+(integration_point_per_element*k)<<std::endl;
                    xg.clear();


                    //loop over the nodes of the grid element
                    for (unsigned int dim = 0; dim < rGeom.WorkingSpaceDimension(); dim++)
                    {
                        for ( unsigned int j = 0; j < rGeom.size(); j ++)
                        {

                            xg[dim] = xg[dim] + shape_functions_values(PointNumber, j) * rGeom[j].Coordinates()[dim];
                        }
                    }
                    //xg[0] = 0.5 * 1.923076923 + 1.923076923 * k;
                    //xg[1] = 1.923076923;

                    //xg[0] = -1.4626;
                    //xg[1] = -0.39545;
                    //std::cout<<"xg "<< xg<<std::endl;
                    p_element -> SetValue(GAUSS_COORD, xg);
                    p_element -> SetValue(MP_NUMBER, integration_point_per_elements);
                    //p_element -> SetValue(MP_CAUCHY_STRESS_VECTOR, MP_CauchyVector);
                    //p_element -> SetValue(MP_ALMANSI_STRAIN_VECTOR, MP_AlmansiVector);
                    //p_element -> SetValue(MP_CONSTITUTIVE_MATRIX, MP_ConstitutiveMatrix);
                    p_element -> SetValue(MP_DENSITY, MP_Density);

                    p_element -> SetValue(MP_MASS, MP_Mass);
                    p_element -> SetValue(MP_VOLUME, MP_Volume);

                    p_element -> SetValue(MP_DISPLACEMENT, MP_Displacement);
                    p_element -> SetValue(MP_VELOCITY, MP_Velocity);
                    //p_element -> SetValue(MP_KINETIC_ENERGY, MP_KineticEnergy);
                    //p_element -> SetValue(MP_STRAIN_ENERGY, MP_StrainEnergy);
                    //push back to model_part
                    mpm_model_part.Elements().push_back(p_element);


                }

                k +=1;

            }

        }
        //define a standard static strategy to be used in the calculation
        if(SolutionType == "StaticSolver")
        {


            typename TSchemeType::Pointer pscheme = typename TSchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace >() );
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );

            double ratio_tolerance = 1e-04;
            double always_converged_norm = 1e-09;
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));

            int MaxIterations = 20;
            bool CalculateReactions = false;
            bool ReformDofAtEachIteration = false;
            bool MoveMeshFlags = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new MPMResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofAtEachIteration,MoveMeshFlags) );
        }

        //define a dynamic strategy to be used in the calculation
        else if(SolutionType == "DynamicSolver")
        {
            double Alpham;
            double Dynamic;
            typename TSchemeType::Pointer pscheme = typename TSchemeType::Pointer( new MPMResidualBasedBossakScheme< TSparseSpace,TDenseSpace >(mr_grid_model_part, Alpham = 0.0, Dynamic=1) );

            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );

            double ratio_tolerance = 0.00005;
            double always_converged_norm = 1e-09;

            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));
            int MaxIterations = 20;
            bool CalculateReactions = false;
            bool ReformDofAtEachIteration = false;
            bool MoveMeshFlags = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new MPMResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofAtEachIteration,MoveMeshFlags) );
        }

        //define a quasi-static strategy to be used in the calculation
        else if(SolutionType == "QuasiStaticSolver")
        {
            double Alpham;
            double Dynamic;
            typename TSchemeType::Pointer pscheme = typename TSchemeType::Pointer( new MPMResidualBasedBossakScheme< TSparseSpace,TDenseSpace >(mr_grid_model_part, Alpham = 0.00, Dynamic=0) );

            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );

            double ratio_tolerance = 0.0001;
            double always_converged_norm = 1e-09;
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));

            int MaxIterations = 100;
            bool CalculateReactions = false;
            bool ReformDofAtEachIteration = false;
            bool MoveMeshFlags = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new MPMResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofAtEachIteration,MoveMeshFlags) );
        }




        initial_model_part.Nodes().clear();
        initial_model_part.Elements().clear();





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
    /*@{ */

    /**
    operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
    virtual void Predict()
    {
    }

    /**
    Initialization of member variables and prior operations
     */
    virtual void Initialize()
    {
    }

    /**
    the problem of interest is solved
     */
    virtual double Solve()
    {
        //check which nodes and elements are ACTIVE and populate the MPM model part
		//std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        
        
        
        this->SearchElement(mr_grid_model_part, mr_mpm_model_part);

        mp_solving_strategy->Initialize();

        //the nodal initial conditions are computed
        mp_solving_strategy->InitializeSolutionStep();

        mp_solving_strategy->Predict();
        //do solution iterations
        mp_solving_strategy->SolveSolutionStep();

        //the nodal solution are mapped on MP
        mp_solving_strategy->FinalizeSolutionStep();
        mp_solving_strategy->Clear();

        //std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();

		//std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() <<std::endl;
		//std::cout << "Time difference to solve (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0 <<std::endl;


        return 0.00;
    }

    virtual Matrix MP16ShapeFunctions()
    {
        double Na1 = 0.33333333333333;
        double Nb1 = 0.45929258829272;
        double Nb2 = 0.08141482341455;
        double Nc1 = 0.17056930775176;
        double Nc2 = 0.65886138449648;

        double Nd1 = 0.05054722831703;
        double Nd2 = 0.89890554336594;

        double Ne1 = 0.26311282963464;
        double Ne2 = 0.72849239295540;
        double Ne3 = 0.00839477740996;

        boost::numeric::ublas::bounded_matrix<double,16,3> MP_ShapeFunctions;// = ZeroMatrix(16,3);
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
        double Na2 = 0.02356522045239;
        double Na1 = 0.488217389773805;

        double Nb2 = 0.120551215411079;
        double Nb1 = 0.43972439229446;

        double Nc2 = 0.457579229975768;
        double Nc1 = 0.271210385012116;

        double Nd2 = 0.744847708916828;
        double Nd1 = 0.127576145541586;

        double Ne2 = 0.957365299093579;
        double Ne1 = 0.021317350453210;

        double Nf1 = 0.115343494534698;
        double Nf2 = 0.275713269685514;
        double Nf3 = 0.608943235779788;

        double Ng1 = 0.022838332222257;
        double Ng2 = 0.281325580989940;
        double Ng3 = 0.695836086787803;

        double Nh1 = 0.025734050548330;
        double Nh2 = 0.116251915907597;
        double Nh3 = 0.858014033544073;
        boost::numeric::ublas::bounded_matrix<double,33,3> MP_ShapeFunctions;// = ZeroMatrix(16,3);

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
     * 1) all the elements are set to be INACTIVE
     * 2) a serching is performed and the grid elements which contain at least a MP are set to be ACTIVE
     *
    */
    virtual void SearchElement(
        ModelPart& grid_model_part,
        ModelPart& mpm_model_part)
    {
        
        
        

        
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(grid_model_part.Elements().size()); ++i){
			
			auto elemItr = grid_model_part.Elements().begin() + i;
			elemItr -> Reset(ACTIVE);
            elemItr ->GetGeometry()[0].Reset(ACTIVE);
            elemItr ->GetGeometry()[1].Reset(ACTIVE);
            elemItr ->GetGeometry()[2].Reset(ACTIVE);
            if (TDim ==3)
            {

                elemItr ->GetGeometry()[3].Reset(ACTIVE);
            }
			
			
		}
        
        

        //******************SEARCH FOR TRIANGLES************************
        

        if (m_GeometryElement == "Triangle")
        {
            const int max_results = 1000;
            array_1d<double, TDim + 1 > N;

            BinBasedFastPointLocator<TDim> SearchStructure(grid_model_part);
            SearchStructure.UpdateSearchDatabase();

			#pragma omp parallel
			{
            typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);



			#pragma omp for
            for(int i = 0; i < static_cast<int>(mpm_model_part.Elements().size()); ++i){

				auto elemItr = mpm_model_part.Elements().begin() + i;
				
				array_1d<double,3> xg = elemItr -> GetValue(GAUSS_COORD);
				//KRATOS_WATCH(xg);
                typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                Element::Pointer pelem;

                //FindPointOnMesh find the element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg, N, pelem, result_begin, max_results);

                if (is_found == true)
                {
                    pelem->Set(ACTIVE);
                    

                    elemItr->GetGeometry()(0) = pelem->GetGeometry()(0);
                    elemItr->GetGeometry()(1) = pelem->GetGeometry()(1);
                    elemItr->GetGeometry()(2) = pelem->GetGeometry()(2);

                    pelem->GetGeometry()[0].Set(ACTIVE);
                    pelem->GetGeometry()[1].Set(ACTIVE);
                    pelem->GetGeometry()[2].Set(ACTIVE);
					
                    if (TDim ==3)
                    {

                        elemItr->GetGeometry()(3) = pelem->GetGeometry()(3);
                        pelem->GetGeometry()[3].Set(ACTIVE);
                    }
                    


                }


			}
		}

            

        }

      


        //******************SEARCH FOR QUADRILATERALS************************
        else if(m_GeometryElement == "Quadrilateral")
        {
            const int max_results = 1000;
            

            QuadBinBasedFastPointLocator<TDim> SearchStructure(grid_model_part);
            SearchStructure.UpdateSearchDatabase();

	    #pragma omp parallel
			{
            typename QuadBinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);



            //loop over the material points
            #pragma omp for
            for(int i = 0; i < static_cast<int>(mpm_model_part.Elements().size()); ++i){

		auto elemItr = mpm_model_part.Elements().begin() + i;
                
                array_1d<double,3> xg = elemItr -> GetValue(GAUSS_COORD);
                typename QuadBinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                Element::Pointer pelem;

                //FindPointOnMesh find the element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg, pelem, result_begin, max_results);
                
                if (is_found == true)
                {
                    pelem->Set(ACTIVE);
                    

                    elemItr->GetGeometry()(0) = pelem->GetGeometry()(0);
                    elemItr->GetGeometry()(1) = pelem->GetGeometry()(1);
                    elemItr->GetGeometry()(2) = pelem->GetGeometry()(2);
                    elemItr->GetGeometry()(3) = pelem->GetGeometry()(3);
                    
                    

                }


            }

        }
		
       }
        
    }




    /**
     * function to perform expensive checks.
     * It is designed to be called ONCE to verify that the input is correct.
     */
    virtual int Check()
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

