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
//#include "strucural_application.h"
#include "particle_mechanics_application.h"

#include "custom_elements/updated_lagrangian.hpp"
#include "custom_strategies/custom_schemes/residual_based_static_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_bossak_scheme.hpp"
#include "custom_strategies/schemes/MPM_residual_based_bossak_scheme.hpp"
#include "custom_strategies/custom_schemes/residual_based_static_scheme.hpp"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
//#include "custom_strategies/schemes/residualbased_predictorcorrector_bossak_scheme.h"
#include "custom_strategies/custom_builders_and_solvers/residual_based_builder_and_solver.hpp"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergence_criteria/displacement_convergence_criterion.hpp"

#include "custom_strategies/strategies/MPM_residual_based_newton_raphson_strategy.hpp"
#include "custom_strategies/strategies/MPM_strategy.h"

#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "linear_solvers/linear_solver.h"

#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
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
    
    MPMStrategy(ModelPart& grid_model_part, ModelPart& mpm_model_part, typename TLinearSolver::Pointer plinear_solver,Element const& NewElement, bool MoveMeshFlag = false, std::string SolutionType = "StaticType")
    : SolvingStrategyType(grid_model_part, MoveMeshFlag), mr_grid_model_part(grid_model_part), mr_mpm_model_part(mpm_model_part)
    {   
        
         //populate for the first time the mpm_model_part
        
        //assigning the nodes to the new model part
        mpm_model_part.Nodes() = grid_model_part.Nodes();
        
        mpm_model_part.SetProcessInfo(grid_model_part.pGetProcessInfo());
        mpm_model_part.SetBufferSize(grid_model_part.GetBufferSize());
        mpm_model_part.SetProperties(grid_model_part.pProperties());
        mpm_model_part.SetConditions(grid_model_part.pConditions());
        
        
        
        
        array_1d<double,3> xg ;
        array_1d<double,3> MP_Displacement;
        array_1d<double,3> MP_Velocity;
        double MP_KineticEnergy = 0.0;
        double MP_StrainEnergy = 0.0;
        Vector MP_CauchyVector = ZeroVector(3);
        Vector MP_AlmansiVector = ZeroVector(3);
        Matrix MP_ConstitutiveMatrix = ZeroMatrix(6,6);
        double MP_Mass;
        double MP_Volume;
        
        unsigned int k = 0;
        const unsigned int number_elements = grid_model_part.NumberOfElements();
        
        
        for (ModelPart::ElementIterator i = grid_model_part.ElementsBegin();
                    i != grid_model_part.ElementsEnd(); i++)
            {   
                if(i->IsDefined(ACTIVE))
                {
                    
                    
                    Properties::Pointer properties = i->pGetProperties();
                    double Density = i->GetProperties()[DENSITY];  
                    Geometry< Node < 3 > >& rGeom = i->GetGeometry(); // current element's connectivity
                                   
                    Matrix shape_functions_values = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_1);
                    
                    //std::cout<<"shape_functions_values "<< shape_functions_values<<std::endl;
                    //const GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints( GeometryData::GI_GAUSS_4); 
                    
                    //INITIAL NUMBER OF MATERIAL POINTS PER ELEMENT
                    unsigned int integration_point_per_elements = shape_functions_values.size1();

                    
                    //evaluation of element area/volume
                    double area = GeometryUtils::CalculateVolume2D(rGeom);
                    
                    //int integration_point_per_element = 3;
                    
                    MP_Mass = area * Density / integration_point_per_elements;
                    
                    MP_Volume = area / integration_point_per_elements;
                    
                    
                    //loop over the material points that fall in each grid element
                    for ( unsigned int PointNumber = 0; PointNumber < integration_point_per_elements; PointNumber++ ) 
                        {   
                                                    
                            Element::Pointer p_element = NewElement.Create((1+PointNumber+number_elements)+(integration_point_per_elements*k), rGeom, properties);
                            //Element::Pointer p_element = NewElement.Create((1+PointNumber+number_nodes)+(integration_point_per_element*k), rGeom, properties);
                            
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
                            p_element -> SetValue(MP_CAUCHY_STRESS_VECTOR, MP_CauchyVector);
                            p_element -> SetValue(MP_ALMANSI_STRAIN_VECTOR, MP_AlmansiVector);
                            p_element -> SetValue(MP_CONSTITUTIVE_MATRIX, MP_ConstitutiveMatrix);
                            
                            
                            p_element -> SetValue(MP_MASS, MP_Mass);
                            p_element -> SetValue(MP_VOLUME, MP_Volume);
                            
                            p_element -> SetValue(MP_DISPLACEMENT, MP_Displacement);
                            p_element -> SetValue(MP_VELOCITY, MP_Velocity);
                            p_element -> SetValue(MP_KINETIC_ENERGY, MP_KineticEnergy);
                            p_element -> SetValue(MP_STRAIN_ENERGY, MP_StrainEnergy);
                            //push back to model_part
                            mpm_model_part.Elements().push_back(p_element);
                        
                        
                        }
                        
                    k +=1;
                    
                }
            
            }       
        //define a standard static strategy to be used in the calculation
        if(SolutionType == "StaticSolver")
        {
            
            
            typename TSchemeType::Pointer pscheme = typename TSchemeType::Pointer( new ResidualBasedStaticScheme< TSparseSpace,TDenseSpace >() );
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            
            double ratio_tolerance = 0.0001;
            double always_converged_norm = 1e-09;
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));
            
            int MaxIterations = 2;
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
           
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            
            double ratio_tolerance = 0.0001;
            double always_converged_norm = 1e-09;
            
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));
            int MaxIterations = 10;
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
            
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            
            double ratio_tolerance = 0.0001;
            double always_converged_norm = 1e-09;
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));
            
            int MaxIterations = 10;
            bool CalculateReactions = false;
            bool ReformDofAtEachIteration = false;
            bool MoveMeshFlags = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new MPMResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofAtEachIteration,MoveMeshFlags) );
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
        
        //std::cout<<" clear the solving strategy"<<std::endl;
        
        
        return 0.00;
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
    
        //Set all the grid elements to be inactive
        for (ModelPart::ElementIterator i = grid_model_part.ElementsBegin();
                i != grid_model_part.ElementsEnd(); ++i)
            {   
                
                i -> Reset(ACTIVE);
                i -> SetValue(COUNTER, 0);     
            }
                   
            
        const int max_results = 1000; 
        array_1d<double, TDim + 1 > N;
        
        BinBasedFastPointLocator<TDim> SearchStructure(grid_model_part);
        SearchStructure.UpdateSearchDatabase();
        
        
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);     
        
        
        
        //loop over the material points
        for (ModelPart::ElementIterator k = mpm_model_part.ElementsBegin(); 
                k != mpm_model_part.ElementsEnd(); ++k)
            {   
                //const unsigned int number_of_nodes = k -> GetGeometry().PointsNumber();
                array_1d<double,3> xg = k -> GetValue(GAUSS_COORD);
                typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
                                    
                Element::Pointer pelem;
                
                //FindPointOnMesh find the element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg, N, pelem, result_begin, max_results);                        
                
                if (is_found == true)
                {
                pelem->Set(ACTIVE);
                //std::cout<<"pelem->Id() "<<pelem->Id()<<std::endl;
                int counter = pelem -> GetValue(COUNTER);            
                counter += 1;
                pelem -> SetValue(COUNTER, counter);
                
                k->GetGeometry()(0) = pelem->GetGeometry()(0);
                k->GetGeometry()(1) = pelem->GetGeometry()(1);
                k->GetGeometry()(2) = pelem->GetGeometry()(2);
                
                //if use quadrilateral add this line
                //k->GetGeometry()(3) = pelem->GetGeometry()(3);
                
                
                
                
                }
                
                
            } 
            
        //loop over grid element to know how many MP fall in each element
        for (ModelPart::ElementIterator i = grid_model_part.ElementsBegin();
                i != grid_model_part.ElementsEnd(); ++i)
            { 
                if(i->IsDefined(ACTIVE))
                {
                    int MP_number_per_element = i->GetValue(COUNTER);
                    for (ModelPart::ElementIterator k = mpm_model_part.ElementsBegin(); 
                        k != mpm_model_part.ElementsEnd(); ++k)
                        {
                            if (k->GetGeometry()(0) == i->GetGeometry()(0) && k->GetGeometry()(1) == i->GetGeometry()(1) && k->GetGeometry()(2) == i->GetGeometry()(2))
                            {
                                k->SetValue(MP_NUMBER, MP_number_per_element);
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

#endif /* KRATOS_NEW_SOLVING_STRATEGY  defined */

