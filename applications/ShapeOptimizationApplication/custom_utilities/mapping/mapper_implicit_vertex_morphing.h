// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//
// ==============================================================================

#ifndef MAPPER_IMPLICIT_VERTEX_MORPHING_H
#define MAPPER_IMPLICIT_VERTEX_MORPHING_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/builtin_timer.h"
#include "spaces/ublas_space.h"
#include "mapper_base.h"
#include "custom_elements/helmholtz_element.h"
#include "custom_elements/helmholtz_vec_element.h"
#include "custom_strategies/strategies/helmholtz_strategy.h"
#include "custom_strategies/strategies/helmholtz_vec_strategy.h"
#include "containers/model.h"
#include "linear_solvers/linear_solver.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/

class MapperImplicitVertexMorphing : public Mapper
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef array_1d<double,3> array_3d;
    typedef Element BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::NodesArrayType NodesArrayType;
    typedef BaseType::PropertiesType PropertiesType;
    typedef BaseType::IndexType IndexType;
    typedef BaseType::SizeType SizeType;    
    typedef BaseType::MatrixType MatrixType;
    typedef BaseType::VectorType VectorType;    
    typedef GeometryData::IntegrationMethod IntegrationMethod;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;    

    /// Pointer definition of MapperImplicitVertexMorphing
    KRATOS_CLASS_POINTER_DEFINITION(MapperImplicitVertexMorphing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperImplicitVertexMorphing( ModelPart& rModelPart, LinearSolverType::Pointer pLinearSolver, Parameters MapperSettings )
        : mrModelPart(rModelPart), mpLinearSystemSolver(pLinearSolver),
          mMapperSettings(MapperSettings)
    {
    }

    /// Destructor.
    virtual ~MapperImplicitVertexMorphing()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Initialize() override
    {
        BuiltinTimer timer;
        KRATOS_INFO("ShapeOpt") << "Starting initialization of mapper..." << std::endl;

        std::string element_type = mMapperSettings["element_type"].GetString();

        // here we create a model part for implicit VM
        mpVMModePart = &(mrModelPart.GetModel().CreateModelPart(mrModelPart.Name()+"_Implicit_VM_Part", 1));

        // initializing vm model nodes and variables
        mpVMModePart->Nodes() = mrModelPart.Nodes();

        // creating vm elements
        ModelPart::ElementsContainerType &rmesh_elements =
            mpVMModePart->Elements();

        // create a new property for the vm
        Properties::Pointer p_vm_property = mpVMModePart->CreateNewProperties(0);        
        p_vm_property->SetValue(HELMHOLTZ_RADIUS,mMapperSettings["filter_radius"].GetDouble());

        for (int i = 0; i < (int)mrModelPart.Elements().size(); i++) {
            ModelPart::ElementsContainerType::iterator it = mrModelPart.ElementsBegin() + i;
            Element::Pointer p_element;
            if (element_type.compare("helmholtz_element") == 0)
                p_element = new HelmholtzElement(it->Id(), it->pGetGeometry(), p_vm_property);
            else
                p_element = new HelmholtzVecElement(it->Id(), it->pGetGeometry(), p_vm_property);                
            rmesh_elements.push_back(p_element);
        }

        // calculate number of neighbour elements for each node.
        CalculateNodeNeighbourCount();


        if (element_type.compare("helmholtz_element") == 0)
            mpHelmholtzStrategy = new HelmholtzStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType> (*mpVMModePart,mpLinearSystemSolver);
        else 
            mpHelmholtzStrategy = new HelmholtzVecStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType> (*mpVMModePart,mpLinearSystemSolver);            
        
        mpHelmholtzStrategy->Initialize();



        ProcessInfo &rCurrentProcessInfo = mpVMModePart->GetProcessInfo();
        rCurrentProcessInfo[COMPUTE_CONTROL_POINTS] = true;

        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_coords = node_i.Coordinates();
            array_3d& r_nodal_variable_hl_source = node_i.FastGetSolutionStepValue(HELMHOLTZ_SOURCE);
            r_nodal_variable_hl_source(0) = r_nodal_coords(0);
            r_nodal_variable_hl_source(1) = r_nodal_coords(1);
            r_nodal_variable_hl_source(2) = r_nodal_coords(2);
        }

        mpHelmholtzStrategy->Solve();

        //set the values of control points
        SetVariable1ToVarible2(HELMHOLTZ_VARS,CONTROL_POINT);

        //set back the COMPUTE_CONTROL_POINTS flag to false
        rCurrentProcessInfo[COMPUTE_CONTROL_POINTS] = false;
        mIsMappingInitialized = true;

        Update();

        KRATOS_INFO("ShapeOpt") << "Finished initialization of mapper in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

        //first we need to multiply with mass matrix 
        SetVariableZero(HELMHOLTZ_VARS);
        SetVariableZero(HELMHOLTZ_SOURCE);
        //now we need to multiply with the mass matrix 
        for(auto& elem_i : mpVMModePart->Elements())
        {
            VectorType origin_values;
            MatrixType mass_matrix;
            CalculateMassMatrix(elem_i,mass_matrix);
            GetElementVariableValuesVector(elem_i,rOriginVariable,origin_values);
            VectorType int_vals = prod(mass_matrix,origin_values);
            AddElementVariableValuesVector(elem_i,HELMHOLTZ_SOURCE,int_vals);
        }

        mpHelmholtzStrategy->Solve();

        //filling the solution
        SetVariable1ToVarible2(HELMHOLTZ_VARS,rDestinationVariable); 

        KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override
    {
        KRATOS_ERROR << "Scalar mapping not possible." << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "Starting inverse mapping of " << rDestinationVariable.Name() << "..." << std::endl;
        
        //filling the source
        SetVariableZero(HELMHOLTZ_VARS);
        SetVariable1ToVarible2(rDestinationVariable,HELMHOLTZ_SOURCE);

        //now solve 
        mpHelmholtzStrategy->Solve();

        //first we set origin values to zero
        SetVariableZero(rOriginVariable);
        //now we need to multiply with the mass matrix with results to be consistence.
        for(auto& elem_i : mpVMModePart->Elements())
        {
            VectorType helmholtz_values;
            MatrixType mass_matrix;
            CalculateMassMatrix(elem_i,mass_matrix);
            GetElementVariableValuesVector(elem_i,HELMHOLTZ_VARS,helmholtz_values);
            VectorType int_vals = prod(mass_matrix,helmholtz_values);
            AddElementVariableValuesVector(elem_i,rOriginVariable,int_vals);
        }

    
        KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) override
    {
        KRATOS_ERROR << "Scalar mapping not possible." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Update() override
    {
        if (mIsMappingInitialized == false)
            KRATOS_ERROR << "Mapping has to be initialized before calling the Update-function!";
    }

    // --------------------------------------------------------------------------

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "MapperImplicitVertexMorphing";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperImplicitVertexMorphing";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    // Initialized by class constructor
    ModelPart& mrModelPart;
    LinearSolverType::Pointer mpLinearSystemSolver = nullptr;
    Parameters mMapperSettings;
    bool mIsMappingInitialized = false;
    ModelPart* mpVMModePart;
    ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType>* mpHelmholtzStrategy;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CalculateMassMatrix(
        const Element& rElement,
        MatrixType& rMassMatrix
        )
    {
        KRATOS_TRY;

        const auto& r_geom = rElement.GetGeometry();
        SizeType dimension = r_geom.WorkingSpaceDimension();
        SizeType number_of_nodes = r_geom.size();
        SizeType mat_size = dimension * number_of_nodes;

        // Clear matrix
        if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
            rMassMatrix.resize( mat_size, mat_size, false );
        rMassMatrix = ZeroMatrix( mat_size, mat_size );

        // CONSISTENT MASS

        Element::GeometryType::JacobiansType J0;
        r_geom.Jacobian(J0,r_geom.GetDefaultIntegrationMethod());
        Matrix InvJ0(dimension,dimension);
        double detJ0;


        const IntegrationMethod integration_method = r_geom.GetDefaultIntegrationMethod();
        const auto& integration_points = r_geom.IntegrationPoints(integration_method);
        const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            //calculating inverse jacobian and jacobian determinant
            MathUtils<double>::InvertMatrix(J0[point_number],InvJ0,detJ0);;
            const double integration_weight = integration_points[point_number].Weight() * detJ0;
            
            const Vector& rN = row(Ncontainer,point_number);

            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                const SizeType index_i = i * dimension;

                for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                    const SizeType index_j = j * dimension;
                    const double NiNj_weight = rN[i] * rN[j] * integration_weight;

                    for ( IndexType k = 0; k < dimension; ++k )
                        rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void GetElementVariableValuesVector(const Element& rElement,
                                        const Variable<array_3d> &rVariable,
                                        VectorType &rValues) const
    {
        const GeometryType &rgeom = rElement.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const unsigned int dimension = rElement.GetGeometry().WorkingSpaceDimension();
        const unsigned int local_size = num_nodes * dimension;

        if (rValues.size() != local_size)
            rValues.resize(local_size, false);

        if (dimension == 2) {
            SizeType index = 0;
            for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
                const array_3d& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable);    
                rValues[index++] = r_nodal_variable[0];
                rValues[index++] = r_nodal_variable[1];
            }
        } else if (dimension == 3) {
            SizeType index = 0;
            for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
                const array_3d& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable);    
                rValues[index++] = r_nodal_variable[0];
                rValues[index++] = r_nodal_variable[1];
                rValues[index++] = r_nodal_variable[2];
            }
        }
    }
    void AddElementVariableValuesVector(Element& rElement,
                                        const Variable<array_3d> &rVariable,
                                        const VectorType &rValues
                                        ) 
    {
        GeometryType &rgeom = rElement.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const unsigned int dimension = rElement.GetGeometry().WorkingSpaceDimension();

        if (dimension == 2) {
            SizeType index = 0;
            for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
                array_3d& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable); 
                r_nodal_variable[0] += rValues[index++];
                r_nodal_variable[1] += rValues[index++];
            }
        } else if (dimension == 3) {
            SizeType index = 0;
            for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
                array_3d& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable);
                r_nodal_variable[0] += rValues[index++];
                r_nodal_variable[1] += rValues[index++];
                r_nodal_variable[2] += rValues[index++];
            }
        }
    }
    
    void SetVariableZero(const Variable<array_3d> &rVariable) 
    {
        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_variable = node_i.FastGetSolutionStepValue(rVariable);
            r_nodal_variable[0] = 0.0;
            r_nodal_variable[1] = 0.0;
            r_nodal_variable[2] = 0.0;                    
        }
    }

    void SetVariable1ToVarible2(const Variable<array_3d> &rVariable1,const Variable<array_3d> &rVariable2) 
    {
        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_variable1 = node_i.FastGetSolutionStepValue(rVariable1);
            array_3d& r_nodal_variable2 = node_i.FastGetSolutionStepValue(rVariable2);
            r_nodal_variable2[0] = r_nodal_variable1[0];
            r_nodal_variable2[1] = r_nodal_variable1[1];
            r_nodal_variable2[2] = r_nodal_variable1[2];                    
        }
    }

    void CalculateNodeNeighbourCount()
    {

        auto& r_nodes = mpVMModePart->Nodes();
        int mNumNodes = r_nodes.size();

        VariableUtils variable_utils;
        variable_utils.SetFlag(STRUCTURE,true,r_nodes);

        // Note: this should not be parallel, the operation is not threadsafe if the variable is uninitialized
        for (auto& r_node : r_nodes)
        {
            r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS,0);
        }

        mNumNodes = mpVMModePart->GetCommunicator().GetDataCommunicator().SumAll(mNumNodes);

        auto& r_elements = mpVMModePart->Elements();
        const int num_elements = r_elements.size();

        #pragma omp parallel for
        for (int i = 0; i < num_elements; i++)
        {
            auto i_elem = r_elements.begin() + i;
            auto& r_geom = i_elem->GetGeometry();
            for (unsigned int i = 0; i < r_geom.PointsNumber(); i++)
            {
                auto& r_node = r_geom[i];
                if (r_node.Is(STRUCTURE))
                {
                    r_node.SetLock();
                    r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS) += 1;
                    r_node.UnSetLock();
                }
            }
        }

        mpVMModePart->GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);

    }            
            

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//      MapperImplicitVertexMorphing& operator=(MapperImplicitVertexMorphing const& rOther);

    /// Copy constructor.
//      MapperImplicitVertexMorphing(MapperImplicitVertexMorphing const& rOther);


    ///@}

}; // Class MapperImplicitVertexMorphing

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_IMPLICIT_VERTEX_MORPHING_H
