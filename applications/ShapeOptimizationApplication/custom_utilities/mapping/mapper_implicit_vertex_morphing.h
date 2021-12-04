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
#include "custom_strategies/strategies/helmholtz_strategy.h"
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
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef array_1d<double,3> array_3d;
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


        //here we create a model part for implicit VM
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
            ModelPart::ElementsContainerType::iterator it =
                mrModelPart.ElementsBegin() + i;
            Element::Pointer p_element = new HelmholtzElement(it->Id(), it->pGetGeometry(), p_vm_property);
            rmesh_elements.push_back(p_element);
        }

        mpHelmholtzStrategy = new HelmholtzStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType> (*mpVMModePart,mpLinearSystemSolver);
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

        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_variable_control = node_i.FastGetSolutionStepValue(CONTROL_POINT);
            array_3d& r_nodal_variable_hl_vars = node_i.FastGetSolutionStepValue(HELMHOLTZ_VARS);
            r_nodal_variable_control(0) = r_nodal_variable_hl_vars(0);
            r_nodal_variable_control(1) = r_nodal_variable_hl_vars(1);
            r_nodal_variable_control(2) = r_nodal_variable_hl_vars(2);
        } 

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

        //filling the source
        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_variable_org = node_i.FastGetSolutionStepValue(rOriginVariable);
            array_3d& r_nodal_variable_hl_source = node_i.FastGetSolutionStepValue(HELMHOLTZ_SOURCE);
            array_3d& r_nodal_variable_hl_vars = node_i.FastGetSolutionStepValue(HELMHOLTZ_VARS);
            r_nodal_variable_hl_source(0) = r_nodal_variable_org(0);
            r_nodal_variable_hl_source(1) = r_nodal_variable_org(1);
            r_nodal_variable_hl_source(2) = r_nodal_variable_org(2);
            r_nodal_variable_hl_vars(0) = 0.0;
            r_nodal_variable_hl_vars(1) = 0.0;
            r_nodal_variable_hl_vars(2) = 0.0;
        }

        mpHelmholtzStrategy->Solve();

        //filling the solution
        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_variable_des = node_i.FastGetSolutionStepValue(rDestinationVariable);
            array_3d& r_nodal_variable_hl_vars = node_i.FastGetSolutionStepValue(HELMHOLTZ_VARS);
            r_nodal_variable_des(0) = r_nodal_variable_hl_vars(0);
            r_nodal_variable_des(1) = r_nodal_variable_hl_vars(1);
            r_nodal_variable_des(2) = r_nodal_variable_hl_vars(2);
        }  

        KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

        KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
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
        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_variable_des = node_i.FastGetSolutionStepValue(rDestinationVariable);
            array_3d& r_nodal_variable_hl_source = node_i.FastGetSolutionStepValue(HELMHOLTZ_SOURCE);
            array_3d& r_nodal_variable_hl_vars = node_i.FastGetSolutionStepValue(HELMHOLTZ_VARS);
            r_nodal_variable_hl_source(0) = r_nodal_variable_des(0);
            r_nodal_variable_hl_source(1) = r_nodal_variable_des(1);
            r_nodal_variable_hl_source(2) = r_nodal_variable_des(2);
            r_nodal_variable_hl_vars(0) = 0.0;
            r_nodal_variable_hl_vars(1) = 0.0;
            r_nodal_variable_hl_vars(2) = 0.0;            
        }

        mpHelmholtzStrategy->Solve();

        //filling the solution
        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_variable_origin = node_i.FastGetSolutionStepValue(rOriginVariable);
            array_3d& r_nodal_variable_hl_vars = node_i.FastGetSolutionStepValue(HELMHOLTZ_VARS);
            r_nodal_variable_origin(0) = r_nodal_variable_hl_vars(0);
            r_nodal_variable_origin(1) = r_nodal_variable_hl_vars(1);
            r_nodal_variable_origin(2) = r_nodal_variable_hl_vars(2);
        }        

        KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "Starting inverse mapping of " << rDestinationVariable.Name() << "..." << std::endl;


        KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Update() override
    {
        if (mIsMappingInitialized == false)
            KRATOS_ERROR << "Mapping has to be initialized before calling the Update-function!";

        BuiltinTimer timer;
        KRATOS_INFO("ShapeOpt") << "Starting to update mapper..." << std::endl;


        KRATOS_INFO("ShapeOpt") << "Finished updating of mapper in " << timer.ElapsedSeconds() << " s." << std::endl;
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
    HelmholtzStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType>* mpHelmholtzStrategy;

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

    // --------------------------------------------------------------------------

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
