// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef IMPLICIT_VERTEX_MORPHING_H
#define IMPLICIT_VERTEX_MORPHING_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------

#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/tetrahedral_mesh_orientation_check.h"
#include "utilities/builtin_timer.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "input_output/vtk_output.h"
#include "containers/model.h"
#include "custom_controls/shape_controls/shape_control.h"
#include "custom_elements/helmholtz_surf_shape_element.h"

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

class KRATOS_API(OPTIMIZATION_APPLICATION) ImplicitVertexMorphing : public ShapeControl
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

    /// Pointer definition of ImplicitVertexMorphing
    KRATOS_CLASS_POINTER_DEFINITION(ImplicitVertexMorphing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ImplicitVertexMorphing( std::string ControlName, Model& rModel, LinearSolverType::Pointer pLinearSolver, Parameters ControlSettings )
        : mpLinearSystemSolver(pLinearSolver), ShapeControl(ControlName,rModel,ControlSettings){}

    /// Destructor.
    virtual ~ImplicitVertexMorphing()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Initialize() override {

        BuiltinTimer timer;
        KRATOS_INFO("ImplicitVertexMorphing:Initialize ") << "Starting initialization of shape control "<<mControlName<<" ..." << std::endl;

        bool only_suf_param =  mControlSettings["technique_settings"]["only_design_surface_parameterization"].GetBool();


        //  std::cout<<mrModel.Info();
        // initializing vm model nodes and variables
        for(auto& control_obj : mControlSettings["controlling_objects"]){
            ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
            ModelPart& root_model_part = r_controlling_object.GetRootModelPart();
            std::string vm_model_part_name =  root_model_part.Name()+"_Implicit_VM_Part";
            ModelPart* p_vm_model_part;

            if (root_model_part.HasSubModelPart(vm_model_part_name))
                p_vm_model_part = &(root_model_part.GetSubModelPart(vm_model_part_name));
            else{
                p_vm_model_part = &(root_model_part.CreateSubModelPart(vm_model_part_name));
                mpVMModelParts.push_back(p_vm_model_part);
            }

            for(auto& node : r_controlling_object.Nodes())
                p_vm_model_part->AddNode(&node);

        }

        for(int i=0;i<mpVMModelParts.size();i++)
        {
            std::cout<<" i = "<<i<<", name : "<<mpVMModelParts[i]->Name()<<std::endl;
        }

        // std::cout<<mrModel.Info();



        // // creating vm elements
        // ModelPart::ElementsContainerType &rmesh_elements = mpVMModePart->Elements();  

        // // creating vm conditions
        // ModelPart::ConditionsContainerType &rmesh_conditions = mpVMModePart->Conditions();         

        // bool only_suf_param =  mMapperSettings["only_design_surface_parameterization"].GetBool();
        // ModelPart& design_surface_sub_model_part = mrModelPart.GetSubModelPart(mMapperSettings["design_surface_sub_model_part_name"].GetString());



        std::exit(0);

        const std::vector<std::string>& controlling_objects = mControlSettings["controlling_objects"].GetStringArray();
    };
    // --------------------------------------------------------------------------
    void Update() override {};  
    // --------------------------------------------------------------------------
    void MapControlUpdate(const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override{};
    // --------------------------------------------------------------------------
    void MapFirstDerivative(const Variable<array_3d> &rDerivativeVariable, const Variable<array_3d> &rMappedDerivativeVariable) override{};  

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
        return "ImplicitVertexMorphing";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ImplicitVertexMorphing";
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
    LinearSolverType::Pointer mpLinearSystemSolver = nullptr;
    std::vector<ModelPart*> mpVMModelParts;
    
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
    // std::vector<Kratos::unique_ptr<MapperVertexMorphing>> mMappers;


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
//      ImplicitVertexMorphing& operator=(ImplicitVertexMorphing const& rOther);

    /// Copy constructor.
//      ImplicitVertexMorphing(ImplicitVertexMorphing const& rOther);


    ///@}

}; // Class ImplicitVertexMorphing

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // IMPLICIT_VERTEX_MORPHING_H
