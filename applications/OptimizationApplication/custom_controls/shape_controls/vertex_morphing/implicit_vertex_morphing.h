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
#include "utilities/variable_utils.h"
#include "custom_controls/shape_controls/shape_control.h"
#include "custom_elements/helmholtz_surf_shape_element.h"
#include "custom_elements/helmholtz_bulk_shape_element.h"
#include "custom_conditions/helmholtz_surf_shape_condition.h"
#include "custom_strategies/strategies/helmholtz_strategy.h"


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
    typedef HelmholtzStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType> StrategyType;

    /// Pointer definition of ImplicitVertexMorphing
    KRATOS_CLASS_POINTER_DEFINITION(ImplicitVertexMorphing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ImplicitVertexMorphing( std::string ControlName, Model& rModel, std::vector<LinearSolverType::Pointer>& rLinearSolvers, Parameters ControlSettings )
        :  ShapeControl(ControlName,rModel,ControlSettings){
            for(int lin_i=0;lin_i<rLinearSolvers.size();lin_i++)
                rLinearSystemSolvers.push_back(rLinearSolvers[lin_i]);
            mTechniqueSettings = ControlSettings["technique_settings"];
        }

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

        CreateModelParts();

        CalculateNodeNeighbourCount();
        
        for(int model_i=0;model_i<mpVMModelParts.size();model_i++){
            StrategyType* mpStrategy = new StrategyType (*mpVMModelParts[model_i],rLinearSystemSolvers[model_i]);            
            mpStrategy->Initialize();
            mpStrategies.push_back(mpStrategy);
        }

        AdjustFilterSizes();

        ComputeInitialControlPoints();

        KRATOS_INFO("ImplicitVertexMorphing:Initialize") << "Finished initialization of shape control "<<mControlName<<" in " << timer.ElapsedSeconds() << " s." << std::endl;

    };
    // --------------------------------------------------------------------------
    void Update() override {
        AdjustFilterSizes();
    };  
    // --------------------------------------------------------------------------
    void MapControlUpdate(const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override{

        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ImplicitVertexMorphing:MapControlUpdate:") << " Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

        for(int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];

            //first we need to multiply with mass matrix 
            SetVariableZero(mpVMModePart,HELMHOLTZ_VARS_SHAPE);
            SetVariableZero(mpVMModePart,HELMHOLTZ_SOURCE_SHAPE);
            //now we need to multiply with the mass matrix 
            for(auto& elem_i : mpVMModePart->Elements())
            {
                VectorType origin_values;
                GetElementVariableValuesVector(elem_i,rOriginVariable,origin_values);
                MatrixType mass_matrix;
                elem_i.Calculate(HELMHOLTZ_MASS_MATRIX,mass_matrix,mpVMModePart->GetProcessInfo());            
                VectorType int_vals = prod(mass_matrix,origin_values);
                AddElementVariableValuesVector(elem_i,HELMHOLTZ_SOURCE_SHAPE,int_vals);
            }

            for(auto& cond_i : mpVMModePart->Conditions())
            {
                VectorType origin_values;
                GetConditionVariableValuesVector(cond_i,rOriginVariable,origin_values);
                MatrixType mass_matrix;
                cond_i.Calculate(HELMHOLTZ_MASS_MATRIX,mass_matrix,mpVMModePart->GetProcessInfo());            
                VectorType int_vals = prod(mass_matrix,origin_values);
                AddConditionVariableValuesVector(cond_i,HELMHOLTZ_SOURCE_SHAPE,int_vals);
            }

            mpStrategies[model_i]->Solve();

            //filling the solution
            SetVariable1ToVarible2(mpVMModePart,HELMHOLTZ_VARS_SHAPE,rDestinationVariable); 

        }

        KRATOS_INFO("ImplicitVertexMorphing:MapControlUpdate:") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;

    };
    // --------------------------------------------------------------------------
    void MapFirstDerivative(const Variable<array_3d> &rDerivativeVariable, const Variable<array_3d> &rMappedDerivativeVariable) override{


        AdjustFilterSizes();

        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ImplicitVertexMorphing:MapFirstDerivative") << "Starting mapping of " << rDerivativeVariable.Name() << "..." << std::endl;

        for(int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];    

            //filling the source
            SetVariableZero(mpVMModePart,HELMHOLTZ_VARS_SHAPE);
            SetVariable1ToVarible2(mpVMModePart,rDerivativeVariable,HELMHOLTZ_SOURCE_SHAPE);

            //now solve 
            mpStrategies[model_i]->Solve();

            SetVariable1ToVarible2(mpVMModePart,HELMHOLTZ_VARS_SHAPE,rMappedDerivativeVariable);
        }

    
        KRATOS_INFO("ImplicitVertexMorphing:MapFirstDerivative") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;

    };  

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
    std::vector<LinearSolverType::Pointer> rLinearSystemSolvers;
    std::vector<StrategyType*>mpStrategies;
    std::vector<ModelPart*> mpVMModelParts;
    std::vector<Properties::Pointer> mpVMModelPartsProperties;
    Parameters mTechniqueSettings;
    
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

    void CalculateNodeNeighbourCount()
    {
        for(int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];
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
    } 

    void CreateModelParts()
    {
        bool only_suf_param =  mTechniqueSettings["only_design_surface_parameterization"].GetBool();
        // creating vm model nodes and variables
        for(auto& control_obj : mControlSettings["controlling_objects"]){
            ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
            ModelPart& root_model_part = r_controlling_object.GetRootModelPart();
            const std::size_t domain_size = root_model_part.GetProcessInfo()[DOMAIN_SIZE];

            // check if control_obj has surface condition and root_model_part has elements
            KRATOS_ERROR_IF_NOT(r_controlling_object.Conditions().size()>0)
            <<"ImplicitVertexMorphing::CreateModelParts: controlling object "<<control_obj.GetString()<<" must have surface conditions !"<<std::endl;
            KRATOS_ERROR_IF_NOT(root_model_part.Elements().size()>0)
            <<"ImplicitVertexMorphing::CreateModelParts: root model of controlling object "<<control_obj.GetString()<<" must have 3D elements !"<<std::endl;
            KRATOS_ERROR_IF_NOT(domain_size == 3)
            << "ImplicitVertexMorphing::CreateModelParts: controlling_object should be a 3D model part " << std::endl;

            std::string vm_model_part_name =  root_model_part.Name()+"_Implicit_VM_Part";
            ModelPart* p_vm_model_part;
            Properties::Pointer p_vm_model_part_property;

            if (root_model_part.HasSubModelPart(vm_model_part_name)){
                p_vm_model_part = &(root_model_part.GetSubModelPart(vm_model_part_name));
                for(int i =0; i<mpVMModelParts.size(); i++)
                    if(mpVMModelParts[i]->Name()==p_vm_model_part->Name())
                        p_vm_model_part_property = mpVMModelPartsProperties[i];
            }
            else{
                p_vm_model_part = &(root_model_part.CreateSubModelPart(vm_model_part_name));
                p_vm_model_part_property = p_vm_model_part->CreateNewProperties(root_model_part.NumberOfProperties()+1);
                mpVMModelPartsProperties.push_back(p_vm_model_part_property);
                mpVMModelParts.push_back(p_vm_model_part);
            }


            p_vm_model_part_property->SetValue(HELMHOLTZ_SURF_RADIUS_SHAPE,mTechniqueSettings["surface_filter_radius"].GetDouble());
            p_vm_model_part_property->SetValue(HELMHOLTZ_BULK_RADIUS_SHAPE,mTechniqueSettings["bulk_filter_radius"].GetDouble()); 
            p_vm_model_part_property->SetValue(HELMHOLTZ_SURF_POISSON_RATIO_SHAPE,mTechniqueSettings["poisson_ratio"].GetDouble()); 
            p_vm_model_part_property->SetValue(HELMHOLTZ_BULK_POISSON_RATIO_SHAPE,mTechniqueSettings["poisson_ratio"].GetDouble());
            
            if(only_suf_param)
                for(auto& node : r_controlling_object.Nodes())
                    p_vm_model_part->AddNode(&node);
            else
                for(auto& node : root_model_part.Nodes())
                    p_vm_model_part->AddNode(&node);  


            // creating vm elements
            ModelPart::ElementsContainerType &rmesh_elements = p_vm_model_part->Elements();  

            // creating vm conditions
            ModelPart::ConditionsContainerType &rmesh_conditions = p_vm_model_part->Conditions(); 

            
            if(only_suf_param){
                for (int i = 0; i < (int)r_controlling_object.Conditions().size(); i++) {
                    ModelPart::ConditionsContainerType::iterator it = r_controlling_object.ConditionsBegin() + i;
                    Element::Pointer p_element = new HelmholtzSurfShapeElement(it->Id(), it->pGetGeometry(), p_vm_model_part_property);
                    rmesh_elements.push_back(p_element);
                }  
            }
            else{
                for (int i = 0; i < (int)root_model_part.Elements().size(); i++) {
                    ModelPart::ElementsContainerType::iterator it = root_model_part.ElementsBegin() + i;
                    Element::Pointer p_element = new HelmholtzBulkShapeElement(it->Id(), it->pGetGeometry(), p_vm_model_part_property);                
                    rmesh_elements.push_back(p_element);
                }            
                for (int i = 0; i < (int)r_controlling_object.Conditions().size(); i++) {
                    ModelPart::ConditionsContainerType::iterator it = r_controlling_object.ConditionsBegin() + i;
                    Condition::Pointer p_condition = new HelmholtzSurfShapeCondition(it->Id(), it->pGetGeometry(), p_vm_model_part_property);               
                    rmesh_conditions.push_back(p_condition);
                }    
                TetrahedralMeshOrientationCheck tetrahedralMeshOrientationCheck(*p_vm_model_part,false,TetrahedralMeshOrientationCheck::ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS);
                tetrahedralMeshOrientationCheck.Execute();                  
            }
        }

        // now add dofs
        for(int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];
            for(auto& node_i : mpVMModePart->Nodes())
            {
                node_i.AddDof(HELMHOLTZ_VARS_SHAPE_X);
                node_i.AddDof(HELMHOLTZ_VARS_SHAPE_Y);
                node_i.AddDof(HELMHOLTZ_VARS_SHAPE_Z);
            }
        }

        const auto& fixed_model_parts =  mTechniqueSettings["fixed_model_parts"];
        const auto& fixed_model_parts_X =  mTechniqueSettings["fixed_model_parts_X"];
        const auto& fixed_model_parts_Y =  mTechniqueSettings["fixed_model_parts_Y"];
        const auto& fixed_model_parts_Z =  mTechniqueSettings["fixed_model_parts_Z"];

        for(int i=0; i<fixed_model_parts.size();i++)
        {
            const auto& model_part = mrModel.GetModelPart(fixed_model_parts[i].GetString());
            for(auto& node_i : model_part.Nodes())
            {
                if(fixed_model_parts_X[i].GetBool())
                    node_i.Fix(HELMHOLTZ_VARS_SHAPE_X);
                if(fixed_model_parts_Y[i].GetBool())
                    node_i.Fix(HELMHOLTZ_VARS_SHAPE_Y);
                if(fixed_model_parts_Z[i].GetBool())
                    node_i.Fix(HELMHOLTZ_VARS_SHAPE_Z);
            }
        }
    }

    void AdjustFilterSizes()
    {
        // now adjust the filter size if its adaptive
        if(mTechniqueSettings["only_design_surface_parameterization"].GetBool() && (mTechniqueSettings["automatic_filter_size"].GetBool() || 
        mTechniqueSettings["adaptive_filter_size"].GetBool()) ){
            for(auto& control_obj : mControlSettings["controlling_objects"]){
                ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
                ModelPart& root_model_part = r_controlling_object.GetRootModelPart();
                std::string vm_model_part_name =  root_model_part.Name()+"_Implicit_VM_Part";
                ModelPart* mpVMModePart = &(root_model_part.GetSubModelPart(vm_model_part_name));
                Properties::Pointer p_vm_model_part_property;
                for(int i =0; i<mpVMModelParts.size(); i++)
                    if(mpVMModelParts[i]->Name()==mpVMModePart->Name())
                        p_vm_model_part_property = mpVMModelPartsProperties[i];

                double max_detJ = 0.0;
                double min_detJ = 1e9;
                for(auto& elem_i : mpVMModePart->Elements()){
                    const auto& r_geom = elem_i.GetGeometry();
                    const IntegrationMethod& integration_method = r_geom.GetDefaultIntegrationMethod();
                    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(integration_method);
                    const unsigned int NumGauss = integration_points.size();
                    Vector GaussPtsJDet = ZeroVector(NumGauss);
                    r_geom.DeterminantOfJacobian(GaussPtsJDet,integration_method);  
                    for(std::size_t i_point = 0; i_point<NumGauss; ++i_point){
                        if(GaussPtsJDet[i_point]>max_detJ)
                            max_detJ = GaussPtsJDet[i_point];
                        if(GaussPtsJDet[i_point]<min_detJ)
                            min_detJ = GaussPtsJDet[i_point];                  
                    }   
                }     

                double surface_filter_size = sqrt(std::pow(max_detJ/min_detJ, 0.5));
                p_vm_model_part_property->SetValue(HELMHOLTZ_SURF_RADIUS_SHAPE,surface_filter_size);

                KRATOS_INFO("ImplicitVertexMorphing:AdjustFilterSizes") << " Surface filter of "<<control_obj.GetString() <<" is adjusted to " << surface_filter_size << std::endl;
            }
        }
        else if(!mTechniqueSettings["only_design_surface_parameterization"].GetBool()){

            for(auto& control_obj : mControlSettings["controlling_objects"]){
                ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
                ModelPart& root_model_part = r_controlling_object.GetRootModelPart();
                std::string vm_model_part_name =  root_model_part.Name()+"_Implicit_VM_Part";
                ModelPart* mpVMModePart = &(root_model_part.GetSubModelPart(vm_model_part_name)); 
                Properties::Pointer p_vm_model_part_property;
                for(int i =0; i<mpVMModelParts.size(); i++)
                    if(mpVMModelParts[i]->Name()==mpVMModePart->Name())
                        p_vm_model_part_property = mpVMModelPartsProperties[i];           

                double max_detJ = 0.0;
                double min_detJ = 1e9;
                double total_vol = 0.0;
                double vol_int_radius = 0.0;
                for(auto& elem_i : mpVMModePart->Elements()){
                    const auto& r_geom = elem_i.GetGeometry(); 
                    total_vol += r_geom.Volume();
                    const IntegrationMethod& integration_method = r_geom.GetDefaultIntegrationMethod();
                    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(integration_method);
                    const unsigned int NumGauss = integration_points.size();
                    Vector GaussPtsJDet = ZeroVector(NumGauss);
                    r_geom.DeterminantOfJacobian(GaussPtsJDet,integration_method);  
                    for(std::size_t i_point = 0; i_point<NumGauss; ++i_point){
                        if(GaussPtsJDet[i_point]>max_detJ)
                            max_detJ = GaussPtsJDet[i_point];
                        if(GaussPtsJDet[i_point]<min_detJ)
                            min_detJ = GaussPtsJDet[i_point];                  
                    }   
                    
                }
                for(auto& elem_i : mpVMModePart->Elements()){
                    const auto& r_geom = elem_i.GetGeometry(); 
                    total_vol += r_geom.Volume();
                    const IntegrationMethod& integration_method = r_geom.GetDefaultIntegrationMethod();
                    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(integration_method);
                    const unsigned int NumGauss = integration_points.size();
                    Vector GaussPtsJDet = ZeroVector(NumGauss);
                    r_geom.DeterminantOfJacobian(GaussPtsJDet,integration_method);  

                    for(std::size_t i_point = 0; i_point<NumGauss; ++i_point){
                        vol_int_radius += GaussPtsJDet[i_point] * integration_points[i_point].Weight() * std::pow(1.0/GaussPtsJDet[i_point],1.5);
                    }
                }

                // std::cout<<"max_detJ: "<<max_detJ<<", min_detJ: "<<min_detJ<<std::endl; 

                double max_cond_par_elem_detJ = 0.0;
                double min_cond_par_elem_detJ = 1e9;
                double max_cond_detJ = 0.0;
                double min_cond_detJ = 1e9;            
                double total_cond_par_elem_vol = 0.0;
                double total_cond_area = 0.0;
                double min_cond_length = 1e9;
                double max_cond_length = 0.0;
                for(auto& cond_i : mpVMModePart->Conditions()){
                    const auto& r_geom = cond_i.GetGeometry(); 
                    double cond_area = r_geom.Area();
                    double cond_length = r_geom.Length();
                    total_cond_area += cond_area;

                    if(cond_length>max_cond_length)
                        max_cond_length = cond_length;
                    if(cond_length<min_cond_length)
                        min_cond_length = cond_length; 
         
                    const IntegrationMethod& cond_integration_method = r_geom.GetDefaultIntegrationMethod();
                    const GeometryType::IntegrationPointsArrayType& cond_integration_points = r_geom.IntegrationPoints(cond_integration_method);
                    const unsigned int cond_NumGauss = cond_integration_points.size();
                    Vector cond_GaussPtsJDet = ZeroVector(cond_NumGauss);
                    r_geom.DeterminantOfJacobian(cond_GaussPtsJDet,cond_integration_method);  
                    for(std::size_t cond_i_point = 0; cond_i_point<cond_NumGauss; ++cond_i_point){
                        if(cond_GaussPtsJDet[cond_i_point]>max_cond_detJ)
                            max_cond_detJ = cond_GaussPtsJDet[cond_i_point];
                        if(cond_GaussPtsJDet[cond_i_point]<min_cond_detJ)
                            min_cond_detJ = cond_GaussPtsJDet[cond_i_point];                  
                    } 

                    
                    auto& parentElement = cond_i.GetValue(NEIGHBOUR_ELEMENTS);
                    const auto& r_elem_geom = parentElement[0].GetGeometry();                                
                    const IntegrationMethod& integration_method = r_elem_geom.GetDefaultIntegrationMethod();
                    const GeometryType::IntegrationPointsArrayType& integration_points = r_elem_geom.IntegrationPoints(integration_method);
                    const unsigned int NumGauss = integration_points.size();
                    Vector GaussPtsJDet = ZeroVector(NumGauss);
                    r_elem_geom.DeterminantOfJacobian(GaussPtsJDet,integration_method);  
                    for(std::size_t i_point = 0; i_point<NumGauss; ++i_point){
                        if(GaussPtsJDet[i_point]>max_cond_par_elem_detJ)
                            max_cond_par_elem_detJ = GaussPtsJDet[i_point];
                        if(GaussPtsJDet[i_point]<min_cond_par_elem_detJ)
                            min_cond_par_elem_detJ = GaussPtsJDet[i_point];                  
                    }   
                }


                double surface_filter_size = p_vm_model_part_property->GetValue(HELMHOLTZ_SURF_RADIUS_SHAPE);
                double bulk_filter_size = (total_cond_area * surface_filter_size * surface_filter_size)/(36.0 * vol_int_radius);
                bulk_filter_size = std::pow(bulk_filter_size,2.0/3.0);
                

                p_vm_model_part_property->SetValue(HELMHOLTZ_BULK_RADIUS_SHAPE,bulk_filter_size);
                p_vm_model_part_property->SetValue(HELMHOLTZ_SURF_RADIUS_SHAPE,surface_filter_size);
                KRATOS_INFO("ImplicitVertexMorphing:AdjustFilterSizes") << " Surface and bulk filter sizes of "<<control_obj.GetString() <<" is adjusted to " << surface_filter_size <<
                " and "<< bulk_filter_size<<" respectively !"<<std::endl;
                // std::exit(0);
            }

        }      
    }  

    void ComputeInitialControlPoints()
    {
        for(int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];
            
            ProcessInfo &rCurrentProcessInfo = (mpVMModePart)->GetProcessInfo();
            
            SetVariableZero(mpVMModePart,HELMHOLTZ_SOURCE_SHAPE);
            SetVariableZero(mpVMModePart,CX);

            // // calculate (K+M) * X element wise, here we treat CX as an auxilary 
            rCurrentProcessInfo[COMPUTE_CONTROL_POINTS_SHAPE] = false;
            for(auto& node_i : mpVMModePart->Nodes())
            {
                array_3d& r_nodal_variable_hl_vars = node_i.FastGetSolutionStepValue(HELMHOLTZ_VARS_SHAPE);
                r_nodal_variable_hl_vars(0) = node_i.X0();
                r_nodal_variable_hl_vars(1) = node_i.Y0();
                r_nodal_variable_hl_vars(2) = node_i.Z0();
            }
            
            for(auto& elem_i : mpVMModePart->Elements())
            {
                VectorType rhs;
                MatrixType lhs;
                elem_i.Initialize(mpVMModePart->GetProcessInfo());
                elem_i.CalculateLocalSystem(lhs,rhs,mpVMModePart->GetProcessInfo());
                AddElementVariableValuesVector(elem_i,CX,rhs,-1.0);
            }        

            for(auto& cond_i : mpVMModePart->Conditions())
            {
                VectorType rhs;
                MatrixType lhs;
                cond_i.Initialize(mpVMModePart->GetProcessInfo());
                cond_i.CalculateLocalSystem(lhs,rhs,mpVMModePart->GetProcessInfo());
                AddConditionVariableValuesVector(cond_i,CX,rhs,-1.0);
            }

            
            // here we fill the RHS of M * S = (K+M) * X
            SetVariable1ToVarible2(mpVMModePart,CX,HELMHOLTZ_SOURCE_SHAPE);
            SetVariableZero(mpVMModePart,HELMHOLTZ_VARS_SHAPE);

            // apply BC on the RHS and unassign BC
            for(auto& node_i : mpVMModePart->Nodes())
            {
                array_3d& r_nodal_variable_hl_vars = node_i.FastGetSolutionStepValue(HELMHOLTZ_VARS_SHAPE);

                if(node_i.IsFixed(HELMHOLTZ_VARS_SHAPE_X))
                    r_nodal_variable_hl_vars(0) = node_i.X0();

                if(node_i.IsFixed(HELMHOLTZ_VARS_SHAPE_Y))
                    r_nodal_variable_hl_vars(1) = node_i.Y0();

                if(node_i.IsFixed(HELMHOLTZ_VARS_SHAPE_Z))
                    r_nodal_variable_hl_vars(2) = node_i.Z0();                
            }

            
            // here we compute S = M-1 * (K+M) * X
            rCurrentProcessInfo[COMPUTE_CONTROL_POINTS_SHAPE] = true;   
            mpStrategies[model_i]->Solve();
            // here we clear/reset the problem for the clearaty's sake
            mpStrategies[model_i]->GetStrategy()->Clear();
            SetVariable1ToVarible2(mpVMModePart,HELMHOLTZ_VARS_SHAPE,CX);
            rCurrentProcessInfo[COMPUTE_CONTROL_POINTS_SHAPE] = false;
        }

    }

    void GetElementVariableValuesVector(const Element& rElement,
                                        const Variable<array_3d> &rVariable,
                                        VectorType &rValues) const
    {
        const GeometryType &rgeom = rElement.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        SizeType dofs_per_node=3;

        const unsigned int local_size = num_nodes * dofs_per_node;

        if (rValues.size() != local_size)
            rValues.resize(local_size, false);

        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
            const array_3d& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable);    
            rValues[dofs_per_node*i_node+0] = r_nodal_variable[0];
            rValues[dofs_per_node*i_node+1] = r_nodal_variable[1];
            rValues[dofs_per_node*i_node+2] = r_nodal_variable[2];
            if(dofs_per_node==6){
                rValues[dofs_per_node*i_node+3] = 0;
                rValues[dofs_per_node*i_node+4] = 0;
                rValues[dofs_per_node*i_node+5] = 0;                  
            }
          
        }
    }
    void GetConditionVariableValuesVector(const Condition& rCondition,
                                        const Variable<array_3d> &rVariable,
                                        VectorType &rValues) const
    {
        const GeometryType &rgeom = rCondition.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const unsigned int dimension = rCondition.GetGeometry().WorkingSpaceDimension();
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
                                        const VectorType &rValues,
                                        const double& rWeight = 1.0
                                        ) 
    {
        GeometryType &rgeom = rElement.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();

        SizeType dofs_per_node=3;

        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
            array_3d& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable);
            r_nodal_variable[0] += rWeight * rValues[dofs_per_node*i_node+0];
            r_nodal_variable[1] += rWeight * rValues[dofs_per_node*i_node+1];
            r_nodal_variable[2] += rWeight * rValues[dofs_per_node*i_node+2];
        }
    }

    void AddConditionVariableValuesVector(Condition& rCondition,
                                        const Variable<array_3d> &rVariable,
                                        const VectorType &rValues,
                                        const double& rWeight = 1.0
                                        ) 
    {
        GeometryType &rgeom = rCondition.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const unsigned int dimension = rCondition.GetGeometry().WorkingSpaceDimension();

        if (dimension == 2) {
            SizeType index = 0;
            for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
                array_3d& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable); 
                r_nodal_variable[0] += rWeight * rValues[index++];
                r_nodal_variable[1] += rWeight * rValues[index++];
            }
        } else if (dimension == 3) {
            SizeType index = 0;
            for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
                array_3d& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable);
                r_nodal_variable[0] += rWeight * rValues[index++];
                r_nodal_variable[1] += rWeight * rValues[index++];
                r_nodal_variable[2] += rWeight * rValues[index++];
            }
        }
    }    
    
    void SetVariableZero(ModelPart* mpVMModePart, const Variable<array_3d> &rVariable) 
    {
        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_variable = node_i.FastGetSolutionStepValue(rVariable);
            r_nodal_variable[0] = 0.0;
            r_nodal_variable[1] = 0.0;
            r_nodal_variable[2] = 0.0;                    
        }
    }

    void SetVariable1ToVarible2(ModelPart* mpVMModePart,const Variable<array_3d> &rVariable1,const Variable<array_3d> &rVariable2) 
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
