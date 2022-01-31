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
#include "processes/tetrahedral_mesh_orientation_check.h"
#include "utilities/builtin_timer.h"
#include "spaces/ublas_space.h"
#include "mapper_base.h"
#include "custom_conditions/helmholtz_condition.h"
#include "custom_elements/helmholtz_element.h"
#include "custom_elements/helmholtz_surf_element.h"
#include "custom_elements/helmholtz_surf_prism_element.h"
#include "custom_elements/helmholtz_vec_element.h"
#include "custom_strategies/strategies/helmholtz_strategy.h"
#include "custom_strategies/strategies/helmholtz_vec_strategy.h"
#include "containers/model.h"
#include "linear_solvers/linear_solver.h"
#include "input_output/vtk_output.h"

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
        ModelPart::ElementsContainerType &rmesh_elements = mpVMModePart->Elements();  

        // creating vm conditions
        ModelPart::ConditionsContainerType &rmesh_conditions = mpVMModePart->Conditions();         

        bool only_suf_param =  mMapperSettings["only_design_surface_parameterization"].GetBool();
        ModelPart& design_surface_sub_model_part = mrModelPart.GetSubModelPart(mMapperSettings["design_surface_sub_model_part_name"].GetString());

        // create a new property for the vm surface       
        p_vm_surf_property = mpVMModePart->CreateNewProperties(mpVMModePart->NumberOfProperties()+1);        
        p_vm_surf_property->SetValue(HELMHOLTZ_RADIUS,mMapperSettings["surface_filter_radius"].GetDouble());   

        // create a new property for the vm elements       
        p_vm_bulk_property = mpVMModePart->CreateNewProperties(mpVMModePart->NumberOfProperties()+1);        
        p_vm_bulk_property->SetValue(HELMHOLTZ_RADIUS,mMapperSettings["bulk_filter_radius"].GetDouble());              

        if(only_suf_param){
            p_vm_surf_property->SetValue(HELMHOLTZ_POISSON_RATIO,mMapperSettings["poisson_ratio"].GetDouble()); 
            for (int i = 0; i < (int)design_surface_sub_model_part.Conditions().size(); i++) {
                ModelPart::ConditionsContainerType::iterator it = design_surface_sub_model_part.ConditionsBegin() + i;
                Element::Pointer p_element;
                p_element = new HelmholtzSurfPrismElement(it->Id(), it->pGetGeometry(), p_vm_surf_property);               
                rmesh_elements.push_back(p_element);
            }  
        }
        else{
            for (int i = 0; i < (int)mrModelPart.Elements().size(); i++) {
                ModelPart::ElementsContainerType::iterator it = mrModelPart.ElementsBegin() + i;
                Element::Pointer p_element;
                if (element_type.compare("helmholtz_element") == 0)
                    p_element = new HelmholtzElement(it->Id(), it->pGetGeometry(), p_vm_bulk_property);
                else
                    p_element = new HelmholtzVecElement(it->Id(), it->pGetGeometry(), p_vm_bulk_property);                
                rmesh_elements.push_back(p_element);
            }            
            for (int i = 0; i < (int)design_surface_sub_model_part.Conditions().size(); i++) {
                ModelPart::ConditionsContainerType::iterator it = design_surface_sub_model_part.ConditionsBegin() + i;
                Condition::Pointer p_condition;
                p_condition = new HelmholtzCondition(it->Id(), it->pGetGeometry(), p_vm_surf_property);               
                rmesh_conditions.push_back(p_condition);
            }    

            TetrahedralMeshOrientationCheck tetrahedralMeshOrientationCheck(*mpVMModePart,false,TetrahedralMeshOrientationCheck::ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS);
            tetrahedralMeshOrientationCheck.Execute();                       
        }

        // calculate number of neighbour elements for each node.
        CalculateNodeNeighbourCount();

        mpHelmholtzStrategy = new HelmholtzVecStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType> (*mpVMModePart,mpLinearSystemSolver);            
        
        mpHelmholtzStrategy->Initialize();

        // now adjust the filter size if its adaptive
        if(mMapperSettings["only_design_surface_parameterization"].GetBool() && (mMapperSettings["automatic_filter_size"].GetBool() || 
        mMapperSettings["adaptive_filter_size"].GetBool()) ){

            double max_length = 0;
            for(auto& elem_i : mpVMModePart->Elements())
                if(elem_i.GetGeometry().Length()>max_length)
                    max_length = elem_i.GetGeometry().Length();

            double surface_filter_size = max_length * 4;
            p_vm_surf_property->SetValue(HELMHOLTZ_RADIUS,surface_filter_size); 
            KRATOS_INFO("ShapeOpt") << " surface filter size is adjusted to " << surface_filter_size << std::endl;
        }
        else if(!mMapperSettings["only_design_surface_parameterization"].GetBool() && element_type.compare("helmholtz_vec_element") == 0){
            double max_elem_weight = 0;
            for(auto& elem_i : mpVMModePart->Elements())
            {
                double max_gp_weight_factor = 0;
                elem_i.Calculate(HELMHOLTZ_POISSON_RATIO,max_gp_weight_factor,mpVMModePart->GetProcessInfo());
                if(max_gp_weight_factor>max_elem_weight)
                    max_elem_weight = max_gp_weight_factor;
            }
            double surface_filter_size = 1.1 * std::sqrt(max_elem_weight);        
            p_vm_surf_property->SetValue(HELMHOLTZ_RADIUS,surface_filter_size);
            KRATOS_INFO("ShapeOpt") << " surface filter size is adjusted to " << surface_filter_size << std::endl;
        }        


        // calculate the initial control points
        ProcessInfo &rCurrentProcessInfo = (mpVMModePart)->GetProcessInfo();
        
        SetVariableZero(HELMHOLTZ_SOURCE);
        SetVariableZero(CONTROL_POINT);

        // // calculate (K+M) * X element wise, here we treat CONTROL_POINT as an auxilary 
        rCurrentProcessInfo[COMPUTE_CONTROL_POINTS] = false;
        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_variable_hl_vars = node_i.FastGetSolutionStepValue(HELMHOLTZ_VARS);
            r_nodal_variable_hl_vars(0) = node_i.X0();
            r_nodal_variable_hl_vars(1) = node_i.Y0();
            r_nodal_variable_hl_vars(2) = node_i.Z0();
        }
        
        for(auto& elem_i : mpVMModePart->Elements())
        {
            VectorType rhs;
            MatrixType lhs;
            elem_i.CalculateLocalSystem(lhs,rhs,mpVMModePart->GetProcessInfo());
            AddElementVariableValuesVector(elem_i,CONTROL_POINT,rhs,-1.0);
        }

        for(auto& cond_i : mpVMModePart->Conditions())
        {
            VectorType rhs;
            MatrixType lhs;
            cond_i.CalculateLocalSystem(lhs,rhs,mpVMModePart->GetProcessInfo());
            AddConditionVariableValuesVector(cond_i,CONTROL_POINT,rhs,-1.0);
        }

        // here we fill the RHS of M * S = (K+M) * X
        SetVariable1ToVarible2(CONTROL_POINT,HELMHOLTZ_SOURCE);
        SetVariableZero(HELMHOLTZ_VARS);

        // apply BC on the RHS and unassign BC
        for(auto& node_i : mpVMModePart->Nodes())
        {
            array_3d& r_nodal_variable_hl_vars = node_i.FastGetSolutionStepValue(HELMHOLTZ_VARS);

            if(node_i.IsFixed(HELMHOLTZ_VARS_X))
                r_nodal_variable_hl_vars(0) = node_i.X0();

            if(node_i.IsFixed(HELMHOLTZ_VARS_Y))
                r_nodal_variable_hl_vars(1) = node_i.Y0();

            if(node_i.IsFixed(HELMHOLTZ_VARS_Z))
                r_nodal_variable_hl_vars(2) = node_i.Z0();                
        }

        

        // here we compute S = M-1 * (K+M) * X
        rCurrentProcessInfo[COMPUTE_CONTROL_POINTS] = true;    
        mpHelmholtzStrategy->Solve();
        // here we clear/reset the problem for the clearaty's sake
        mpHelmholtzStrategy->GetStrategy()->Clear();
        SetVariable1ToVarible2(HELMHOLTZ_VARS,CONTROL_POINT);
        rCurrentProcessInfo[COMPUTE_CONTROL_POINTS] = false;
               
        //now export the linear system if needed
        if(mMapperSettings["export_linear_system"].GetBool()){            

            //set coords to the HELMHOLTZ_SOURCE to be exported for out processing
            for(auto& node_i : mpVMModePart->Nodes())
            {
                array_3d& r_nodal_variable_hl_source = node_i.FastGetSolutionStepValue(HELMHOLTZ_SOURCE);
                r_nodal_variable_hl_source(0) = node_i.X0();
                r_nodal_variable_hl_source(1) = node_i.Y0();
                r_nodal_variable_hl_source(2) = node_i.Z0();
            }

            SetVariableZero(HELMHOLTZ_VARS);
            mpHelmholtzStrategy->ExportSystem();
            SetVariableZero(HELMHOLTZ_VARS);
            SetVariableZero(HELMHOLTZ_SOURCE);
        }
       
        
        mIsMappingInitialized = true;

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
            GetElementVariableValuesVector(elem_i,rOriginVariable,origin_values);
            MatrixType mass_matrix;
            elem_i.Calculate(HELMHOLTZ_MASS_MATRIX,mass_matrix,mpVMModePart->GetProcessInfo());            
            VectorType int_vals = prod(mass_matrix,origin_values);
            AddElementVariableValuesVector(elem_i,HELMHOLTZ_SOURCE,int_vals);
        }

        for(auto& cond_i : mpVMModePart->Conditions())
        {
            VectorType origin_values;
            GetConditionVariableValuesVector(cond_i,rOriginVariable,origin_values);
            MatrixType mass_matrix;
            cond_i.Calculate(HELMHOLTZ_MASS_MATRIX,mass_matrix,mpVMModePart->GetProcessInfo());            
            VectorType int_vals = prod(mass_matrix,origin_values);
            AddConditionVariableValuesVector(cond_i,HELMHOLTZ_SOURCE,int_vals);
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
        // SetVariableZero(rOriginVariable);
        // //now we need to multiply with the mass matrix with results to be consistence.
        // for(auto& elem_i : mpVMModePart->Elements())
        // {
        //     VectorType helmholtz_values;
        //     elem_i.GetValuesVector(helmholtz_values);
        //     MatrixType mass_matrix;
        //     elem_i.Calculate(HELMHOLTZ_MASS_MATRIX,mass_matrix,mpVMModePart->GetProcessInfo());
        //     VectorType int_vals = prod(mass_matrix,helmholtz_values);
        //     AddElementVariableValuesVector(elem_i,rOriginVariable,int_vals);
        // }

        // for(auto& cond_i : mpVMModePart->Conditions())
        // {
        //     VectorType helmholtz_values;
        //     cond_i.GetValuesVector(helmholtz_values);
        //     MatrixType mass_matrix;
        //     cond_i.Calculate(HELMHOLTZ_MASS_MATRIX,mass_matrix,mpVMModePart->GetProcessInfo());
        //     VectorType int_vals = prod(mass_matrix,helmholtz_values);
        //     AddConditionVariableValuesVector(cond_i,rOriginVariable,int_vals);
        // }

        SetVariable1ToVarible2(HELMHOLTZ_VARS,rOriginVariable);
    
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
        TetrahedralMeshOrientationCheck tetrahedralMeshOrientationCheck(*mpVMModePart,false);
        tetrahedralMeshOrientationCheck.Execute();

        // now adjust the filter size if it is needed
        if(mMapperSettings["only_design_surface_parameterization"].GetBool() && mMapperSettings["adaptive_filter_size"].GetBool() && 
        !mMapperSettings["formulate_on_the_undeformed_configuration"].GetBool()){

            double max_length = 0;
            for(auto& elem_i : mpVMModePart->Elements())
                if(elem_i.GetGeometry().Length()>max_length)
                    max_length = elem_i.GetGeometry().Length();

            double surface_filter_size = max_length * 4;
            p_vm_surf_property->SetValue(HELMHOLTZ_RADIUS,surface_filter_size); 
            KRATOS_INFO("ShapeOpt") << " surface filter size is adjusted to " << surface_filter_size << std::endl;

        }
        // else if(!mMapperSettings["only_design_surface_parameterization"].GetBool() && 
        //         mMapperSettings["element_type"].GetString().compare("helmholtz_vec_element") == 0 &&
        //         !mMapperSettings["formulate_on_the_undeformed_configuration"].GetBool()){
        //         double max_elem_weight = 0;
        //         for(auto& elem_i : mpVMModePart->Elements())
        //         {
        //             double max_gp_weight_factor = 0;
        //             elem_i.Calculate(HELMHOLTZ_POISSON_RATIO,max_gp_weight_factor,mpVMModePart->GetProcessInfo());
        //             if(max_gp_weight_factor>max_elem_weight)
        //                 max_elem_weight = max_gp_weight_factor;
        //         }
        //         double surface_filter_size = 1.1 * std::sqrt(max_elem_weight);        
        //         p_vm_surf_property->SetValue(HELMHOLTZ_RADIUS,surface_filter_size);
        //         KRATOS_INFO("ShapeOpt") << " surface filter size is adjusted to " << surface_filter_size << std::endl;
        // }  

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
    HelmholtzVecStrategy<SparseSpaceType, LocalSpaceType,LinearSolverType>* mpHelmholtzStrategy;
    Properties::Pointer p_vm_bulk_property;
    Properties::Pointer p_vm_surf_property;

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
        const unsigned int dimension = rElement.GetGeometry().WorkingSpaceDimension();

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
