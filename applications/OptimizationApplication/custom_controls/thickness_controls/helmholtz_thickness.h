//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

#ifndef HELMHOLTZ_THICKNESS_H
#define HELMHOLTZ_THICKNESS_H

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
#include "custom_controls/thickness_controls/thickness_control.h"
#include "custom_elements/helmholtz_surf_thickness_element.h"
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

class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzThickness : public ThicknessControl
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

    /// Pointer definition of HelmholtzThickness
    KRATOS_CLASS_POINTER_DEFINITION(HelmholtzThickness);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HelmholtzThickness( std::string ControlName, Model& rModel, std::vector<LinearSolverType::Pointer>& rLinearSolvers, Parameters ControlSettings )
        :  ThicknessControl(ControlName,rModel,ControlSettings){
            for(long unsigned int lin_i=0;lin_i<rLinearSolvers.size();lin_i++)
                rLinearSystemSolvers.push_back(rLinearSolvers[lin_i]);
            mTechniqueSettings = ControlSettings["technique_settings"];
        }

    /// Destructor.
    virtual ~HelmholtzThickness()
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
        KRATOS_INFO("HelmholtzThickness:Initialize ") << "Starting initialization of thickness control "<<mControlName<<" ..." << std::endl;

        CreateHelmholtzThicknessModelParts();

        CalculateNodeNeighbourCount();
        
        for(long unsigned int model_i=0;model_i<mpVMModelParts.size();model_i++){
            StrategyType* mpStrategy = new StrategyType (*mpVMModelParts[model_i],rLinearSystemSolvers[model_i]);            
            mpStrategy->Initialize();
            mpStrategies.push_back(mpStrategy);
        }

        physical_thicknesses =  mTechniqueSettings["physical_thicknesses"].GetVector();
        initial_thickness = mTechniqueSettings["initial_thickness"].GetDouble();
        SIMP_pow_fac = mTechniqueSettings["SIMP_power_fac"].GetInt();
        beta = mTechniqueSettings["beta_settings"]["initial_value"].GetDouble();
        adaptive_beta = mTechniqueSettings["beta_settings"]["adaptive"].GetBool();
        beta_fac = mTechniqueSettings["beta_settings"]["increase_fac"].GetDouble();
        max_beta = mTechniqueSettings["beta_settings"]["max_value"].GetDouble();
        beta_update_period = mTechniqueSettings["beta_settings"]["update_period"].GetInt();

        filtered_thicknesses.resize(physical_thicknesses.size());        
        for(long unsigned int i=0;i<physical_thicknesses.size();i++)
            filtered_thicknesses[i] = i;

        double initial_filtered_thickness = ProjectBackward(initial_thickness,filtered_thicknesses,physical_thicknesses,beta);

        for(long unsigned int model_i=0;model_i<mpVMModelParts.size();model_i++){
            SetVariable(mpVMModelParts[model_i],CT,initial_filtered_thickness); 
            SetVariable(mpVMModelParts[model_i],FT,initial_filtered_thickness); 
            SetVariable(mpVMModelParts[model_i],PT,initial_thickness);
        } 

        const auto& fixed_model_parts =  mTechniqueSettings["fixed_model_parts"];
        const auto& fixed_model_parts_thicknesses = mTechniqueSettings["fixed_model_parts_thicknesses"].GetVector();

        for(long unsigned int i=0; i<fixed_model_parts.size();i++)
        {
            const auto& model_part = mrModel.GetModelPart(fixed_model_parts[i].GetString());
            auto model_part_phyisical_thick = fixed_model_parts_thicknesses[i];
            double model_part_filtered_thick = ProjectBackward(model_part_phyisical_thick,filtered_thicknesses,physical_thicknesses,beta);
            #pragma omp parallel for
            for(auto& node_i : model_part.Nodes()){
                auto& control_thickness = node_i.FastGetSolutionStepValue(CT);
                auto& filtered_thickness = node_i.FastGetSolutionStepValue(FT);
                auto& physical_thickness = node_i.FastGetSolutionStepValue(PT);
                control_thickness = model_part_filtered_thick;
                filtered_thickness = model_part_filtered_thick;
                physical_thickness = model_part_phyisical_thick;
            }                
        }

        KRATOS_INFO("HelmholtzThickness:Initialize") << "Finished initialization of thickness control "<<mControlName<<" in " << timer.ElapsedSeconds() << " s." << std::endl;

    };
    // --------------------------------------------------------------------------
    void Update() override {
        opt_itr++;

        if(adaptive_beta){
            if (opt_itr % beta_update_period == 0 && beta <max_beta){
                beta *= beta_fac;
                if(beta>max_beta)
                    beta = max_beta;                  
                KRATOS_INFO("HelmholtzThickness:Update") << "beta is updated to " <<beta<< std::endl;
            }                          
        }

        ComputeFilteredThickness();
        ComputePhyiscalThickness();
    };  
    // --------------------------------------------------------------------------
    void MapControlUpdate(const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override{};
    // --------------------------------------------------------------------------
    void MapFirstDerivative(const Variable<double> &rDerivativeVariable, const Variable<double> &rMappedDerivativeVariable) override{


        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("HelmholtzThickness:MapFirstDerivative") << "Starting mapping of " << rDerivativeVariable.Name() << "..." << std::endl;

        for(long unsigned int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];
            SetVariable1ToVarible2(mpVMModePart,rDerivativeVariable,HELMHOLTZ_SOURCE_THICKNESS);
            SetVariable(mpVMModePart,HELMHOLTZ_VAR_THICKNESS,0.0);

            //now solve 
            mpStrategies[model_i]->Solve();
            SetVariable1ToVarible2(mpVMModePart,HELMHOLTZ_VAR_THICKNESS,rMappedDerivativeVariable);
        }
    
        KRATOS_INFO("HelmholtzThickness:MapFirstDerivative") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;

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
        return "HelmholtzThickness";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "HelmholtzThickness";
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
    double beta;
    int SIMP_pow_fac;
    bool adaptive_beta;
    double beta_fac;
    double max_beta;
    int beta_update_period;
    int opt_itr = 0;
    Vector physical_thicknesses;
    Vector filtered_thicknesses;
    double initial_thickness; 
    
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
        for(long unsigned int model_i =0;model_i<mpVMModelParts.size();model_i++)
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

    void CreateHelmholtzThicknessModelParts()
    {
        // creating vm model nodes and variables
        for(auto& control_obj : mControlSettings["controlling_objects"]){
            ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
            ModelPart& root_model_part = r_controlling_object.GetRootModelPart();
            std::string vm_model_part_name =  root_model_part.Name()+"_Helmholtz_Thickness_Part";
            ModelPart* p_vm_model_part;
            Properties::Pointer p_vm_model_part_property;

            if (root_model_part.HasSubModelPart(vm_model_part_name)){
                p_vm_model_part = &(root_model_part.GetSubModelPart(vm_model_part_name));
                for(long unsigned int i =0; i<mpVMModelParts.size(); i++)
                    if(mpVMModelParts[i]->Name()==p_vm_model_part->Name())
                        p_vm_model_part_property = mpVMModelPartsProperties[i];
            }
            else{
                p_vm_model_part = &(root_model_part.CreateSubModelPart(vm_model_part_name));
                p_vm_model_part_property = p_vm_model_part->CreateNewProperties(root_model_part.NumberOfProperties()+1);
                mpVMModelPartsProperties.push_back(p_vm_model_part_property);
                mpVMModelParts.push_back(p_vm_model_part);
            }

            p_vm_model_part_property->SetValue(HELMHOLTZ_RADIUS_THICKNESS,mTechniqueSettings["filter_radius"].GetDouble());

            for(auto& node : r_controlling_object.Nodes())
                p_vm_model_part->AddNode(&node);

            // creating elements
            ModelPart::ElementsContainerType &rmesh_elements = p_vm_model_part->Elements();   

            //check if the controlling model part has elements which have thickness value
            if(!(r_controlling_object.Elements().size()>0))
                KRATOS_ERROR << "HelmholtzThickness:CreateModelParts : controlling model part " <<control_obj.GetString()<<" does not have elements"<<std::endl;


            for (int i = 0; i < (int)r_controlling_object.Elements().size(); i++) {
                ModelPart::ElementsContainerType::iterator it = r_controlling_object.ElementsBegin() + i;
                const Properties& elem_i_prop = it->GetProperties();
                Properties::Pointer elem_i_new_prop = r_controlling_object.CreateNewProperties(r_controlling_object.NumberOfProperties()+1);
                *elem_i_new_prop = elem_i_prop;
                it->SetProperties(elem_i_new_prop);
                Element::Pointer p_element = new HelmholtzSurfThicknessElement(it->Id(), it->pGetGeometry(), p_vm_model_part_property);
                rmesh_elements.push_back(p_element);
            }   
        }

        // now add dofs
        for(long unsigned int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];
            for(auto& node_i : mpVMModePart->Nodes())
            {
                node_i.AddDof(HELMHOLTZ_VAR_THICKNESS);
            }
        }

        // now apply dirichlet BC
        const auto& fixed_model_parts =  mTechniqueSettings["fixed_model_parts"];

        for(long unsigned int i=0; i<fixed_model_parts.size();i++)
        {
            const auto& model_part = mrModel.GetModelPart(fixed_model_parts[i].GetString());
            for(auto& node_i : model_part.Nodes())
                node_i.Fix(HELMHOLTZ_VAR_THICKNESS);
        }   


    }

    void ComputeFilteredThickness(){   

        for(long unsigned int model_i=0;model_i<mpVMModelParts.size();model_i++){

            //first update control thickness
            AddVariable1ToVarible2(mpVMModelParts[model_i],D_CT,CT);

            //now filter nodal control thickness 
            //first we need to multiply with mass matrix 
            SetVariable(mpVMModelParts[model_i],HELMHOLTZ_VAR_THICKNESS,0.0);
            SetVariable(mpVMModelParts[model_i],HELMHOLTZ_SOURCE_THICKNESS,0.0);
            //now we need to multiply with the mass matrix
            #pragma omp parallel for 
            for(auto& elem_i : mpVMModelParts[model_i]->Elements())
            {
                VectorType origin_values;
                GetElementVariableValuesVector(elem_i,CT,origin_values);
                MatrixType mass_matrix;
                elem_i.Calculate(HELMHOLTZ_MASS_MATRIX,mass_matrix,mpVMModelParts[model_i]->GetProcessInfo());           
                VectorType int_vals = prod(mass_matrix,origin_values);
                AddElementVariableValuesVector(elem_i,HELMHOLTZ_SOURCE_THICKNESS,int_vals);
            }

            //now we need to apply BCs
            #pragma omp parallel for 
            for(auto& node_i : mpVMModelParts[model_i]->Nodes())
            {
                auto& r_nodal_variable_hl_var = node_i.FastGetSolutionStepValue(HELMHOLTZ_VAR_THICKNESS);
                if(node_i.IsFixed(HELMHOLTZ_VAR_THICKNESS))
                    r_nodal_variable_hl_var = node_i.FastGetSolutionStepValue(FT);           
            }
            mpStrategies[model_i]->Solve();
            SetVariable1ToVarible2(mpVMModelParts[model_i],HELMHOLTZ_VAR_THICKNESS,FT);
        }        
    }     

    void ComputePhyiscalThickness(){
        for(long unsigned int model_i=0;model_i<mpVMModelParts.size();model_i++){
            //now do the projection and then set the PT
            #pragma omp parallel for
            for(auto& node_i : mpVMModelParts[model_i]->Nodes()){
                const auto& filtered_thick = node_i.FastGetSolutionStepValue(FT);
                auto& thickness = node_i.FastGetSolutionStepValue(PT);
                auto& p_thickness = node_i.FastGetSolutionStepValue(PPT);
                auto& thickness_der = node_i.FastGetSolutionStepValue(D_PT_D_FT);
                auto& p_thickness_der = node_i.FastGetSolutionStepValue(D_PPT_D_FT);
                thickness = ProjectForward(filtered_thick,filtered_thicknesses,physical_thicknesses,beta,1);
                p_thickness = ProjectForward(filtered_thick,filtered_thicknesses,physical_thicknesses,beta,SIMP_pow_fac);
                thickness_der = ProjectionDerivative(filtered_thick,filtered_thicknesses,physical_thicknesses,beta,1);
                p_thickness_der = ProjectionDerivative(filtered_thick,filtered_thicknesses,physical_thicknesses,beta,SIMP_pow_fac);
            }
        }

        // update elements' t
        for(auto& control_obj : mControlSettings["controlling_objects"]){
            ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
            #pragma omp parallel for
            for (int i = 0; i < (int)r_controlling_object.Elements().size(); i++) {
                ModelPart::ElementsContainerType::iterator it = r_controlling_object.ElementsBegin() + i;
                double elem_i_thicknes = 0.0;
                double elem_i_p_thicknes = 0.0;
                for(unsigned int node_element = 0; node_element<it->GetGeometry().size(); node_element++){
                    elem_i_thicknes += it->GetGeometry()[node_element].FastGetSolutionStepValue(PT);
                    elem_i_p_thicknes += it->GetGeometry()[node_element].FastGetSolutionStepValue(PPT);
                }
                    
                elem_i_thicknes /= it->GetGeometry().size();
                elem_i_p_thicknes /= it->GetGeometry().size();

                it->GetProperties().SetValue(THICKNESS,elem_i_p_thicknes);
                it->GetProperties().SetValue(PT,elem_i_thicknes);
            }
        }
    }


    double ProjectForward(double x,Vector x_limits,Vector y_limits,double beta, int penal_fac = 1){

        double x1=0,x2=0,y1=0,y2=0;
        if(x>=x_limits[x_limits.size()-1]){
            x1=x_limits[x_limits.size()-2];
            x2=x_limits[x_limits.size()-1];
            y1=y_limits[y_limits.size()-2];
            y2=y_limits[y_limits.size()-1];
        }
        else if(x<=x_limits[0]){
            x1=x_limits[0];
            x2=x_limits[1];
            y1=y_limits[0];
            y2=y_limits[1];
        }
        else{
            for(long unsigned int i=0;i<x_limits.size()-1;i++)
                if((x>=x_limits[i]) && (x<=x_limits[i+1]))
                {
                    y1 = y_limits[i];
                    y2 = y_limits[i+1];
                    x1 = x_limits[i];
                    x2 = x_limits[i+1];
                    break;
                }            
        }        
        
        double pow_val = -2.0*beta*(x-(x1+x2)/2);

        if(pow_val>700)
            pow_val = 700;
        if(pow_val<-700) 
            pow_val = -700;        

        return (y2-y1)/(std::pow(1+std::exp(pow_val),penal_fac)) + y1;
    }


    double ProjectBackward(double y,Vector x_limits,Vector y_limits,double beta){
        
        double x = 0;
        if(y>=y_limits[y_limits.size()-1])
            x = x_limits[y_limits.size()-1];
        else if(y<=y_limits[0])
            x = x_limits[0];
        else{
            for(long unsigned int i=0;i<y_limits.size()-1;i++)
                if((y>y_limits[i]) && (y<y_limits[i+1]))
                {
                    double y1 = y_limits[i];
                    double y2 = y_limits[i+1];
                    double x1 = x_limits[i];
                    double x2 = x_limits[i+1];
                    x = ((x2+x1)/2.0) + (1.0/(-2.0*beta)) * std::log(((y2-y1)/(y-y1))-1);
                    break;
                } 
                else if(y==y_limits[i]){
                    x = x_limits[i];
                    break;
                }
                else if(y==y_limits[i+1]){
                    x = x_limits[i+1];
                    break;
                }                           
        }
        return x;
    }

    double ProjectionDerivative(double x,Vector x_limits,Vector y_limits,double beta,int penal_fac = 1){

        double x1=0.0,x2=0.0,y1=0.0,y2=0.0;
        if(x>=x_limits[x_limits.size()-1]){
            x1=x_limits[x_limits.size()-2];
            x2=x_limits[x_limits.size()-1];
            y1=y_limits[y_limits.size()-2];
            y2=y_limits[y_limits.size()-1];
        }
        else if(x<=x_limits[0]){
            x1=x_limits[0];
            x2=x_limits[1];
            y1=y_limits[0];
            y2=y_limits[1];
        }
        else{
            for(long unsigned int i=0;i<x_limits.size()-1;i++)
                if((x>=x_limits[i]) && (x<=x_limits[i+1]))
                {
                    y1 = y_limits[i];
                    y2 = y_limits[i+1];
                    x1 = x_limits[i];
                    x2 = x_limits[i+1];
                    break;
                }            
        }

        double pow_val = -2.0*beta*(x-(x1+x2)/2);

        if(pow_val>700)
            pow_val = 700;
        if(pow_val<-700) 
            pow_val = -700;

        double dydx = (y2-y1) * (1.0/std::pow(1+std::exp(pow_val),penal_fac+1)) * penal_fac * 2.0 * beta * std::exp(pow_val);

        // if (y2<y1)
        //     dydx *=-1;

        return dydx;

    }  

    void GetElementVariableValuesVector(const Element& rElement,
                                        const Variable<double> &rVariable,
                                        VectorType &rValues) const
    {
        const GeometryType &rgeom = rElement.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();

        const unsigned int local_size = num_nodes;

        if (rValues.size() != local_size)
            rValues.resize(local_size, false);

        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
            const auto& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable);    
            rValues[i_node] = r_nodal_variable; 
        }
    }
    void GetConditionVariableValuesVector(const Condition& rCondition,
                                        const Variable<double> &rVariable,
                                        VectorType &rValues) const
    {
        const GeometryType &rgeom = rCondition.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const unsigned int local_size = num_nodes;

        if (rValues.size() != local_size)
            rValues.resize(local_size, false);

        SizeType index = 0;
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
            const auto& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable);    
            rValues[index++] = r_nodal_variable;
        }
    }    
    void AddElementVariableValuesVector(Element& rElement,
                                        const Variable<double> &rVariable,
                                        const VectorType &rValues,
                                        const double& rWeight = 1.0
                                        ) 
    {
        GeometryType &rgeom = rElement.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();

        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
            #pragma omp atomic
            rgeom[i_node].FastGetSolutionStepValue(rVariable) += (rWeight * rValues[i_node]);
        }
    }

    void AddConditionVariableValuesVector(Condition& rCondition,
                                        const Variable<double> &rVariable,
                                        const VectorType &rValues,
                                        const double& rWeight = 1.0
                                        ) 
    {
        GeometryType &rgeom = rCondition.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();

        SizeType index = 0;
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
            auto& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable);
            r_nodal_variable += rWeight * rValues[index++];
        }
    }    
    
    void SetVariable(ModelPart* mpVMModePart, const Variable<double> &rVariable, const double value) 
    {
        #pragma omp parallel for
        for(auto& node_i : mpVMModePart->Nodes())
            node_i.FastGetSolutionStepValue(rVariable) = value;
    }

    void SetVariable1ToVarible2(ModelPart* mpVMModePart,const Variable<double> &rVariable1,const Variable<double> &rVariable2) 
    {
        #pragma omp parallel for
        for(auto& node_i : mpVMModePart->Nodes())
            node_i.FastGetSolutionStepValue(rVariable2) = node_i.FastGetSolutionStepValue(rVariable1);
    }    

    void AddVariable1ToVarible2(ModelPart* mpVMModePart,const Variable<double> &rVariable1,const Variable<double> &rVariable2) 
    {
        #pragma omp parallel for
        for(auto& node_i : mpVMModePart->Nodes())
            node_i.FastGetSolutionStepValue(rVariable2) += node_i.FastGetSolutionStepValue(rVariable1);
    }    

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//      HelmholtzThickness& operator=(HelmholtzThickness const& rOther);

    /// Copy constructor.
//      HelmholtzThickness(HelmholtzThickness const& rOther);


    ///@}

}; // Class HelmholtzThickness

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // HELMHOLTZ_THICKNESS_H
