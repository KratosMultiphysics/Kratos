// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef HELMHOLTZ_MATERIAL_H
#define HELMHOLTZ_MATERIAL_H

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
#include "custom_controls/material_controls/material_control.h"
#include "custom_elements/helmholtz_bulk_topology_element.h"
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

class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzMaterial : public MaterialControl
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

    /// Pointer definition of HelmholtzMaterial
    KRATOS_CLASS_POINTER_DEFINITION(HelmholtzMaterial);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HelmholtzMaterial( std::string ControlName, Model& rModel, std::vector<LinearSolverType::Pointer>& rLinearSolvers, Parameters ControlSettings )
        :  MaterialControl(ControlName,rModel,ControlSettings){
            for(int lin_i=0;lin_i<rLinearSolvers.size();lin_i++)
                rLinearSystemSolvers.push_back(rLinearSolvers[lin_i]);
            mTechniqueSettings = ControlSettings["technique_settings"];
        }

    /// Destructor.
    virtual ~HelmholtzMaterial()
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
        KRATOS_INFO("HelmholtzMaterial:Initialize ") << "Starting initialization of material control "<<mControlName<<" ..." << std::endl;

        CreateModelParts();

        CalculateNodeNeighbourCount();
        
        for(int model_i=0;model_i<mpVMModelParts.size();model_i++){
            StrategyType* mpStrategy = new StrategyType (*mpVMModelParts[model_i],rLinearSystemSolvers[model_i]);            
            mpStrategy->Initialize();
            mpStrategies.push_back(mpStrategy);
        }

        physical_densities =  mTechniqueSettings["physical_densities"].GetVector();
        initial_density = mTechniqueSettings["initial_density"].GetDouble();
        youngs_modules =  mTechniqueSettings["youngs_modules"].GetVector();
        beta = mTechniqueSettings["beta"].GetDouble();        
        filtered_densities.resize(physical_densities.size());
        
        for(int i=0;i<physical_densities.size();i++){
            filtered_densities[i] = i;
        }
        double initial_filtered_density = ProjectBackward(initial_density,filtered_densities,physical_densities,beta);
        double initial_control_density = initial_filtered_density;
        double initial_young_modulus = ProjectForward(initial_filtered_density,filtered_densities,youngs_modules,beta);

        std::cout<<"initial_filtered_density : "<<initial_filtered_density<<std::endl;
        std::cout<<"initial_young_modulus : "<<initial_young_modulus<<std::endl;


        for(int model_i=0;model_i<mpVMModelParts.size();model_i++){
            SetVariable(mpVMModelParts[model_i],CD,initial_filtered_density); 
            SetVariable(mpVMModelParts[model_i],FD,initial_filtered_density); 
            SetVariable(mpVMModelParts[model_i],PD,initial_density); 
            SetVariable(mpVMModelParts[model_i],YOUNG_MODULUS,initial_young_modulus);
        }             

        KRATOS_INFO("HelmholtzMaterial:Initialize") << "Finished initialization of material control "<<mControlName<<" in " << timer.ElapsedSeconds() << " s." << std::endl;

    };
    // --------------------------------------------------------------------------
    void Update() override {

        opt_itr++;

        if (opt_itr % 20 == 0 && beta <20.0)
            beta *=1.5;
        if(beta>20.0)
            beta = 20.0;
        
        std::cout<<"++++++++++++++++++++++ beta : "<<beta<<" ++++++++++++++++++++++"<<std::endl;

        ComputeFilteredDensity();
        ComputePhyiscalDensity();
        ComputeYoungModulus();
        

    };  
    // --------------------------------------------------------------------------
    void MapControlUpdate(const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override{};
    // --------------------------------------------------------------------------
    void MapFirstDerivative(const Variable<double> &rDerivativeVariable, const Variable<double> &rMappedDerivativeVariable) override{


        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("HelmholtzMaterial:MapFirstDerivative") << "Starting mapping of " << rDerivativeVariable.Name() << "..." << std::endl;

        for(int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];
            //do the inverse projection and penalization
            for(auto& node_i : mpVMModelParts[model_i]->Nodes()){
                const auto& filtered_density = node_i.FastGetSolutionStepValue(FD);
                const auto& derivative = node_i.FastGetSolutionStepValue(rDerivativeVariable);
                auto& helmholtz_source = node_i.FastGetSolutionStepValue(HELMHOLTZ_SOURCE_DENSITY);
                if(rDerivativeVariable.Name()=="D_MASS_D_PD")
                    helmholtz_source = ProjectFirstDerivative(derivative,filtered_density,filtered_densities,physical_densities,beta);
                else
                    helmholtz_source = ProjectFirstDerivative(derivative,filtered_density,filtered_densities,youngs_modules,beta);
            }

            SetVariable(mpVMModePart,HELMHOLTZ_VAR_DENSITY,0.0);

            //now solve 
            mpStrategies[model_i]->Solve();
            SetVariable1ToVarible2(mpVMModePart,HELMHOLTZ_VAR_DENSITY,rMappedDerivativeVariable);
        }
        KRATOS_INFO("HelmholtzMaterial:MapFirstDerivative") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
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
        return "HelmholtzMaterial";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "HelmholtzMaterial";
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
    int opt_itr = 0;
    Vector physical_densities;
    Vector filtered_densities;
    double initial_density;
    Vector youngs_modules;    
    
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
        // creating vm model nodes and variables
        for(auto& control_obj : mControlSettings["controlling_objects"]){
            ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
            ModelPart& root_model_part = r_controlling_object.GetRootModelPart();
            std::string vm_model_part_name =  root_model_part.Name()+"_HELMHOLTZ_MATERIAL_Part";
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

            p_vm_model_part_property->SetValue(HELMHOLTZ_RADIUS_DENSITY,mTechniqueSettings["filter_radius"].GetDouble());

            for(auto& node : r_controlling_object.Nodes())
                p_vm_model_part->AddNode(&node);

            // creating elements
            ModelPart::ElementsContainerType &rmesh_elements = p_vm_model_part->Elements();   

            //check if the controlling model part has elements which have desnity value
            if(!r_controlling_object.Elements().size()>0)
                KRATOS_ERROR << "HelmholtzMaterial:CreateModelParts : controlling model part " <<control_obj.GetString()<<" does not have elements"<<std::endl;

            for (int i = 0; i < (int)r_controlling_object.Elements().size(); i++) {
                ModelPart::ElementsContainerType::iterator it = r_controlling_object.ElementsBegin() + i;
                const Properties& elem_i_prop = it->GetProperties();
                if(!elem_i_prop.Has(DENSITY))
                    KRATOS_ERROR << "HelmholtzMaterial:CreateModelParts : element " << it->Id()<<" of controlling model part "<<control_obj.GetString()<<" does not have desnity property !"<<std::endl;
                Element::Pointer p_element = new HelmholtzBulkTopologyElement(it->Id(), it->pGetGeometry(), p_vm_model_part_property);
                rmesh_elements.push_back(p_element);
            }   
        }

        // now add dofs
        for(int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];
            for(auto& node_i : mpVMModePart->Nodes())
            {
                node_i.AddDof(HELMHOLTZ_VAR_DENSITY);
            }
        }

    }

    void ComputeFilteredDensity(){   

        for(int model_i=0;model_i<mpVMModelParts.size();model_i++){

            //first update control density
            AddVariable1ToVarible2(mpVMModelParts[model_i],D_CD,CD);

            //now filter nodal control desnity 
            //first we need to multiply with mass matrix 
            SetVariable(mpVMModelParts[model_i],HELMHOLTZ_VAR_DENSITY,0.0);
            SetVariable(mpVMModelParts[model_i],HELMHOLTZ_SOURCE_DENSITY,0.0);
            //now we need to multiply with the mass matrix 
            for(auto& elem_i : mpVMModelParts[model_i]->Elements())
            {
                VectorType origin_values;
                GetElementVariableValuesVector(elem_i,CD,origin_values);
                MatrixType mass_matrix;
                elem_i.Calculate(HELMHOLTZ_MASS_MATRIX,mass_matrix,mpVMModelParts[model_i]->GetProcessInfo());           
                VectorType int_vals = prod(mass_matrix,origin_values);
                AddElementVariableValuesVector(elem_i,HELMHOLTZ_SOURCE_DENSITY,int_vals);
            }
            mpStrategies[model_i]->Solve();
            SetVariable1ToVarible2(mpVMModelParts[model_i],HELMHOLTZ_VAR_DENSITY,FD);
        }        
    } 

    void ComputePhyiscalDensity(){            
        for(int model_i=0;model_i<mpVMModelParts.size();model_i++){
            //now do the projection and then set the PD
            for(auto& node_i : mpVMModelParts[model_i]->Nodes()){
                const auto& filtered_density = node_i.FastGetSolutionStepValue(FD);
                auto& physical_density = node_i.FastGetSolutionStepValue(PD);
                physical_density = ProjectForward(filtered_density,filtered_densities,physical_densities,beta);
            }
        }
    }

    void ComputeYoungModulus(){      
        for(int model_i=0;model_i<mpVMModelParts.size();model_i++){
            //now do the projection and then set the PD
            for(auto& node_i : mpVMModelParts[model_i]->Nodes()){
                const auto& filtered_density = node_i.FastGetSolutionStepValue(FD);
                auto& youngs_modulus = node_i.FastGetSolutionStepValue(YOUNG_MODULUS);
                youngs_modulus = ProjectForward(filtered_density,filtered_densities,youngs_modules,beta);
            }
        }
    }

    double ProjectForward(double x,Vector x_limits,Vector y_limits,double beta){

        double x1,x2,y1,y2;
        int index_x1 = 0;
        if(x>=x_limits[x_limits.size()-1]){
            x1=x_limits[x_limits.size()-2];
            index_x1 = x_limits.size()-2;
            x2=x_limits[x_limits.size()-1];
            y1=y_limits[y_limits.size()-2];
            y2=y_limits[y_limits.size()-1];
        }
        else if(x<=x_limits[0]){
            x1=x_limits[0];
            index_x1 = 0;
            x2=x_limits[1];
            y1=y_limits[0];
            y2=y_limits[1];
        }
        else{
            for(int i=0;i<x_limits.size()-1;i++)
                if((x>=x_limits[i]) && (x<=x_limits[i+1]))
                {
                    y1 = y_limits[i];
                    y2 = y_limits[i+1];
                    x1 = x_limits[i];
                    index_x1 = i;
                    x2 = x_limits[i+1];
                    break;
                }            
        }        
        
        double pow_val = -2.0*beta*(x-(x1+x2)/2);

        if(index_x1>0){
            double prev_x1,prev_x2,prev_y1,prev_y2;
            prev_x1 = x_limits[index_x1-1];
            prev_x2 = x_limits[index_x1];
            prev_y1 = y_limits[index_x1-1];
            prev_y2 = y_limits[index_x1];
            double prev_pow_val = -2.0*beta*(x1-(prev_x1+prev_x2)/2);
            y1 = (prev_y2-prev_y1)/(1+std::exp(prev_pow_val)) + prev_y1;     
        }

        // if(pow_val<-600)
        //     pow_val = -600;
        // if(pow_val>600)
        //     pow_val = 600;

        return (y2-y1)/(1+std::exp(pow_val)) + y1;
    }


    double ProjectBackward(double y,Vector x_limits,Vector y_limits,double beta){
        
        double x = 0;
        if(y>=y_limits[y_limits.size()-1])
            x = x_limits[y_limits.size()-1];
        else if(y<=y_limits[0])
            x = x_limits[0];
        else{
            for(int i=0;i<y_limits.size()-1;i++)
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

    double ProjectFirstDerivative(double dfdy,double x,Vector x_limits,Vector y_limits,double beta){

        double dfdx = 0;
        double x1,x2,y1,y2;
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
            for(int i=0;i<x_limits.size()-1;i++)
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
        // if(pow_val<-600)
        //     pow_val = -600;
        // if(pow_val>600)
        //     pow_val = 600;

        double dydx = (1.0/(1+std::exp(pow_val))) * (1.0/(1+std::exp(pow_val))) * 2.0 * beta * std::exp(pow_val);

        if (y2<y1)
            dydx *=-1;

        return dfdy * dydx;

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
            auto& r_nodal_variable = rgeom[i_node].FastGetSolutionStepValue(rVariable);
            r_nodal_variable += rWeight * rValues[i_node];
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
        for(auto& node_i : mpVMModePart->Nodes())
        {
            auto& r_nodal_variable = node_i.FastGetSolutionStepValue(rVariable);
            r_nodal_variable = value;                   
        }
    }

    void SetVariable1ToVarible2(ModelPart* mpVMModePart,const Variable<double> &rVariable1,const Variable<double> &rVariable2) 
    {
        for(auto& node_i : mpVMModePart->Nodes())
        {
            auto& r_nodal_variable1 = node_i.FastGetSolutionStepValue(rVariable1);
            auto& r_nodal_variable2 = node_i.FastGetSolutionStepValue(rVariable2);
            r_nodal_variable2 = r_nodal_variable1;                  
        }
    }    

    void AddVariable1ToVarible2(ModelPart* mpVMModePart,const Variable<double> &rVariable1,const Variable<double> &rVariable2) 
    {
        for(auto& node_i : mpVMModePart->Nodes())
        {
            auto& r_nodal_variable1 = node_i.FastGetSolutionStepValue(rVariable1);
            auto& r_nodal_variable2 = node_i.FastGetSolutionStepValue(rVariable2);
            r_nodal_variable2 += r_nodal_variable1;                  
        }
    }    

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//      HelmholtzMaterial& operator=(HelmholtzMaterial const& rOther);

    /// Copy constructor.
//      HelmholtzMaterial(HelmholtzMaterial const& rOther);


    ///@}

}; // Class HelmholtzMaterial

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // HELMHOLTZ_MATERIAL_H
