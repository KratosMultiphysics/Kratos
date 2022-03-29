// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef HELMHOLTZ_TOPOLOGY_H
#define HELMHOLTZ_TOPOLOGY_H

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
#include "custom_controls/topology_controls/topology_control.h"
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

class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzTopology : public TopologyControl
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

    /// Pointer definition of HelmholtzTopology
    KRATOS_CLASS_POINTER_DEFINITION(HelmholtzTopology);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HelmholtzTopology( std::string ControlName, Model& rModel, std::vector<LinearSolverType::Pointer>& rLinearSolvers, Parameters ControlSettings )
        :  TopologyControl(ControlName,rModel,ControlSettings){
            for(int lin_i=0;lin_i<rLinearSolvers.size();lin_i++)
                rLinearSystemSolvers.push_back(rLinearSolvers[lin_i]);
            mTechniqueSettings = ControlSettings["technique_settings"];
        }

    /// Destructor.
    virtual ~HelmholtzTopology()
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
        KRATOS_INFO("HelmholtzTopology:Initialize ") << "Starting initialization of topology control "<<mControlName<<" ..." << std::endl;

        CreateModelParts();

        CalculateNodeNeighbourCount();
        
        for(int model_i=0;model_i<mpVMModelParts.size();model_i++){
            StrategyType* mpStrategy = new StrategyType (*mpVMModelParts[model_i],rLinearSystemSolvers[model_i]);            
            mpStrategy->Initialize();
            mpStrategies.push_back(mpStrategy);
        }

        AdjustFilterSizes();

        //initialize control topology
        for(int model_i=0;model_i<mpVMModelParts.size();model_i++)
            SetVariable(mpVMModelParts[model_i],CD,0.5);

        //now compute physical density
        beta = mTechniqueSettings["beta"].GetDouble();
        sigmoid_projection = mTechniqueSettings["sigmoid_projection"].GetBool();
        penalization = mTechniqueSettings["penalization"].GetBool();
        opt_itr = 0;
        ComputePhyiscalDensity();

        KRATOS_INFO("HelmholtzTopology:Initialize") << "Finished initialization of desnity control "<<mControlName<<" in " << timer.ElapsedSeconds() << " s." << std::endl;

    };
    // --------------------------------------------------------------------------
    void Update() override {
        ComputePhyiscalDensity();
    };  
    // --------------------------------------------------------------------------
    void MapControlUpdate(const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override{};
    // --------------------------------------------------------------------------
    void MapFirstDerivative(const Variable<double> &rDerivativeVariable, const Variable<double> &rMappedDerivativeVariable) override{


        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("HelmholtzTopology:MapFirstDerivative") << "Starting mapping of " << rDerivativeVariable.Name() << "..." << std::endl;

        for(int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];
            //do the inverse projection and penalization
            for(auto& node_i : mpVMModelParts[model_i]->Nodes()){
                const auto& filtered_density = node_i.FastGetSolutionStepValue(FD);
                const auto& derivative = node_i.FastGetSolutionStepValue(rDerivativeVariable);
                auto& helmholtz_source = node_i.FastGetSolutionStepValue(HELMHOLTZ_SOURCE_DENSITY);
                if(sigmoid_projection)
                {
                    double value = -2*beta*(filtered_density-0.5);
                    if(value<-600)
                        value = -600;
                    if(value>600)
                        value = 600;
                    helmholtz_source = derivative * (2.0*beta*std::exp(value)) * std::pow(1.0/(1+std::exp(value)),2);

                }
                else if(penalization)
                {
                    helmholtz_source = 3 * filtered_density * filtered_density * derivative;
                }
                else
                    helmholtz_source = derivative;
            }

            SetVariable(mpVMModePart,HELMHOLTZ_VAR_DENSITY,0.0);

            //now solve 
            mpStrategies[model_i]->Solve();

            SetVariable1ToVarible2(mpVMModePart,HELMHOLTZ_VAR_DENSITY,rMappedDerivativeVariable);
        }

    
        KRATOS_INFO("HelmholtzTopology:MapFirstDerivative") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;

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
        return "HelmholtzTopology";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "HelmholtzTopology";
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
    int opt_itr;
    bool sigmoid_projection;
    bool penalization;
    
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
            std::string vm_model_part_name =  root_model_part.Name()+"_HELMHOLTZ_TOPOLOGY_Part";
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
                KRATOS_ERROR << "HelmholtzTopology:CreateModelParts : controlling model part " <<control_obj.GetString()<<" does not have elements"<<std::endl;

            for (int i = 0; i < (int)r_controlling_object.Elements().size(); i++) {
                ModelPart::ElementsContainerType::iterator it = r_controlling_object.ElementsBegin() + i;
                const Properties& elem_i_prop = it->GetProperties();
                if(!elem_i_prop.Has(DENSITY))
                    KRATOS_ERROR << "HelmholtzTopology:CreateModelParts : element " << it->Id()<<" of controlling model part "<<control_obj.GetString()<<" does not have desnity property !"<<std::endl;
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

    void AdjustFilterSizes()
    {
        // now adjust the filter size if its adaptive
        if(mTechniqueSettings["automatic_filter_size"].GetBool()){
            for(auto& control_obj : mControlSettings["controlling_objects"]){
                ModelPart& r_controlling_object = mrModel.GetModelPart(control_obj.GetString());
                ModelPart& root_model_part = r_controlling_object.GetRootModelPart();
                std::string vm_model_part_name =  root_model_part.Name()+"_HELMHOLTZ_TOPOLOGY_Part";
                ModelPart* mpVMModePart = &(root_model_part.GetSubModelPart(vm_model_part_name));

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
                for(auto& elem_i : mpVMModePart->Elements())
                    elem_i.GetProperties().SetValue(HELMHOLTZ_RADIUS_DENSITY,surface_filter_size);

                KRATOS_INFO("HelmholtzTopology:AdjustFilterSizes") << " Surface filter of "<<control_obj.GetString() <<" is adjusted to " << surface_filter_size << std::endl;
            }
        }    
    }

    void ComputePhyiscalDensity(){            

        //now initialize control and physical density
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

            //now do the projection and then set the PD
            for(auto& node_i : mpVMModelParts[model_i]->Nodes()){
                const auto& filtered_density = node_i.FastGetSolutionStepValue(FD);
                auto& physical_density = node_i.FastGetSolutionStepValue(PD);
                auto& density = node_i.FastGetSolutionStepValue(DENSITY);
                if (sigmoid_projection)
                {
                    double value = -2*beta*(filtered_density-0.5); 
                    if(value<-600)
                        value = -600;
                    if(value>600)
                        value = 600;                    
                    physical_density = (1.0/(1+std::exp(value)));
                    if(physical_density<0.001)
                        physical_density = 0.001;

                }
                else if(penalization){
                    physical_density = std::pow(filtered_density,3);
                }
                else
                    physical_density = filtered_density;

                density = physical_density;

            }
        }
        opt_itr++;
        
        // if(opt_itr==20)
        //     beta *=1.5;
        // if(opt_itr==40)
        //     beta *=1.5;
        // if(opt_itr==80)
        //     beta *=1.5;
        // if(opt_itr==160)
        //     beta *=1.5; 
        // if(opt_itr==300)
        //     beta *=1.5; 
        // if(opt_itr==600)
        //     beta *=1.5;           

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
//      HelmholtzTopology& operator=(HelmholtzTopology const& rOther);

    /// Copy constructor.
//      HelmholtzTopology(HelmholtzTopology const& rOther);


    ///@}

}; // Class HelmholtzTopology

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // HELMHOLTZ_TOPOLOGY_H
