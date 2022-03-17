// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

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
            for(int lin_i=0;lin_i<rLinearSolvers.size();lin_i++)
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

        CreateModelParts();

        CalculateNodeNeighbourCount();
        
        for(int model_i=0;model_i<mpVMModelParts.size();model_i++){
            StrategyType* mpStrategy = new StrategyType (*mpVMModelParts[model_i],rLinearSystemSolvers[model_i]);            
            mpStrategy->Initialize();
            mpStrategies.push_back(mpStrategy);
        }

        AdjustFilterSizes();

        //initialize control thicknesses with 
        t_min = mTechniqueSettings["min_thickness"].GetInt();
        t_max = mTechniqueSettings["max_thickness"].GetInt();
        t_initial = mTechniqueSettings["initial_thickness"].GetDouble();
        for(int model_i=0;model_i<mpVMModelParts.size();model_i++)
            SetVariable(mpVMModelParts[model_i],CT,t_initial);

        //now compute physical thicknesses
        beta = mTechniqueSettings["beta"].GetDouble();
        sigmoid_projection = mTechniqueSettings["sigmoid_projection"].GetBool();
        penalization = mTechniqueSettings["penalization"].GetBool();
        ComputePhyiscalThickness();

        KRATOS_INFO("HelmholtzThickness:Initialize") << "Finished initialization of thickness control "<<mControlName<<" in " << timer.ElapsedSeconds() << " s." << std::endl;

    };
    // --------------------------------------------------------------------------
    void Update() override {
        ComputePhyiscalThickness();
    };  
    // --------------------------------------------------------------------------
    void MapControlUpdate(const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override{};
    // --------------------------------------------------------------------------
    void MapFirstDerivative(const Variable<double> &rDerivativeVariable, const Variable<double> &rMappedDerivativeVariable) override{


        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("HelmholtzThickness:MapFirstDerivative") << "Starting mapping of " << rDerivativeVariable.Name() << "..." << std::endl;

        for(int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];
            //do the inverse projection and penalization
            for(auto& node_i : mpVMModelParts[model_i]->Nodes()){
                const auto& filtered_thickness = node_i.FastGetSolutionStepValue(FT);
                const auto& derivative = node_i.FastGetSolutionStepValue(rDerivativeVariable);
                auto& helmholtz_source = node_i.FastGetSolutionStepValue(HELMHOLTZ_SOURCE_THICKNESS);
                if(sigmoid_projection)
                {

                    double reference_thickness = 0;
                    if (filtered_thickness<=t_min)
                        reference_thickness = t_min;
                    else if (filtered_thickness>=t_max)
                        reference_thickness = t_max-1;
                    else
                    {
                        for(int t=t_min;t<t_max;t++){                        
                            if (filtered_thickness>=t && filtered_thickness<=t+1)
                            {
                                reference_thickness = t;
                                break;
                            }
                        }
                    }

                    double value = -2*beta*(filtered_thickness-(reference_thickness+0.5));
                    if(value<-600)
                        value = -600;
                    if(value>600)
                        value = 600;
                    helmholtz_source = derivative * (2.0*beta*std::exp(value)) * std::pow(1.0/(1+std::exp(value)),2);

                }
                else if(penalization)
                {
                    helmholtz_source = 3 * filtered_thickness * filtered_thickness * derivative;
                }
                else
                    helmholtz_source = derivative;

            }

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
    Parameters mTechniqueSettings;
    double beta;
    bool sigmoid_projection;
    bool penalization;
    double t_min;
    double t_max;
    double t_initial;
    
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
            std::string vm_model_part_name =  root_model_part.Name()+"_Helmholtz_Thickness_Part";
            ModelPart* p_vm_model_part;

            if (root_model_part.HasSubModelPart(vm_model_part_name))
                p_vm_model_part = &(root_model_part.GetSubModelPart(vm_model_part_name));
            else{
                p_vm_model_part = &(root_model_part.CreateSubModelPart(vm_model_part_name));
                mpVMModelParts.push_back(p_vm_model_part);
            }

            for(auto& node : r_controlling_object.Nodes())
                p_vm_model_part->AddNode(&node);

            // creating elements
            ModelPart::ElementsContainerType &rmesh_elements = p_vm_model_part->Elements();   

            //check if the controlling model part has elements which have thickness value
            if(!r_controlling_object.Elements().size()>0)
                KRATOS_ERROR << "HelmholtzThickness:CreateModelParts : controlling model part " <<control_obj.GetString()<<" does not have elements"<<std::endl;

            for (int i = 0; i < (int)r_controlling_object.Elements().size(); i++) {
                ModelPart::ElementsContainerType::iterator it = r_controlling_object.ElementsBegin() + i;
                const Properties& elem_i_prop = it->GetProperties();
                if(!elem_i_prop.Has(THICKNESS))
                    KRATOS_ERROR << "HelmholtzThickness:CreateModelParts : element " << it->Id()<<" of controlling model part "<<control_obj.GetString()<<" does not have thickness property !"<<std::endl;
                Properties::Pointer model_part_new_prop = r_controlling_object.CreateNewProperties(r_controlling_object.NumberOfProperties()+1);
                *model_part_new_prop = elem_i_prop;
                model_part_new_prop->SetValue(HELMHOLTZ_RADIUS_THICKNESS,mTechniqueSettings["filter_radius"].GetDouble());
                it->SetProperties(model_part_new_prop);
                Element::Pointer p_element = new HelmholtzSurfThicknessElement(it->Id(), it->pGetGeometry(), model_part_new_prop);
                rmesh_elements.push_back(p_element);
            }  
        }

        // now add dofs
        for(int model_i =0;model_i<mpVMModelParts.size();model_i++)
        {
            ModelPart* mpVMModePart = mpVMModelParts[model_i];
            for(auto& node_i : mpVMModePart->Nodes())
            {
                node_i.AddDof(HELMHOLTZ_VAR_THICKNESS);
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
                std::string vm_model_part_name =  root_model_part.Name()+"_Helmholtz_Thickness_Part";
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
                    elem_i.GetProperties().SetValue(HELMHOLTZ_RADIUS_THICKNESS,surface_filter_size);

                KRATOS_INFO("HelmholtzThickness:AdjustFilterSizes") << " Surface filter of "<<control_obj.GetString() <<" is adjusted to " << surface_filter_size << std::endl;
            }
        }    
    }

    void ComputePhyiscalThickness(){

        //now initialize control and physical thicknesses
        for(int model_i=0;model_i<mpVMModelParts.size();model_i++){
            //first update control thicknesses
            AddVariable1ToVarible2(mpVMModelParts[model_i],D_CT,CT);

            //now filter nodal control thickness 
            //first we need to multiply with mass matrix 
            SetVariable(mpVMModelParts[model_i],HELMHOLTZ_VAR_THICKNESS,0.0);
            SetVariable(mpVMModelParts[model_i],HELMHOLTZ_SOURCE_THICKNESS,0.0);
            //now we need to multiply with the mass matrix 
            for(auto& elem_i : mpVMModelParts[model_i]->Elements())
            {
                VectorType origin_values;
                GetElementVariableValuesVector(elem_i,CT,origin_values);
                MatrixType mass_matrix;
                elem_i.Calculate(HELMHOLTZ_MASS_MATRIX,mass_matrix,mpVMModelParts[model_i]->GetProcessInfo());            
                VectorType int_vals = prod(mass_matrix,origin_values);
                AddElementVariableValuesVector(elem_i,HELMHOLTZ_SOURCE_THICKNESS,int_vals);
            }

            mpStrategies[model_i]->Solve();

            SetVariable1ToVarible2(mpVMModelParts[model_i],HELMHOLTZ_VAR_THICKNESS,FT);

            //now do the projection and then set the PT
            for(auto& node_i : mpVMModelParts[model_i]->Nodes()){
                const auto& filtered_thickness = node_i.FastGetSolutionStepValue(FT);
                auto& physical_thickness = node_i.FastGetSolutionStepValue(PT);
                if (sigmoid_projection)
                {
                    double reference_thickness = 0;
                    if (filtered_thickness<=t_min)
                        reference_thickness = t_min;
                    else if (filtered_thickness>=t_max)
                        reference_thickness = t_max-1;
                    else
                    {
                        for(int t=t_min;t<t_max;t++){                        
                            if (filtered_thickness>=t && filtered_thickness<=t+1)
                            {
                                reference_thickness = t;
                                break;
                            }
                        }
                    }

                    double value = -2*beta*(filtered_thickness-(reference_thickness+0.5)); 
                    if(value<-600)
                        value = -600;
                    if(value>600)
                        value = 600;                    
                    physical_thickness = (1.0/(1+std::exp(value)))+reference_thickness;
                    if(physical_thickness<0.001)
                        physical_thickness = 0.001;

                }
                else if(penalization){
                    physical_thickness = std::pow(filtered_thickness,3);
                }
                else
                    physical_thickness = filtered_thickness;
            }
            //now compute element thicknesses
            for(auto& elem_i : mpVMModelParts[model_i]->Elements()){

                const auto& r_geom = elem_i.GetGeometry();	
                const auto& integration_method = r_geom.GetDefaultIntegrationMethod();
                const auto& integration_points = r_geom.IntegrationPoints(integration_method);
                const unsigned int NumGauss = integration_points.size();
                Vector GaussPtsJDet = ZeroVector(NumGauss);
                r_geom.DeterminantOfJacobian(GaussPtsJDet, integration_method);
                const auto& Ncontainer = r_geom.ShapeFunctionsValues(integration_method); 

                double elem_area = 0.0; 
                for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
                    elem_area += integration_points[i_point].Weight() * GaussPtsJDet[i_point];

                double element_thickness = 0;
                for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
                {
                	const double IntToReferenceWeight = integration_points[i_point].Weight() * GaussPtsJDet[i_point];
                	const auto& rN = row(Ncontainer,i_point);
                	int node_index = 0;
                	for (auto& node_i : r_geom){
                        const auto& nodal_thickness = node_i.FastGetSolutionStepValue(PT);
                        element_thickness += nodal_thickness * rN[node_index] * IntToReferenceWeight / elem_area;
                		node_index++;
                	}						
                }
                elem_i.GetProperties().SetValue(THICKNESS,element_thickness);                                    
            }

        }

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
