// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef MASS_OPT_RESPONSE_H
#define MASS_OPT_RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "containers/model.h"
#include "includes/model_part.h"
#include "custom_responses/response.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "utilities/variable_utils.h"
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "spatial_containers/spatial_containers.h"
#include "processes/find_conditions_neighbours_process.h"

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

class KRATOS_API(OPTIMIZATION_APPLICATION) MassOptResponse : public Response
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

    /// Pointer definition of MassOptResponse
    KRATOS_CLASS_POINTER_DEFINITION(MassOptResponse);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MassOptResponse(std::string ResponseName, Model& rModel, Parameters& ResponseSettings )
        : Response(ResponseName,"mass",rModel, ResponseSettings){
            for(int i=0;i<mrResponseSettings["control_types"].size();i++){
                auto control_type = mrResponseSettings["control_types"][i].GetString();
                if(control_type=="shape"){
                    std::string gradient_mode = mrResponseSettings["gradient_settings"]["gradient_mode"].GetString();
                    if (gradient_mode.compare("finite_differencing") == 0)
                    {
                        double delta = mrResponseSettings["gradient_settings"]["step_size"].GetDouble();
                        mDelta = delta;
                    }
                    else
                        KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: finite_differencing" << std::endl;                    
                }
            }
        }

    /// Destructor.
    virtual ~MassOptResponse()
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
        for(int i=0;i<mrResponseSettings["evaluated_objects"].size();i++){
            auto eval_obj = mrResponseSettings["evaluated_objects"][i].GetString();
            ModelPart& eval_model_part = mrModel.GetModelPart(eval_obj);
            auto controlled_obj = mrResponseSettings["controlled_objects"][i].GetString();
            ModelPart& controlled_model_part = mrModel.GetModelPart(controlled_obj);
            auto control_type = mrResponseSettings["control_types"][i].GetString();

            KRATOS_ERROR_IF_NOT(eval_model_part.Elements().size()>0)
            <<"MassOptResponse::Initialize: evaluated object "<<eval_obj<<" must have elements !"<<std::endl;

            KRATOS_ERROR_IF_NOT(controlled_model_part.Elements().size()>0)
                <<"MassOptResponse::Initialize: controlled object "<<controlled_obj<<" for "<<control_type<<" sensitivity must have elements !"<<std::endl;
                   
        }
    };
    // --------------------------------------------------------------------------
    double CalculateValue() override {
        double total_mass = 0.0;
        for(auto& eval_obj : mrResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            const std::size_t domain_size = r_eval_object.GetProcessInfo()[DOMAIN_SIZE];
			for (auto& elem_i : r_eval_object.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active)
                    total_mass += CalculateElementMass(elem_i);
			}
        }
        return total_mass;
    };    

    double CalculateElementMass(Element& elem_i){
        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();

        // We copy the current coordinates and move the coordinates to the initial configuration
        std::vector<array_1d<double, 3>> current_coordinates(number_of_nodes);
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(current_coordinates[i_node]) = r_this_geometry[i_node].Coordinates();
            noalias(r_this_geometry[i_node].Coordinates()) = r_this_geometry[i_node].GetInitialPosition().Coordinates();
        }

        double element_mass = 0.0;
        VectorType element_nodal_thickness_density;
        GetElementsNodalThicknessDensity(elem_i,element_nodal_thickness_density);

        const auto& integrationPoints = r_this_geometry.IntegrationPoints(r_this_geometry.GetDefaultIntegrationMethod());
        const unsigned int numberOfIntegrationPoints = integrationPoints.size();
        const Matrix& N_container = r_this_geometry.ShapeFunctionsValues(r_this_geometry.GetDefaultIntegrationMethod());                    

        for ( unsigned int pointNumber = 0; pointNumber < numberOfIntegrationPoints; pointNumber++ )
        {
            // Get FEM-shape-function-value for current integration point
            Vector N_FEM_GPi = row( N_container, pointNumber);
            double integration_weight = integrationPoints[pointNumber].Weight() * r_this_geometry.DeterminantOfJacobian(pointNumber,r_this_geometry.GetDefaultIntegrationMethod());
            element_mass += integration_weight * inner_prod(N_FEM_GPi,element_nodal_thickness_density);
        }

        // We restore the current configuration
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(r_this_geometry[i_node].Coordinates()) = current_coordinates[i_node];
        }

        return element_mass;
    }

    // --------------------------------------------------------------------------
    void CalculateGradient() override {

		KRATOS_TRY;

        for(int i=0;i<mrResponseSettings["controlled_objects"].size();i++){
            auto controlled_obj = mrResponseSettings["controlled_objects"][i].GetString();
            ModelPart& controlled_model_part = mrModel.GetModelPart(controlled_obj);
            auto control_type = mrResponseSettings["control_types"][i].GetString();
            if(control_type=="shape")
                VariableUtils().SetHistoricalVariableToZero(D_MASS_D_X, controlled_model_part.Nodes());
            else if(control_type=="material")
                VariableUtils().SetHistoricalVariableToZero(D_MASS_D_FD, controlled_model_part.Nodes());
            else if(control_type=="thickness")
                VariableUtils().SetHistoricalVariableToZero(D_MASS_D_PT, controlled_model_part.Nodes());

			for (auto& elem_i : controlled_model_part.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active){
                    if(control_type=="shape")
                        CalculateElementShapeGradients(elem_i);
                    if(control_type=="material")
                        CalculateElementMaterialGradients(elem_i);
                    if(control_type=="thickness")
                        CalculateElementThicknessGradients(elem_i);                                                
                }
            }

            
        }

		KRATOS_CATCH("");
 
    };  

    void CalculateElementShapeGradients(Element& elem_i){

        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();

        // We copy the current coordinates and move the coordinates to the initial configuration
        std::vector<array_1d<double, 3>> current_coordinates(number_of_nodes);
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(current_coordinates[i_node]) = r_this_geometry[i_node].Coordinates();
            noalias(r_this_geometry[i_node].Coordinates()) = r_this_geometry[i_node].GetInitialPosition().Coordinates();
        }
        double mass_before_fd = CalculateElementMass(elem_i);

        for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
            auto& node_grads = r_this_geometry[i_node].FastGetSolutionStepValue(D_MASS_D_X);	

			r_this_geometry[i_node].X() += mDelta;
			r_this_geometry[i_node].X0() += mDelta;
            double mass_after_fd = CalculateElementMass(elem_i);
			node_grads[0] += (mass_after_fd - mass_before_fd) / mDelta;
			r_this_geometry[i_node].X() -= mDelta;
			r_this_geometry[i_node].X0() -= mDelta;

			r_this_geometry[i_node].Y() += mDelta;
			r_this_geometry[i_node].Y0() += mDelta;
            mass_after_fd = CalculateElementMass(elem_i);
			node_grads[1] += (mass_after_fd - mass_before_fd) / mDelta;
			r_this_geometry[i_node].Y() -= mDelta;
			r_this_geometry[i_node].Y0() -= mDelta;   

			r_this_geometry[i_node].Z() += mDelta;
			r_this_geometry[i_node].Z0() += mDelta;
            mass_after_fd = CalculateElementMass(elem_i);
			node_grads[2] += (mass_after_fd - mass_before_fd) / mDelta;
			r_this_geometry[i_node].Z() -= mDelta;
			r_this_geometry[i_node].Z0() -= mDelta;                                  
        }

        // We restore the current configuration
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(r_this_geometry[i_node].Coordinates()) = current_coordinates[i_node];
        }


    };

    void CalculateElementMaterialGradients(Element& elem_i){

        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();

        // We copy the current coordinates and move the coordinates to the initial configuration
        std::vector<array_1d<double, 3>> current_coordinates(number_of_nodes);
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(current_coordinates[i_node]) = r_this_geometry[i_node].Coordinates();
            noalias(r_this_geometry[i_node].Coordinates()) = r_this_geometry[i_node].GetInitialPosition().Coordinates();
        }

        VectorType element_nodal_thickness;
        GetElementsNodalThickness(elem_i,element_nodal_thickness);

        const auto& integrationPoints = r_this_geometry.IntegrationPoints(r_this_geometry.GetDefaultIntegrationMethod());
        const unsigned int numberOfIntegrationPoints = integrationPoints.size();
        const Matrix& N_container = r_this_geometry.ShapeFunctionsValues(r_this_geometry.GetDefaultIntegrationMethod());                    

        for ( unsigned int pointNumber = 0; pointNumber < numberOfIntegrationPoints; pointNumber++ )
        {
            // Get FEM-shape-function-value for current integration point
            Vector N_FEM_GPi = row( N_container, pointNumber);
            double integration_weight = integrationPoints[pointNumber].Weight() * r_this_geometry.DeterminantOfJacobian(pointNumber,r_this_geometry.GetDefaultIntegrationMethod());
            for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
                auto& node_val = r_this_geometry[i_node].FastGetSolutionStepValue(D_MASS_D_FD);
                auto& d_pd_d_fd = r_this_geometry[i_node].FastGetSolutionStepValue(D_PD_D_FD);
                node_val += integration_weight * N_FEM_GPi[i_node] * element_nodal_thickness[i_node] * d_pd_d_fd;
            }
        }

        // We restore the current configuration
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(r_this_geometry[i_node].Coordinates()) = current_coordinates[i_node];
        }
    };        

    void CalculateElementThicknessGradients(Element& elem_i){
        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();

        // We copy the current coordinates and move the coordinates to the initial configuration
        std::vector<array_1d<double, 3>> current_coordinates(number_of_nodes);
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(current_coordinates[i_node]) = r_this_geometry[i_node].Coordinates();
            noalias(r_this_geometry[i_node].Coordinates()) = r_this_geometry[i_node].GetInitialPosition().Coordinates();
        }

        VectorType element_nodal_density;
        GetElementsNodalDensity(elem_i,element_nodal_density);

        const auto& integrationPoints = r_this_geometry.IntegrationPoints(r_this_geometry.GetDefaultIntegrationMethod());
        const unsigned int numberOfIntegrationPoints = integrationPoints.size();
        const Matrix& N_container = r_this_geometry.ShapeFunctionsValues(r_this_geometry.GetDefaultIntegrationMethod());                    

        for ( unsigned int pointNumber = 0; pointNumber < numberOfIntegrationPoints; pointNumber++ )
        {
            // Get FEM-shape-function-value for current integration point
            Vector N_FEM_GPi = row( N_container, pointNumber);
            double integration_weight = integrationPoints[pointNumber].Weight() * r_this_geometry.DeterminantOfJacobian(pointNumber,r_this_geometry.GetDefaultIntegrationMethod());
            for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
                auto& node_val = r_this_geometry[i_node].FastGetSolutionStepValue(D_MASS_D_PT);
                node_val += integration_weight * N_FEM_GPi[i_node] * element_nodal_density[i_node];
            }
        }

        // We restore the current configuration
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(r_this_geometry[i_node].Coordinates()) = current_coordinates[i_node];
        }        
    };

    void GetElementsNodalThicknessDensity(const Element& rElement, VectorType &rValues){
        const GeometryType &rgeom = rElement.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const Properties& this_properties = rElement.GetProperties();
        const unsigned int local_size = num_nodes;    
        if (rValues.size() != local_size)
            rValues.resize(local_size, false);   

        for (SizeType i_node = 0; i_node < num_nodes; ++i_node){
            if(rgeom[i_node].SolutionStepsDataHas(PD)){
                rValues[i_node] = rgeom[i_node].FastGetSolutionStepValue(PD);;
            }                  
            else if(this_properties.Has(DENSITY))
                rValues[i_node] = this_properties[DENSITY];
            else
                rValues[i_node] = 1;
        }

        for (SizeType i_node = 0; i_node < num_nodes; ++i_node){
            if(rgeom[i_node].SolutionStepsDataHas(PT))
                rValues[i_node] *= rgeom[i_node].FastGetSolutionStepValue(PT);  
            else if(this_properties.Has(THICKNESS))
                rValues[i_node] *= this_properties[THICKNESS];
            else
                rValues[i_node] *= 1;
        }

    }; 

    void GetElementsNodalThickness(const Element& rElement, VectorType &rValues){
        const GeometryType &rgeom = rElement.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const Properties& this_properties = rElement.GetProperties();
        const unsigned int local_size = num_nodes;    
        if (rValues.size() != local_size)
            rValues.resize(local_size, false);   

        for (SizeType i_node = 0; i_node < num_nodes; ++i_node){
            if(rgeom[i_node].SolutionStepsDataHas(PT))
                rValues[i_node] = rgeom[i_node].FastGetSolutionStepValue(PT);  
            else if(this_properties.Has(THICKNESS))
                rValues[i_node] = this_properties[THICKNESS];
            else
                rValues[i_node] = 1;
        }
    }; 

    void GetElementsNodalDensity(const Element& rElement, VectorType &rValues){
        const GeometryType &rgeom = rElement.GetGeometry();
        const SizeType num_nodes = rgeom.PointsNumber();
        const Properties& this_properties = rElement.GetProperties();
        const unsigned int local_size = num_nodes;    
        if (rValues.size() != local_size)
            rValues.resize(local_size, false);   

        for (SizeType i_node = 0; i_node < num_nodes; ++i_node){
            if(rgeom[i_node].SolutionStepsDataHas(PD))
                rValues[i_node] = rgeom[i_node].FastGetSolutionStepValue(PD);  
            else if(this_properties.Has(DENSITY))
                rValues[i_node] = this_properties[DENSITY];
            else
                rValues[i_node] = 1;
        }
    };

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
        return "MassOptResponse";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MassOptResponse";
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
    double mDelta;

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
//      MassOptResponse& operator=(MassOptResponse const& rOther);

    /// Copy constructor.
//      MassOptResponse(MassOptResponse const& rOther);


    ///@}

}; // Class MassOptResponse

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MASS_OPT_RESPONSE_H
