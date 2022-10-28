// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef LINEAR_STRAIN_ENERGY_OPT_RESPONSE_H
#define LINEAR_STRAIN_ENERGY_OPT_RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------


// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "custom_responses/response.h"

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

class KRATOS_API(OPTIMIZATION_APPLICATION) LinearStrainEnergyOptResponse : public Response
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearStrainEnergyOptResponse
    KRATOS_CLASS_POINTER_DEFINITION(LinearStrainEnergyOptResponse);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LinearStrainEnergyOptResponse(std::string ResponseName, Model& rModel, Parameters& ResponseSettings )
        : Response(ResponseName,"linear_strain_energy",rModel, ResponseSettings){
            for(int i=0;i<mrResponseSettings["control_types"].size();i++){
                auto control_type = mrResponseSettings["control_types"][i].GetString();
                if(control_type=="shape"){
                    std::string gradient_mode = mrResponseSettings["gradient_settings"]["gradient_mode"].GetString();
                    if (gradient_mode.compare("semi_analytic") == 0)
                    {
                        double delta = mrResponseSettings["gradient_settings"]["step_size"].GetDouble();
                        mDelta = delta;
                    }
                    else
                        KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;                    
                }
            }
        }

    /// Destructor.
    virtual ~LinearStrainEnergyOptResponse()
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
            <<"LinearStrainEnergyOptResponse::Initialize: evaluated object "<<eval_obj<<" must have elements !"<<std::endl;

            KRATOS_ERROR_IF_NOT(controlled_model_part.Elements().size()>0)
                <<"LinearStrainEnergyOptResponse::Initialize: controlled object "<<controlled_obj<<" for "<<control_type<<" sensitivity must have elements !"<<std::endl;
                   
        }
    };
    // --------------------------------------------------------------------------
    double CalculateValue() override {
        double total_strain_energy = 0.0;
        for(auto& eval_obj : mrResponseSettings["evaluated_objects"]){
            ModelPart& r_eval_object = mrModel.GetModelPart(eval_obj.GetString());
            const ProcessInfo &CurrentProcessInfo = r_eval_object.GetProcessInfo();
            const std::size_t domain_size = r_eval_object.GetProcessInfo()[DOMAIN_SIZE];
            // Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
            for (auto& elem_i : r_eval_object.Elements())
            {
                const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
                if(element_is_active)
                {
                    Matrix LHS;
                    Vector RHS;
                    Vector u;

                    // Get state solution relevant for energy calculation
                    const auto& rConstElemRef = elem_i;
                    rConstElemRef.GetValuesVector(u,0);

                    elem_i.CalculateLocalSystem(LHS,RHS,CurrentProcessInfo);

                    // Compute strain energy
                    total_strain_energy += 0.5 * inner_prod(u,prod(LHS,u));
                }
            }
        }
        return total_strain_energy;
    };    


    // --------------------------------------------------------------------------
    void CalculateGradient() override {

		KRATOS_TRY;

        for(int i=0;i<mrResponseSettings["controlled_objects"].size();i++){
            auto controlled_obj = mrResponseSettings["controlled_objects"][i].GetString();
            ModelPart& controlled_model_part = mrModel.GetModelPart(controlled_obj);
            const ProcessInfo &CurrentProcessInfo = controlled_model_part.GetProcessInfo();
            auto control_type = mrResponseSettings["control_types"][i].GetString();
            std::string grad_field_name;
            if(control_type=="shape"){
                grad_field_name = mrResponseSettings["gradient_settings"]["shape_gradient_field_name"].GetString();
                VariableUtils().SetHistoricalVariableToZero(KratosComponents<Variable<array_1d<double,3>>>::Get(grad_field_name), controlled_model_part.Nodes());
            }                
            else if(control_type=="material"){
                grad_field_name = mrResponseSettings["gradient_settings"]["material_gradient_field_name"].GetString();
                VariableUtils().SetHistoricalVariableToZero(KratosComponents<Variable<double>>::Get(grad_field_name), controlled_model_part.Nodes());                
            }
            else if(control_type=="thickness"){
                grad_field_name = mrResponseSettings["gradient_settings"]["thickness_gradient_field_name"].GetString();
                VariableUtils().SetHistoricalVariableToZero(KratosComponents<Variable<double>>::Get(grad_field_name), controlled_model_part.Nodes());                 
            }

            //elems
			for (auto& elem_i : controlled_model_part.Elements()){
				const bool element_is_active = elem_i.IsDefined(ACTIVE) ? elem_i.Is(ACTIVE) : true;
				if(element_is_active){
                    if(control_type=="shape")
                        CalculateElementShapeGradients(elem_i,grad_field_name,CurrentProcessInfo);
                    if(control_type=="material")
                        CalculateElementMaterialGradients(elem_i,grad_field_name,CurrentProcessInfo);
                    if(control_type=="thickness")
                        CalculateElementThicknessGradients(elem_i,grad_field_name,CurrentProcessInfo);                                                
                }
            }

            //conds
			for (auto& cond_i : controlled_model_part.Conditions()){
				const bool cond_is_active = cond_i.IsDefined(ACTIVE) ? cond_i.Is(ACTIVE) : true;
				if(cond_is_active){
                    if(control_type=="shape")
                        CalculateConditionShapeGradients(cond_i,grad_field_name,CurrentProcessInfo);                                                               
                }
            }
            
        }

		KRATOS_CATCH("");
 
    };  

    void CalculateElementShapeGradients(Element& elem_i, std::string shape_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();

        Vector u;
        Vector lambda;
        Vector RHS;

        // Get state solution
        const auto& rConstElemRef = elem_i;
        rConstElemRef.GetValuesVector(u,0);

        // Get adjoint variables (Corresponds to 1/2*u)
        lambda = 0.5*u;

        // Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
        elem_i.CalculateRightHandSide(RHS, rCurrentProcessInfo);
        for (auto& node_i : elem_i.GetGeometry())
        {
            array_3d gradient_contribution(3, 0.0);
            Vector RHS_perturbed = Vector(RHS.size());
            Vector derived_RHS = Vector(RHS.size());
            
            // x-direction
            node_i.GetInitialPosition()[0] += mDelta;
            node_i.Coordinates()[0] += mDelta;
            elem_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[0] -= mDelta;
            node_i.Coordinates()[0] -= mDelta;
            gradient_contribution[0] = inner_prod(lambda, derived_RHS);


            // y-direction
            node_i.GetInitialPosition()[1] += mDelta;
            node_i.Coordinates()[1] += mDelta;
            elem_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[1] -= mDelta;
            node_i.Coordinates()[1] -= mDelta;
            gradient_contribution[1] = inner_prod(lambda, derived_RHS);

            // z-direction
            node_i.GetInitialPosition()[2] += mDelta;
            node_i.Coordinates()[2] += mDelta;
            elem_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[2] -= mDelta;
            node_i.Coordinates()[2] -= mDelta;
            gradient_contribution[2] = inner_prod(lambda, derived_RHS);

            // Assemble sensitivity to node
            noalias(node_i.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name))) += gradient_contribution;
        }

    };

    void CalculateConditionShapeGradients(Condition& cond_i, std::string shape_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        // We get the element geometry
        auto& r_this_geometry = cond_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();

        Vector u;
        Vector lambda;
        Vector RHS;

        // Get state solution
        const auto& rConstCondRef = cond_i;
        rConstCondRef.GetValuesVector(u,0);

        // Get adjoint variables (Corresponds to 1/2*u)
        lambda = 0.5*u;

        // Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
        cond_i.CalculateRightHandSide(RHS, rCurrentProcessInfo);
        for (auto& node_i : cond_i.GetGeometry())
        {
            array_3d gradient_contribution(3, 0.0);
            Vector RHS_perturbed = Vector(RHS.size());
            Vector derived_RHS = Vector(RHS.size());
            
            // x-direction
            node_i.GetInitialPosition()[0] += mDelta;
            node_i.Coordinates()[0] += mDelta;
            cond_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[0] -= mDelta;
            node_i.Coordinates()[0] -= mDelta;
            gradient_contribution[0] = inner_prod(lambda, derived_RHS);


            // y-direction
            node_i.GetInitialPosition()[1] += mDelta;
            node_i.Coordinates()[1] += mDelta;
            cond_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[1] -= mDelta;
            node_i.Coordinates()[1] -= mDelta;
            gradient_contribution[1] = inner_prod(lambda, derived_RHS);

            // z-direction
            node_i.GetInitialPosition()[2] += mDelta;
            node_i.Coordinates()[2] += mDelta;
            cond_i.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
            noalias(derived_RHS) = (RHS_perturbed - RHS) / mDelta;
            node_i.GetInitialPosition()[2] -= mDelta;
            node_i.Coordinates()[2] -= mDelta;
            gradient_contribution[2] = inner_prod(lambda, derived_RHS);

            // Assemble sensitivity to node
            noalias(node_i.FastGetSolutionStepValue(KratosComponents<Variable<array_1d<double,3>>>::Get(shape_gradien_name))) += gradient_contribution;
        }

    };

    void CalculateElementMaterialGradients(Element& elem_i, std::string material_gradien_name, const ProcessInfo &rCurrentProcessInfo){

        // We get the element geometry
        auto& r_this_geometry = elem_i.GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();

        Vector u;
        Vector lambda;
        
        // Get state solution
        const auto& rConstElemRef = elem_i;
        rConstElemRef.GetValuesVector(u,0);

        // Get adjoint variables (Corresponds to 1/2*u)
        lambda = 0.5*u;

        Vector d_RHS_d_E;
        double e_max = elem_i.GetProperties().GetValue(E_MAX);
        double e_min = elem_i.GetProperties().GetValue(E_MIN);
        double e_pr = elem_i.GetProperties().GetValue(E_PR);
        double e_pe = elem_i.GetProperties().GetValue(E_PE);
        elem_i.GetProperties().SetValue(YOUNG_MODULUS,1.0);
        elem_i.CalculateRightHandSide(d_RHS_d_E,rCurrentProcessInfo);
        elem_i.GetProperties().SetValue(YOUNG_MODULUS,e_pe);


        Vector d_RHS_d_D;
        double curr_density = elem_i.GetProperties().GetValue(DENSITY);
        elem_i.GetProperties().SetValue(DENSITY,1.0);
        elem_i.CalculateRightHandSide(d_RHS_d_D,rCurrentProcessInfo);
        elem_i.GetProperties().SetValue(DENSITY,curr_density);


        for (SizeType i_node = 0; i_node < number_of_nodes; ++i_node){
            const auto& d_pe_d_fd = r_this_geometry[i_node].FastGetSolutionStepValue(D_PE_D_FD);
            r_this_geometry[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(material_gradien_name)) += d_pe_d_fd * (e_min + 3 * std::pow((e_pr-e_min)/(e_max-e_min),2.0) * (e_max-e_min)) * inner_prod(d_RHS_d_E,lambda) / number_of_nodes;

            // const auto& d_pd_d_fd = r_this_geometry[i_node].FastGetSolutionStepValue(D_PD_D_FD);
            // r_this_geometry[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get(material_gradien_name)) += d_pd_d_fd * inner_prod(d_RHS_d_D,lambda);
        }

    };        


    void CalculateElementThicknessGradients(Element& elem_i, std::string thickness_gradien_name, const ProcessInfo &rCurrentProcessInfo){
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
            if(rgeom[i_node].Has(PD))
                rValues[i_node] = rgeom[i_node].FastGetSolutionStepValue(PD);  
            else if(this_properties.Has(DENSITY))
                rValues[i_node] = this_properties[DENSITY];
            else
                rValues[i_node] = 1;
        }

        for (SizeType i_node = 0; i_node < num_nodes; ++i_node){
            if(rgeom[i_node].Has(PT))
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
            if(rgeom[i_node].Has(PT))
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
            if(rgeom[i_node].Has(PD))
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
        return "LinearStrainEnergyOptResponse";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "LinearStrainEnergyOptResponse";
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
//      LinearStrainEnergyOptResponse& operator=(LinearStrainEnergyOptResponse const& rOther);

    /// Copy constructor.
//      LinearStrainEnergyOptResponse(LinearStrainEnergyOptResponse const& rOther);


    ///@}

}; // Class LinearStrainEnergyOptResponse

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // LINEAR_STRAIN_ENERGY_OPT_RESPONSE_H
