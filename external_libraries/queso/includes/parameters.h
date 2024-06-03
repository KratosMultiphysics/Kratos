// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef PARAMETERS_INCLUDE_H
#define PARAMETERS_INCLUDE_H

////Project includes
#include "queso/includes/define.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/containers/variant_data_container.hpp"

namespace queso {


///@}
///@name  QuESo Classes
///@{

/**
 * @class  Container to store all global parameters.
 * @author Manuel Messmer
**/
class GlobalParameters : public VariantDataContainer {
public:
    ///@}
    ///@name Typedef's
    ///@{

    typedef VariantDataContainer BaseType;
    using BaseType::Get;
    using BaseType::Set;

    ///@}
    ///@name Life cycle
    ///@{

    /// Default constructor
    GlobalParameters() : VariantDataContainer() {
        AddDefaults();
        CheckComponents();
    }

    /// Constructor for list of components.
    GlobalParameters(BaseType::ComponentVectorType Component) : VariantDataContainer(Component)
    {
        AddDefaults();
        EnsureValidValues();
        CheckComponents();
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Provides default parameters to base class.
    /// @return const BaseType::ComponentVectorType&
    const BaseType::ComponentVectorType& GetDefaultComponents() const override {
        return mDefaultComponents;
    }

    /// @brief Provides aivailable parameteters to base class.
    /// @return const BaseType::AvailableComponentVectorType&
    const BaseType::AvailableComponentVectorType& GetAvailableComponents() const override {
        return mAllAvailableComponents;
    }

    /// @brief Adjust values that are not feasible.
    void EnsureValidValues() override {
        // Make sure this value is not numerically zero.
        const double value = Get<double>("min_element_volume_ratio");
        if( value < 0.9e-10 ){
            BaseType::Set<double>("min_element_volume_ratio", 1e-10);
        }

        std::string input_type = Get<std::string>("input_type");
        if( !(input_type == "stl_file") && !(input_type == "kratos_modelpart") ){
            QuESo_ERROR << "Parameter 'input_type': '" << input_type << "' is not valid. Available options are "
                << "'stl_file' and 'kratos_modelpart'." << std::endl;
        }
    }

private:
    ///@}
    ///@name Private Member Variables
    ///@{

    /// Default components
    inline static const BaseType::ComponentVectorType mDefaultComponents = {
        Component("echo_level", 1UL),
        Component("embedding_flag", true),
        Component("input_type", std::string("stl_file")),
        Component("output_directory_name", std::string("queso_output")),
        Component("initial_triangle_edge_length", 1.0),
        Component("min_num_boundary_triangles", 100UL),
        Component("moment_fitting_residual", 1.0e-10),
        Component("min_element_volume_ratio", 1.0e-3),
        Component("b_spline_mesh", true),
        Component("knot_vector_type", std::string("open_knot_vector")),
        Component("lower_bound_uvw", PointType{-1.0, -1.0, -1.0}),
        Component("upper_bound_uvw", PointType{1.0, 1.0, 1.0}),
        Component("polynomial_order", Vector3i{2UL, 2UL, 2UL} ),
        Component("integration_method", IntegrationMethod::Gauss),
        Component("use_customized_trimmed_points", false) };

    /// All available components
    inline static const BaseType::AvailableComponentVectorType mAllAvailableComponents = {
        std::make_pair<std::string, const std::type_info*>("input_filename", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("output_directory_name", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("input_type", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("input_kratos_modelpart_name", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("postprocess_filename", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("echo_level", &typeid(unsigned long) ),
        std::make_pair<std::string, const std::type_info*>("embedding_flag", &typeid(bool) ),
        std::make_pair<std::string, const std::type_info*>("lower_bound_xyz", &typeid(PointType) ),
        std::make_pair<std::string, const std::type_info*>("upper_bound_xyz", &typeid(PointType) ),
        std::make_pair<std::string, const std::type_info*>("lower_bound_uvw", &typeid(PointType) ),
        std::make_pair<std::string, const std::type_info*>("upper_bound_uvw", &typeid(PointType) ),
        std::make_pair<std::string, const std::type_info*>("b_spline_mesh", &typeid(bool) ),
        std::make_pair<std::string, const std::type_info*>("knot_vector_type", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("polynomial_order", &typeid(Vector3i) ),
        std::make_pair<std::string, const std::type_info*>("number_of_elements", &typeid(Vector3i) ),
        std::make_pair<std::string, const std::type_info*>("initial_triangle_edge_length", &typeid(double) ),
        std::make_pair<std::string, const std::type_info*>("min_num_boundary_triangles", &typeid(unsigned long) ),
        std::make_pair<std::string, const std::type_info*>("moment_fitting_residual", &typeid(double) ),
        std::make_pair<std::string, const std::type_info*>("min_element_volume_ratio", &typeid(double) ),
        std::make_pair<std::string, const std::type_info*>("integration_method", &typeid(IntegrationMethodType) ),
        std::make_pair<std::string, const std::type_info*>("use_customized_trimmed_points", &typeid(bool) ) };
    ///@}
}; // End class GlobalParameters

/**
 * @class  Container to store all condition parameters.
 * @author Manuel Messmer
**/
class ConditionParameters : public VariantDataContainer {
public:
    ///@}
    ///@name Typedef's
    ///@{

    typedef VariantDataContainer BaseType;
    using BaseType::Get;
    using BaseType::Set;

    /// Default constructor
    ConditionParameters(std::string type) : VariantDataContainer() {
        AddDefaults();
        Set<std::string>("type", type);
        CreateAvailableComponentList();
        CheckComponents();
    }

    ///@}
    ///@name Life cycle
    ///@{

    /// Constructor
    ConditionParameters(BaseType::ComponentVectorType Component) : VariantDataContainer(Component)
    {
        AddDefaults();
        CreateAvailableComponentList();
        EnsureValidValues();
        CheckComponents();
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Creates list of available components for all types of Conditions. Must be called in Constructor.
    void CreateAvailableComponentList() {
        if( Get<std::string>("type") == "PenaltySupportCondition" ) {
            mAvailableComponents.insert(mAvailableComponents.end(), mAvailableComponentsPenalty.begin(), mAvailableComponentsPenalty.end());
        } else if( Get<std::string>("type") == "LagrangeSupportCondition" ) {
            mAvailableComponents.insert(mAvailableComponents.end(), mAvailableComponentsLagrange.begin(), mAvailableComponentsLagrange.end());
        } else if( Get<std::string>("type") == "SurfaceLoadCondition" ) {
            mAvailableComponents.insert(mAvailableComponents.end(), mAvailableComponentsSurfaceLoad.begin(), mAvailableComponentsSurfaceLoad.end());
        } else if( Get<std::string>("type") == "PressureLoadCondition" ) {
            mAvailableComponents.insert(mAvailableComponents.end(), mAvailableComponentsPressureLoad.begin(), mAvailableComponentsPressureLoad.end());
        } else {
            QuESo_ERROR << "Condition type '" << Get<std::string>("type") << "' is not available. Available types are: "
                << "'PenaltySupportCondition', 'LagrangeSupportCondition', 'SurfaceLoadCondition', 'PressureLoadCondition'\n";
        }
    }

    /// @brief Provides default parameters to base class.
    /// @return const BaseType::ComponentVectorType&
    const BaseType::ComponentVectorType& GetDefaultComponents() const override {
        return mDefaultComponents;
    }

    /// @brief Provides aivailable parameteters to base class.
    /// @return const BaseType::AvailableComponentVectorType&
    const BaseType::AvailableComponentVectorType& GetAvailableComponents() const override {
        return mAvailableComponents;
    }

    /// @brief Adjust values that are not feasible.
    void EnsureValidValues() override {
        std::string input_type = Get<std::string>("input_type");
        if( !(input_type == "stl_file") && !(input_type == "kratos_modelpart") ){
            QuESo_ERROR << "Parameter 'input_type': " << input_type << " is not valid. Available options are "
                << "'stl_file' and 'kratos_modelpart'" << std::endl;
        }
    }

private:
    ///@}
    ///@name Private Member Variables
    ///@{

    /// Default components.
    inline static const BaseType::ComponentVectorType mDefaultComponents = {
        Component("input_type", std::string("stl_file")) };

    /// Components available in all conditions.
    BaseType::AvailableComponentVectorType mAvailableComponents = {
        std::make_pair<std::string, const std::type_info*>("type", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("input_type", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("input_filename", &typeid(std::string) ),
        std::make_pair<std::string, const std::type_info*>("input_kratos_modelpart_name", &typeid(std::string) ) };

    /// Components available in penalty support conditions.
    inline static const BaseType::AvailableComponentVectorType mAvailableComponentsPenalty = {
        std::make_pair<std::string, const std::type_info*>("penalty_factor", &typeid(double) ),
        std::make_pair<std::string, const std::type_info*>("value", &typeid(PointType) ) };

    /// Components available in lagrange support conditions.
    inline static const BaseType::AvailableComponentVectorType mAvailableComponentsLagrange = {
        std::make_pair<std::string, const std::type_info*>("value", &typeid(PointType) ) };

    /// Components available in surface load conditions.
    inline static const BaseType::AvailableComponentVectorType mAvailableComponentsSurfaceLoad = {
        std::make_pair<std::string, const std::type_info*>("modulus", &typeid(double) ),
        std::make_pair<std::string, const std::type_info*>("direction", &typeid(PointType) ) };

    /// Components available in pressure load conditions.
    inline static const BaseType::AvailableComponentVectorType mAvailableComponentsPressureLoad = {
        std::make_pair<std::string, const std::type_info*>("modulus", &typeid(double) ) };
};

/**
 * @class  Parameters
 * @author Manuel Messmer
 * @brief  Dynamic container for all available parameters. Parameters are split into global parameters and a vector of condition parameters:
 *         one for each condition.
 * @see GlobalParameters and ConditionParameters.
**/
class Parameters {
public:
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    Parameters() {
    }

    /// Constructor for GlobalParameters
    Parameters(const GlobalParameters& rGlobalParameters) : mGlobalParameters(rGlobalParameters)
    {
    }

    /// Constructor for list of Components
    Parameters(std::vector<Component> Components) : mGlobalParameters(Components)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Interface to mGlobalParameters.Get<>()
    /// @tparam type
    /// @param rName
    /// @return onst type&
    template<typename type>
    const type& Get( const std::string& rName ) const {
        return mGlobalParameters.Get<type>(rName);
    }

    /// @brief Interface to mGlobalParameters.Set<>()
    /// @tparam type
    /// @param rName
    /// @param rValue
    template<typename type>
    void Set(const std::string& rName, const type& rValue) {
        mGlobalParameters.Set(rName, rValue);
    }

    /// @brief Interface to mGlobalParameters.Has<>()
    /// @tparam type
    /// @param rName
    /// @return bool
    template<typename type>
    bool Has( const std::string& rName ) const {
        return mGlobalParameters.Has<type>(rName);
    }

    /// @brief Adds global settings to parameter container.
    /// @param rGlobalParameters
    void AddGlobalSettings(const GlobalParameters& rGlobalParameters) {
        mGlobalParameters = rGlobalParameters;
    }

    /// @brief Ad condition settings to Parameters.
    /// @param rConditionParameter
    void AddConditionSettings( const ConditionParameters rConditionParameter ){
        mConditionParameters.push_back(rConditionParameter);
    }

    /// @brief Returns number of conditions. Including Neumann and Dirichlet.
    /// @return IndexType
    IndexType NumberOfConditions(){
        return mConditionParameters.size();
    }

    /// @brief Returns vector of conditions parameter settings.
    /// @return std::vector<ConditionParameters>&
    const std::vector<ConditionParameters> & GetConditionsSettingsVector() const {
        return mConditionParameters;
    }

    /// @brief returns global parameters container.
    /// @return const GlobalParameters&
    const GlobalParameters& GetGlobalSettings() const {
        return mGlobalParameters;
    }

    /// @brief Print all paramters from mGLobalParameters and mConditionParameters: (Name, Value).
    /// @param rOStream
    void PrintInfo(std::ostream& rOStream) const {
        rOStream << "\n----------------- Global Settings -----------------\n";
        mGlobalParameters.PrintInfo(rOStream);
        for( const auto& r_condition_settings : mConditionParameters ){
            rOStream << "--------------- Condition Settings ----------------\n";
            r_condition_settings.PrintInfo(rOStream);
        }
        rOStream << "---------------------------------------------------\n\n";
    }

    /// @brief Checks parameter inputs.
    void Check() const {
        // Orders
        Vector3i order = Order();
        IndexType min_order =  Math::Min( order );
        IndexType max_order =  Math::Max( order );
        QuESo_ERROR_IF(min_order < 1) << "Invalid Input. The polynomial order must be p > 0. \n";
        QuESo_INFO_IF(max_order > 4) << "Warning :: QuESo is designed to construct efficient quadrature rules for 1 <= p <= 4. "
            << "For higher polynomial degrees, the process might become slow. It is recommended to use quadratic bases. Generally, they offer the best performance and accuracy.\n";

        QuESo_INFO_IF(min_order == 1 && EchoLevel() > 0) << "Info :: When using LINEAR finite elements in combination with rather complex geometries, it can be beneficial to employ quadratic quadrature rules within cut elements. "
            << "Linear quadrature rules will converge to the correct solution when using fine discretizations. However, quadratic quadrature rules "
            << "are simply more suited to capture complex cut domains and can hence provide better results for coarse meshes. "
            << "Thus, if you are integrating LINEAR finite elements, consider using '\"polynomial_order\" : [2, 2, 2]' and '\"integration_method\" : \"Gauss_Reduced1\"'. This will generate quadratic quadrature rules in all cut elements "
            << "and linear Gauss rules for all full/interior elements.\n";

        // Number of elements
        Vector3i num_elements = NumberOfElements();
        IndexType tot_num_elements = num_elements[0]*num_elements[1]*num_elements[2];
        QuESo_INFO_IF( tot_num_elements < 2 && EchoLevel() > 0 ) << "You are using only one single element.\n";

        // GGQ rules
        bool ggq_rules_used = GGQRuleIsUsed();
        QuESo_ERROR_IF(ggq_rules_used && max_order > 4) << "Generalized Gauss Quadrature (GGQ) rules are only available for p <= 4.\n";

        bool b_spline_mesh = Get<bool>("b_spline_mesh");
        QuESo_ERROR_IF(ggq_rules_used && !b_spline_mesh) << "Generalized Gauss Quadrature (GGQ) rules are only applicable to B-Spline meshes with C^(p-1) continuity.\n";

        QuESo_ERROR_IF(ggq_rules_used && min_order < 2) << "Generalized Gauss Quadrature (GGQ) rules are only applicable to B-Spline meshes with at least p=2.\n";
    }

    ///////////////////////////////////////////////////
    // Direct Getter Functions for mGlobalParameters //
    ///////////////////////////////////////////////////

    const PointType& LowerBoundXYZ() const {
        return Get<PointType>("lower_bound_xyz");
    }

    const PointType& UpperBoundXYZ() const {
        return Get<PointType>("upper_bound_xyz");
    }

    const PointType& LowerBoundUVW() const {
        return Get<PointType>("lower_bound_uvw");
    }

    const PointType& UpperBoundUVW() const {
        return Get<PointType>("upper_bound_uvw");
    }

    const Vector3i& Order() const {
        return Get<Vector3i>("polynomial_order");
    }

    const IntegrationMethodType& IntegrationMethod() const {
        return Get<IntegrationMethodType>("integration_method");
    }

    const Vector3i& NumberOfElements() const {
        return Get<Vector3i>("number_of_elements");;
    }

    double InitialTriangleEdgeLength() const {
        /// deprecated
        return Get<double>("initial_triangle_edge_length");
    }

    IndexType MinimumNumberOfTriangles() const {
        /// deprecated
        return static_cast<IndexType>( Get<unsigned long>("min_num_boundary_triangles") );
    }

    double MomentFittingResidual() const {
        return Get<double>("moment_fitting_residual");
    }

    IndexType EchoLevel() const {
        return static_cast<IndexType>( Get<unsigned long>("echo_level") );
    }

    bool UseCustomizedTrimmedPositions() const{
        /// Only used for testing.
        return Get<bool>("use_customized_trimmed_points");
    }

    bool GGQRuleIsUsed() const {
        return IntegrationMethod() >= 3;
    }

private:

    ///@}
    ///@name Private Member Variables
    ///@{

    GlobalParameters mGlobalParameters{};
    std::vector<ConditionParameters> mConditionParameters{};

    ///@}

}; // End class Parameters
///@} End QuESo classes

///@}

std::ostream& operator<< (std::ostream& rOStream, const Parameters& rThis);
} // End namespace queso

#endif // PARAMETERS_INCLUDE_H
