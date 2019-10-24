//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Me
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_h_process.h"

// Application includes
#include "interface_curvature.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
InterfaceCurvature::InterfaceCurvature(
    ModelPart& rModelPart,
    const bool dummy)
    : Process(),
      mrModelPart(rModelPart)
{
    // Member variables initialization
    mDummyBool = dummy;

    // Initialize the EMBEDDED_IS_ACTIVE variable flag to 0
    this->TestFunction();
}

InterfaceCurvature::InterfaceCurvature(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    //this->CheckDefaultsAndProcessSettings(rParameters);

    // Initialize the EMBEDDED_IS_ACTIVE variable flag to 0
    this->TestFunction();
}

InterfaceCurvature::InterfaceCurvature(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    //this->CheckDefaultsAndProcessSettings(rParameters);

    // Initialize the EMBEDDED_IS_ACTIVE variable flag to 0
    //this->TestFunction();
    int dummy = 0;
}

void InterfaceCurvature::TestFunction()
{
    int iiii = 0;
    for (int i_node = 0; i_node < static_cast<int>(mrModelPart.NumberOfNodes()); ++i_node) {
        auto it_node = mrModelPart.NodesBegin() + i_node;
        it_node->SetValue(EMBEDDED_IS_ACTIVE, 0);
        iiii += 1;
    }

    KRATOS_INFO("Example") << "TESTING INTERFACE CURVATURE, number of node: " << iiii << std::endl;
}

void InterfaceCurvature::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    Parameters default_parameters( R"(
    {
        "model_part_name"                             : ""
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    //mIsModified = false;
    //mDistanceThreshold = rParameters["distance_threshold"].GetDouble();
    //mContinuousDistance = rParameters["continuous_distance"].GetBool();
    //mCheckAtEachStep = rParameters["check_at_each_time_step"].GetBool();
    //mAvoidAlmostEmptyElements = rParameters["avoid_almost_empty_elements"].GetBool();
    //mNegElemDeactivation = rParameters["deactivate_full_negative_elements"].GetBool();
    //mRecoverOriginalDistance = rParameters["recover_original_distance_at_each_step"].GetBool();
    //if (mNegElemDeactivation) {
    //    this->CheckAndStoreVariablesList(rParameters["full_negative_elements_fixed_variables_list"].GetStringArray());
    //}
}

void InterfaceCurvature::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();

    KRATOS_INFO("Example") << "TESTING INTERFACE CURVATURE is now EXECUTED." << std::endl;

    //rOStream << "interface_curvature: NOT WORKING!";
}

void InterfaceCurvature::ExecuteInitialize()
{
    KRATOS_TRY;

    //rOStream << "interface_curvature: NOT WORKING!";

    //if (mContinuousDistance){
        // Continuous distance modification historical variables check
        //const auto& r_node = *mrModelPart.NodesBegin();
        //KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);

        // Compute NODAL_H (used for computing the distance tolerance)
        //FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable> nodal_h_calculator(mrModelPart);
        //nodal_h_calculator.Execute();
    //}

    KRATOS_CATCH("");
}

void InterfaceCurvature::ExecuteBeforeSolutionLoop() {
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void InterfaceCurvature::ExecuteInitializeSolutionStep() {

    //if(!mIsModified){
        // Modify the nodal distance values to avoid bad intersections
        //if (mContinuousDistance) {
            // Modify the continuous distance field
            //this->ModifyDistance();
        //} else {
            // Modify the discontinuous distance field
            //this->ModifyDiscontinuousDistance();
        //}
        //mIsModified = true;

        // If proceeds (depending on the formulation), perform the deactivation
        // Deactivates the full negative elements and sets the inner values to 0
        //if (mNegElemDeactivation) {
        //    this->DeactivateFullNegativeElements();
        //}
    //}
    mDummyBool = false;
}

void InterfaceCurvature::ExecuteFinalizeSolutionStep() {

    // Restore the initial state if the distance is checked at each time step
    //if (mCheckAtEachStep){
        // Restore the is modified flag to false
        //mIsModified = false;
        //if (mNegElemDeactivation){
            // Recover the state previous to the element deactivation
         //   this->RecoverDeactivationPreviousState();
        //}
    //}

    // Recover the original distance values
    //if(mRecoverOriginalDistance) {
        //if (mContinuousDistance){
            //this->RecoverOriginalDistance();
        //} else {
            //this->RecoverOriginalDiscontinuousDistance();
        //}
    //}
    mDummyBool = false;
}

/* Protected functions ****************************************************/

void InterfaceCurvature::SetContinuousDistanceToSplitFlag()
{
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem) {
        auto it_elem = mrModelPart.ElementsBegin() + i_elem;
        auto &r_geom = it_elem->GetGeometry();
        std::vector<double> elem_dist;
        for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
            elem_dist.push_back(r_geom[i_node].FastGetSolutionStepValue(DISTANCE));
            mDummyBool = false;
        }
        //this->SetElementToSplitFlag(*it_elem, elem_dist);
    }
    mDummyBool = true;
}

void InterfaceCurvature::SetDiscontinuousDistanceToSplitFlag()
{
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem) {
        auto it_elem = mrModelPart.ElementsBegin() + i_elem;
        const auto &r_elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);
        //this->SetElementToSplitFlag(*it_elem, r_elem_dist);
    }
}

void InterfaceCurvature::CheckAndStoreVariablesList(const std::vector<std::string>& rVariableStringArray)
{
    //const auto& r_node = *mrModelPart.NodesBegin();
    //for (std::size_t i_variable=0; i_variable < rVariableStringArray.size(); i_variable++){
        //if (KratosComponents<Variable<double>>::Has(rVariableStringArray[i_variable])) {
            //const auto& r_double_var  = KratosComponents<Variable<double>>::Get(rVariableStringArray[i_variable]);
            //KRATOS_CHECK_DOF_IN_NODE(r_double_var, r_node);
            //KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_double_var, r_node)

            //mDoubleVariablesList.push_back(&r_double_var);
        //}
        //else if (KratosComponents<ComponentType>::Has(rVariableStringArray[i_variable])){
            //const auto& r_component_var  = KratosComponents<ComponentType>::Get(rVariableStringArray[i_variable]);
            //KRATOS_CHECK_DOF_IN_NODE(r_component_var, r_node);
            //KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_component_var, r_node)

            //mComponentVariablesList.push_back(&r_component_var);
        //}
        //else if (KratosComponents<Variable<array_1d<double,3>>>::Has(rVariableStringArray[i_variable])){
            //// Checking vector variable in nodal data
            //const auto& r_vector_var = KratosComponents<Variable<array_1d<double,3>>>::Get(rVariableStringArray[i_variable]);
            //KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_vector_var, r_node)

            //// Checking and storing the component variables
            //const auto& r_component_var_x = KratosComponents<ComponentType>::Get(rVariableStringArray[i_variable]+"_X");
            //KRATOS_CHECK_DOF_IN_NODE(r_component_var_x, r_node);
            //mComponentVariablesList.push_back(&r_component_var_x);

            //const auto& r_component_var_y = KratosComponents<ComponentType>::Get(rVariableStringArray[i_variable]+"_Y");
            //KRATOS_CHECK_DOF_IN_NODE(r_component_var_y, r_node);
            //mComponentVariablesList.push_back(&r_component_var_y);

            //if (mrModelPart.GetProcessInfo()[DOMAIN_SIZE] == 3) {
                //const auto& r_component_var_z = KratosComponents<ComponentType>::Get(rVariableStringArray[i_variable]+"_Z");
                //KRATOS_CHECK_DOF_IN_NODE(r_component_var_z, r_node);
                //mComponentVariablesList.push_back(&r_component_var_z);
            //}
        //}
        //else {
            //KRATOS_ERROR << "The variable defined in the list is not a double variable nor a component variable. Given variable: " << rVariableStringArray[i_variable] << std::endl;
        //}
    //}
    mDummyBool = true;
}

}; // namespace Kratos
