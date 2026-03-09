//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Josep Maria Carbonell, Massimiliano Zecchetto
//  Collaborators:
//
//-------------------------------------------------------------
//

#if !defined(KRATOS_UPDATE_CONDITIONS_ON_FREE_SURFACE_PROCESS_H_INCLUDED)
#define KRATOS_UPDATE_CONDITIONS_ON_FREE_SURFACE_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "utilities/builtin_timer.h"
#include "utilities/variable_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class UpdateConditionsOnFreeSurfaceProcess
 * @ingroup PfemFluidDynamicsApplication
 * @brief This process updates the conditions applied to the free surface after the remeshing
 * @author Josep Maria Carbonell, Massimiliano Zecchetto
 * @param rModelPart The PfemFluid model part
 * @param rParameters Parameter object containing the list of sub model parts to be updated and related conditions type
 */
class UpdateConditionsOnFreeSurfaceProcess : public Process {
  public:
	///@name Type Definitions
	///@{

	KRATOS_CLASS_POINTER_DEFINITION(UpdateConditionsOnFreeSurfaceProcess);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	UpdateConditionsOnFreeSurfaceProcess(ModelPart& rModelPart, Parameters rParameters) : mrModelPart(rModelPart) {

		KRATOS_TRY

		// Include only validation with c++11 since raw_literals do not exist in c++03
		Parameters default_parameters(R"(
            {
                "update_conditions"        : true,
                "sub_model_part_list"      : [],
                "reference_condition_list" : []
            }  )");

		// Validate against defaults
		rParameters.ValidateAndAssignDefaults(default_parameters);

		// Each sub model part to be updated must have its reference condition type
		if (rParameters["sub_model_part_list"].size() != rParameters["reference_condition_list"].size()) {
			KRATOS_ERROR << " sub_model_part_list and reference_condition_list must have the same dimension " << std::endl;
		}

		for (unsigned int i = 0; i < rParameters["sub_model_part_list"].size(); i++) {
			mListOfSubModelParts.push_back(rParameters["sub_model_part_list"][i].GetString());
			mListOfReferenceConditions.push_back(rParameters["reference_condition_list"][i].GetString());
		}

		// Set echo level
		mEchoLevel = 0;

		KRATOS_CATCH("");
	}

	/// Destructor.
	virtual ~UpdateConditionsOnFreeSurfaceProcess() {}

	///@}
	///@name Operators
	///@{

	/// This operator is provided to call the process as a function and simply calls the Execute method.
	void operator()() { Execute(); }

	///@}
	///@name Operations
	///@{

	/// Execute method is used to execute the Process algorithms.
	void Execute() override {

		KRATOS_TRY

		this->ClearAllConditions();
		this->UpdateConditionsAfterRemesh();

		KRATOS_CATCH("")
	}

	///@}
	///@name Access
	///@{

	///@}
	///@name Inquiry
	///@{

	///@}
	///@name Input and output
	///@{

	///@}
	///@name Friends
	///@{

	///@}

  private:
	///@name Static Member Variables
	///@{

	///@}
	///@name Static Member Variables
	///@{

	ModelPart& mrModelPart;
	std::vector<std::string> mListOfSubModelParts;
	std::vector<std::string> mListOfReferenceConditions;
	std::size_t mEchoLevel;

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	void UpdateConditionsAfterRemesh() const {

		KRATOS_TRY

		// Get computing model part
		std::string computing_model_part_name;
		for (ModelPart::SubModelPartIterator i_mp = mrModelPart.SubModelPartsBegin(); i_mp != mrModelPart.SubModelPartsEnd(); i_mp++) {
			if ((i_mp->Is(ACTIVE) && i_mp->Is(SOLID)) || (i_mp->Is(ACTIVE) && i_mp->Is(FLUID))) {
				computing_model_part_name = i_mp->Name();
			}
		}
		ModelPart& r_computing_model_part = mrModelPart.GetSubModelPart(computing_model_part_name);

		// Starting free surface condition ID number
		unsigned int condition_id = r_computing_model_part.NumberOfConditions();

		// Loop over free surface sub model parts
		for (unsigned int i = 0; i < mListOfSubModelParts.size(); i++)
		{
			// Free surface sub model part to be updated
			ModelPart& r_sub_model_part = mrModelPart.GetSubModelPart(mListOfSubModelParts[i]);

			// Remove nodes from sub model part
			VariableUtils().SetFlag(TO_ERASE, true, r_sub_model_part.Nodes());
			r_sub_model_part.RemoveNodes(TO_ERASE);
			VariableUtils().SetFlag(TO_ERASE, false, mrModelPart.Nodes()); // this is suspicious

			// Condition type to be created
			const Condition& r_reference_condition = KratosComponents<Condition>::Get(mListOfReferenceConditions[i]);

			// Property to be used in the creation of conditions (assuming r_sub_model_part has only one property)
			Properties::Pointer p_property = mrModelPart.GetParentModelPart().pGetProperties(0);
			for (auto i_prop(r_sub_model_part.PropertiesBegin()); i_prop != r_sub_model_part.PropertiesEnd(); ++i_prop) {
				if (!i_prop->IsEmpty()) {
					const int prop_id = i_prop->Id();
					p_property = r_sub_model_part.GetParentModelPart().pGetProperties(prop_id);
				}
			}

			// For the first sub model part we need to find the boundaries and creating the "skin"
			if ((i == 0))
			{
				BuiltinTimer time_elapsed;

				// Loop over all elements
				for (auto i_elem(r_computing_model_part.ElementsBegin()); i_elem != r_computing_model_part.ElementsEnd(); ++i_elem)
				{
					Geometry<Node>& r_geometry = i_elem->GetGeometry();
					DenseMatrix<unsigned int> lpofa; // connectivities of points defining faces
					DenseVector<unsigned int> lnofa; // number of points defining faces
					r_geometry.NodesInFaces(lpofa);
					r_geometry.NumberNodesInFaces(lnofa);
					ElementWeakPtrVectorType& number_of_neighbour_elements = i_elem->GetValue(NEIGHBOUR_ELEMENTS);

					// Loop over neighbour elements of an element
					unsigned int iface = 0;
					for (auto& i_nelem : number_of_neighbour_elements)
					{
						unsigned int number_of_nodes_in_face = lnofa[iface];

						if (i_nelem.Id() == i_elem->Id())
						{
							// Get face nodes
							Condition::NodesArrayType face_nodes;
							face_nodes.reserve(number_of_nodes_in_face);
							unsigned int number_of_rigid_nodes = 0;
							for (unsigned int j = 1; j <= number_of_nodes_in_face; ++j) {
								face_nodes.push_back(r_geometry(lpofa(j, iface)));
								if (face_nodes[j - 1].Is(RIGID)) {
									number_of_rigid_nodes += 1;
								}
							}

							// Check if it is a free surface
							if (number_of_rigid_nodes != number_of_nodes_in_face)
							{
								// Create new condition on free surface
								Condition::Pointer p_condition = r_reference_condition.Create(condition_id++, face_nodes, p_property);
								mrModelPart.Conditions().push_back(p_condition);
								r_sub_model_part.Conditions().push_back(p_condition);

								// Add nodes to free surface
								for (unsigned int j = 1; j <= number_of_nodes_in_face; ++j) {
									bool add = true;
									for (auto i_node(r_sub_model_part.NodesBegin()); i_node != r_sub_model_part.NodesEnd(); ++i_node)
										if (i_node->Id() == r_geometry(lpofa(j, iface))->Id())
											add = false;
									if (add)
										r_sub_model_part.Nodes().push_back(r_geometry(lpofa(j, iface)));
								}
							}
						}
						iface += 1;
					}
				}
				KRATOS_INFO_IF("UpdateConditionsAfterRemeshProcess ", mEchoLevel > 0)
				    << " SubModelPart: " << mListOfSubModelParts[i] << " Rebuild time: " << time_elapsed.ElapsedSeconds() << std::endl;
			}

			// For the other sub model parts we clone the nodes/conditions from the first sub model part
			else
			{
				BuiltinTimer time_elapsed;

				// Get first free surface sub model part
				ModelPart& r_origin_part = mrModelPart.GetSubModelPart(mListOfSubModelParts[0]);

				// Copy nodes from first free surface sub model part
				r_sub_model_part.AddNodes(r_origin_part.NodesBegin(), r_origin_part.NodesEnd());

				// Initialize vector of new conditions
				ModelPart::ConditionsContainerType temp_conditions;
				temp_conditions.reserve(r_origin_part.NumberOfConditions());

				// Loop over conditions created for first free surface sub model part
				for (auto i_cond = r_origin_part.ConditionsBegin(); i_cond != r_origin_part.ConditionsEnd(); ++i_cond)
				{
					// Create condition (reuse the geometry of the old element to save memory)
					Condition::Pointer p_condition = r_reference_condition.Create(condition_id++, i_cond->pGetGeometry(), p_property);
					temp_conditions.push_back(p_condition);
				}

				// Add conditions
				r_sub_model_part.AddConditions(temp_conditions.begin(), temp_conditions.end());
				mrModelPart.AddConditions(temp_conditions.begin(), temp_conditions.end());

				KRATOS_INFO_IF("UpdateConditionsAfterRemeshProcess ", mEchoLevel > 0)
				    << " SubModelPart: " << mListOfSubModelParts[i] << " Rebuild time: " << time_elapsed.ElapsedSeconds() << std::endl;
			}
		}
		KRATOS_CATCH("");
	}

	void ClearAllConditions() const {
		// Remove conditions only from free surface submodel parts
		for (unsigned int i = 0; i < mListOfSubModelParts.size(); i++)
		{
			ModelPart& r_sub_model_part = mrModelPart.GetSubModelPart(mListOfSubModelParts[i]);
			VariableUtils().SetFlag(TO_ERASE, true, r_sub_model_part.Conditions());
			r_sub_model_part.RemoveConditionsFromAllLevels(TO_ERASE);
		}
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
	UpdateConditionsOnFreeSurfaceProcess& operator=(UpdateConditionsOnFreeSurfaceProcess const& rOther);

	/// this function is a private function

	/// Copy constructor.
	// Process(Process const& rOther);

	///@}

}; // Class Process

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, UpdateConditionsOnFreeSurfaceProcess& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const UpdateConditionsOnFreeSurfaceProcess& rThis) {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_UPDATE_CONDITIONS_ON_FREE_SURFACE_PROCESS_H_INCLUDED  defined
