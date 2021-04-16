#if !defined(KRATOS_BINBASED_INTERPOLATION_UTILITY)
#define  KRATOS_BINBASED_INTERPOLATION_UTILITY

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_flags.h"
#include "geometries/geometry.h"
#include "geometries/triangle_2d_3.h"

//Database includes
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/search/point_point_search.h"

#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

namespace Kratos
{
template <std::size_t TDim>
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) BinBasedInterpolationUtility
{

public:

typedef ModelPart::ElementsContainerType                      ElementsArrayType;
typedef ElementsArrayType::ContainerType                      ResultElementsContainerType;
typedef std::vector<ResultElementsContainerType>              VectorResultElementsContainerType;
typedef ModelPart::ElementsContainerType::iterator            ElementIteratorType;

typedef ModelPart::NodesContainerType                         NodesArrayType;
typedef NodesArrayType::ContainerType                         ResultNodesContainerType;
typedef std::vector<ResultNodesContainerType>                 VectorResultNodesContainerType;
typedef std::vector<Node<3>::Pointer>                         NodalPointersContainerType;
typedef ModelPart::NodesContainerType::iterator               NodeIteratorType;

typedef std::size_t                                           ListIndexType;
typedef SpatialSearch::DistanceType                           DistanceType;
typedef SpatialSearch::VectorDistanceType                     VectorDistanceType;

/// Pointer definition of BinBasedInterpolationUtility
typedef BinBasedInterpolationUtility<TDim> BinBasedInterpolationUtility_TDim;
KRATOS_CLASS_POINTER_DEFINITION(BinBasedInterpolationUtility_TDim);

BinBasedInterpolationUtility()
{ 
    mVariables = VariablesContainer();
}

/// Destructor.
virtual ~BinBasedInterpolationUtility() {}

//AddFluidVariable is not needed here
template<class TDataType>
void AddFluidVariable(Variable<TDataType> const& r_variable){
    std::string variable_list_identifier = "Fluid";
    std::string coupling_variable_description = "Interpolated fluid phase variable";
    AddCouplingVariable<TDataType>(r_variable, variable_list_identifier, coupling_variable_description);
}

template<class TDataType>
void AddVariablesToImposeProjection(Variable<TDataType> const& r_variable){
    std::string variable_list_identifier = "VariableToImpose";
    std::string coupling_variable_description = "Fluid-phase variables with imposed fluid-interpolated values";
    AddCouplingVariable<TDataType>(r_variable, variable_list_identifier, coupling_variable_description);
}

template<class TDataType>
void AddCouplingVariable(Variable<TDataType> const& r_variable, std::string variable_list_identifier, std::string coupling_variable_description){

    if (std::is_same<TDataType, double>::value){
        mVariables.Add(r_variable, variable_list_identifier, "Scalar");
    }

    else if (std::is_same<TDataType, array_1d<double, 3> >::value){
        mVariables.Add(r_variable, variable_list_identifier, "Vector");
    }

    else {
        KRATOS_ERROR << "Variable " << r_variable.Name() << "'s type (" << typeid(r_variable).name()
					 << ") is currently not available as a" << coupling_variable_description << "."
                     << "Please implement." << std::endl;
    }
}

void InterpolateFromFluidMeshCD(ModelPart& reference_model_part, ModelPart& destination_model_part, Parameters& parameters, BinBasedFastPointLocator<TDim>& bin_of_objects_fluid);
/// Turn back information as a stemplate<class T, std::size_t dim> tring.
virtual std::string Info() const
{
    return "";
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
virtual void PrintData(std::ostream& rOStream) const {}

private:

// This class is a simple database that uses two indices (strings) to encode
// different particle lists. The two indices correspond to
// 		i) the set of variables (fluid variables,
//		DEM variables, fluid variables for time filtering etc.) and
// 		ii) the type (Scalar, Vector, etc.)
class VariablesContainer
{
	public:
	std::map<std::set<std::string>, VariablesList> mCouplingVariablesMap;

	VariablesList& GetVariablesList(std::string first_criterion, std::string second_criterion=""){
		std::set<std::string> identifier_set;
		identifier_set.insert("");
		identifier_set.insert(first_criterion);
		identifier_set.insert(second_criterion);
		if (mCouplingVariablesMap.find(identifier_set) == mCouplingVariablesMap.end()){
			VariablesList new_list;
			mCouplingVariablesMap[identifier_set] = new_list;
		}
		return mCouplingVariablesMap[identifier_set];
	}

	// bool Is(VariableData const& rVariable, std::string first_criterion, std::string second_criterion=""){
	// 	return GetVariablesList(first_criterion, second_criterion).Has(rVariable);
	// }

	void Add(VariableData const& rVariable, std::string list_identifier, std::string type_id=""){
		GetVariablesList(list_identifier, type_id).Add(rVariable);
		GetVariablesList(type_id).Add(rVariable);
		GetVariablesList(list_identifier).Add(rVariable);
		GetVariablesList("", "").Add(rVariable);
	}
};

VariablesContainer mVariables;

//***************************************************************************************************************
//***************************************************************************************************************
void Project(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node,  const VariableData *r_destination_variable);
void Interpolate(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<array_1d<double, 3> >& r_origin_variable, const Variable<array_1d<double, 3> >& r_destination_variable);
void Interpolate(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<double>& r_origin_variable, const Variable<double>& r_destination_variable);
void ResetVariables(ModelPart& destination_model_part);
inline void ClearVariable(const NodeIteratorType& node_it, const VariableData& var);

//***************************************************************************************************************
//***************************************************************************************************************

/// Assignment operator.
BinBasedInterpolationUtility& operator=(BinBasedInterpolationUtility const& rOther);

///@}

}; // Class BinBasedInterpolationUtility

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BinBasedInterpolationUtility<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_BINBASED_INTERPOLATION_UTILITY  defined
