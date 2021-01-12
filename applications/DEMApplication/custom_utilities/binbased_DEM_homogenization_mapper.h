//
//   Project Name:        Kratos
//

#if !defined(KRATOS_BINBASED_DEM_HOMOGENIZATION_MAPPER)
#define  KRATOS_BINBASED_DEM_FLUID_HOMOGENIZATION_MAPPER

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_flags.h"
#include "geometries/geometry.h"
#include "geometries/triangle_2d_3.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"

//Database includes
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/search/point_point_search.h"

#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"

// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

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

/// This class allows the interpolation between non-matching simplicial meshes in 2D and 3D with linear shape functions. it is designed for DEM-CFD coupling problems
/** @author  Guillermo Casas Gonzalez <gcasas@cimne.upc.edu>
*
* For every node of the destination model part it is checked in which element of the origin model part it is
* contained and a linear interpolation is performed
*
* The data structure used by default is a bin,
*
* For a more general tool that allows the mapping between 2 and 3D non-matching meshes, please see /kratos/applications/MeshingApplication/custom_utilities/projection.h
*/

//***************************************************************************************************************
//***************************************************************************************************************


template <std::size_t TDim, typename ParticleType>
class KRATOS_API(DEM_APPLICATION) BinBasedDEMHomogenizationMapper
{
public:
///@name Type Definitions
///@{
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

/// Pointer definition of BinBasedDEMHomogenizationMapper
typedef BinBasedDEMHomogenizationMapper<TDim, ParticleType> BinBasedDEMHomogenizationMapper_TDim_TypeOfParticle;
KRATOS_CLASS_POINTER_DEFINITION(BinBasedDEMHomogenizationMapper_TDim_TypeOfParticle);

///@}
///@name Life Cycle
///@{

/// Default constructor.
//--------------------------------------
//  Key for homogenization_type
//--------------------------------------
//   Averaged variables
//--------------------------------------
// 0:   Constant
// 1:   Linear
// 2:   Filtered
//--------------------------------------

BinBasedDEMHomogenizationMapper(Parameters& rParameters)
                             : mMustCalculateMaxNodalArea(true),
                               mLastCouplingFromDEMTime(0.0),
                               mMaxNodalAreaInv(0.0),
                               mNumberOfDEMSamplesSoFarInTheCurrentStep(0)
{
    Parameters default_parameters( R"(
        {
            "homogenization_type": "constant"
        }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);
    mHomogenizationType = mHomogenizationTypeStringToEnumMap[rParameters["homogenization_type"].GetString()];
    //mTimeAveragingType = rParameters["forward_coupling"]["time_averaging_type"].GetInt();

    mGravity = ZeroVector(3);
    mVariables = VariablesContainer();
}

/// Destructor.
virtual ~BinBasedDEMHomogenizationMapper() {}

///@}
///@name Operators
///@{

///@}
///@name Operations
///@{

template<class TDataType>
void AddHomogenizationCouplingVariable(Variable<TDataType> const& r_variable){
    std::string variable_list_identifier = "Homogenization";
    std::string coupling_variable_description = "DEM-interpolated homogenization variable";
    AddCouplingVariable<TDataType>(r_variable, variable_list_identifier, coupling_variable_description);
}

template<class TDataType>
void AddHomogenizationVariableToBeTimeFiltered(Variable<TDataType> const& r_variable, const double time_constant){
	mAlphas[r_variable] = time_constant;
    mIsFirstTimeFiltering[r_variable] = true;
    std::string variable_list_identifier = "HomogenizationTimeFiltered";
    std::string coupling_variable_description = "Homogenization variables to be time-filtered";
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

void InterpolateFromDEMMesh(ModelPart& r_dem_model_part, ModelPart& r_homogenization_model_part, BinBasedFastPointLocator<TDim>& bin_of_objects_homogenization); // this is a bin of objects which contains the FLUID model part
//void VariingRadiusHomogenizeFromDEMMesh(ModelPart& r_dem_model_part, ModelPart& r_homogenization_model_part, const double& search_radius, const double& shape_factor, bool must_search = true);
//void HomogenizeFromDEMMesh(ModelPart& r_dem_model_part, ModelPart& r_homogenization_model_part, const double& search_radius, const double& shape_factor, bool must_search = true);

///@}
///@name Access
///@{

///@}
///@name Inquiry
///@{

///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.
virtual std::string Info() const
{
    return "";
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
virtual void PrintData(std::ostream& rOStream) const {}

///@}
///@name Friends
///@{
DenseVector<unsigned int>& GetElementPartition()
{
  return (mElementsPartition);
}

DenseVector<unsigned int>& GetNodePartition()
{
  return (mNodesPartition);
}

ElementsArrayType::iterator GetElementPartitionBegin(ModelPart& r_model_part, unsigned int k)
{
  ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
  return (pElements.ptr_begin() + this->GetElementPartition()[k]);
}

ElementsArrayType::iterator GetElementPartitionEnd(ModelPart& r_model_part, unsigned int k)
{
  ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
  return (pElements.ptr_begin() + this->GetElementPartition()[k + 1]);
}

NodesArrayType::iterator GetNodePartitionBegin(ModelPart& r_model_part, unsigned int k)
{
  NodesArrayType& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();
  return (pNodes.ptr_begin() + this->GetNodePartition()[k]);
}

NodesArrayType::iterator GetNodePartitionEnd(ModelPart& r_model_part, unsigned int k)
{
  NodesArrayType& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();
  return (pNodes.ptr_begin() + this->GetNodePartition()[k + 1]);
}
///@}

protected:

DenseVector<unsigned int> mElementsPartition;
DenseVector<unsigned int> mNodesPartition;

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

	bool Is(VariableData const& rVariable, std::string first_criterion, std::string second_criterion=""){
		return GetVariablesList(first_criterion, second_criterion).Has(rVariable);
	}

	void Add(VariableData const& rVariable, std::string list_identifier, std::string type_id=""){
		GetVariablesList(list_identifier, type_id).Add(rVariable);
		GetVariablesList(type_id).Add(rVariable);
		GetVariablesList(list_identifier).Add(rVariable);
		GetVariablesList("", "").Add(rVariable);
	}
};

enum HomogenizationType {constant, linear, filtered};
inline static std::map<std::string, HomogenizationType> CreateMap(){
  static std::map<std::string, HomogenizationType> m;
  m["constant"] = constant;
  m["linear"] = linear;
  m["filtered"] = filtered;
  return m;
}
bool mMustCalculateMaxNodalArea;
double mLastCouplingFromDEMTime;
double mMaxNodalAreaInv;
static std::map<std::string, HomogenizationType> mHomogenizationTypeStringToEnumMap;
HomogenizationType mHomogenizationType;
int mTimeAveragingType;
int mNumberOfDEMSamplesSoFarInTheCurrentStep;
array_1d<double, 3> mGravity;

VariablesContainer mVariables;
std::map<VariableData, double> mAlphas;
std::map<VariableData, bool> mIsFirstTimeFiltering;
PointPointSearch::Pointer mpPointPointSearch;

// neighbour lists (for mHomogenizationType = 3)
std::vector<double>  mSearchRadii; // list of nodal search radii (filter radii). It is a vector since spatial search is designed for varying radius
VectorResultNodesContainerType mVectorsOfNeighNodes; // list of arrays of pointers to the particle's nodal neighbours
VectorDistanceType mVectorsOfDistances; // list of arrays of distances to the particle's neighbours
VectorDistanceType mVectorsOfRadii;

//***************************************************************************************************************
//***************************************************************************************************************
void ApplyExponentialTimeFiltering(ModelPart& r_model_part, const VariableData& r_current_variable);
void ApplyExponentialTimeFiltering(ModelPart& r_model_part, const Variable<double>& r_current_variable, const Variable<double>& r_previous_averaged_variable = TIME_AVERAGED_DOUBLE);
void ApplyExponentialTimeFiltering(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& r_current_variable, const Variable<array_1d<double, 3> >& r_previous_averaged_variable = TIME_AVERAGED_ARRAY_3);
void CopyValues(ModelPart& r_model_part, VariableData const& r_origin_variable);
void CopyValues(ModelPart& r_model_part, const Variable<double>& r_origin_variable, const Variable<double>& r_destination_variable = TIME_AVERAGED_DOUBLE);
void CopyValues(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& r_origin_variable, const Variable<array_1d<double, 3> >& r_destination_variable = TIME_AVERAGED_ARRAY_3);
void InterpolateHomogenizationPart(ModelPart& r_dem_model_part, ModelPart& r_homogenization_model_part, BinBasedFastPointLocator<TDim>& bin_of_objects_homogenization); // this is a bin of objects which contains the FLUID model part
void InterpolateOtherHomogenizationVariables(ModelPart& r_dem_model_part, ModelPart& r_homogenization_model_part, BinBasedFastPointLocator<TDim>& bin_of_objects_homogenization); // this is a bin of objects which contains the FLUID model part
void CalculateNodalHomogenizationPartWithLinearWeighing(Element::Pointer p_elem, const Vector& N, ParticleType& particle);
void Distribute(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node,const VariableData *r_destination_variable);
void DistributeDimensionalContributionToHomogenizationPart(Element::Pointer p_elem, const Vector& N, ParticleType& particle);
void CalculatePorosityProjected(ModelPart& r_homogenization_model_part);
void TransferWithConstantWeighing(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<array_1d<double, 3> >& r_destination_variable, const Variable<array_1d<double, 3> >& r_origin_variable);
void TransferWithLinearWeighing(Element::Pointer p_elem, const array_1d<double,TDim + 1>& N, Node<3>::Pointer p_node, const Variable<array_1d<double, 3> >& r_destination_variable, const Variable<array_1d<double, 3> >& r_origin_variable);
void CalculateNodalHomogenizationPartWithConstantWeighing(Element::Pointer p_elem, const Vector& N, ParticleType& particle);
void ResetHomogenizationVariables(ModelPart& r_homogenization_model_part);
void SetToZero(ModelPart& r_model_part, const VariableData& r_variable);
inline void ClearVariable(const NodeIteratorType& node_it, const VariableData& var);
inline unsigned int GetNearestNode(const Vector& N);
double inline GetAlpha(const VariableData& r_variable);
  //***************************************************************************************************************
//***************************************************************************************************************


///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{

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
BinBasedDEMHomogenizationMapper& operator=(BinBasedDEMHomogenizationMapper const& rOther);

///@}

}; // Class BinBasedDEMHomogenizationMapper

/// output stream function
template<std::size_t TDim, typename ParticleType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BinBasedDEMHomogenizationMapper<TDim, ParticleType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_BINBASED_DEM_FLUID_COUPLED_MAPPING  defined
