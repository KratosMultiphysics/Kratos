//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Guillermo Casas, gcasas@cimne.upc.edu$
//   Date:                $Date: 2014-03-08 08:56:42 $
//
//

#if !defined(KRATOS_BINBASED_DEM_FLUID_COUPLED_MAPPING)
#define  KRATOS_BINBASED_DEM_FLUID_COUPLED_MAPPING

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
#include "mollification/density_function_polynomial.h"
#include "custom_functions.h"
#include "fields/velocity_field.h"
#include "fields/fluid_field_utility.h"

//Database includes
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/search/point_point_search.h"

#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"
#include "derivative_recovery.h"

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

// Some function definitions
//***************************************************************************************************************
//***************************************************************************************************************

void ModifyViscosityLikeEinstein(double & viscosity, const double solid_fraction);

void ModifyViscosityLikeLiu(double & viscosity, const double solid_fraction);

//***************************************************************************************************************
//***************************************************************************************************************

template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
class KRATOS_API(SWIMMING_DEM_APPLICATION) BinBasedDEMFluidCoupledMapping
{
public:
///@name Type Definitions
///@{
typedef ModelPart::ElementsContainerType                      ElementsArrayType;
typedef ElementsArrayType::ContainerType                      ResultElementsContainerType;
typedef std::vector<ResultElementsContainerType>              VectorResultElementsContainerType;
typedef ModelPart::ElementsContainerType::iterator            ElementIteratorType;
typedef SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>  ParticleType;

typedef ModelPart::NodesContainerType                         NodesArrayType;
typedef NodesArrayType::ContainerType                         ResultNodesContainerType;
typedef std::vector<ResultNodesContainerType>                 VectorResultNodesContainerType;
typedef std::vector<Node<3>::Pointer>                         NodalPointersContainerType;
typedef ModelPart::NodesContainerType::iterator               NodeIteratorType;

typedef std::size_t                                           ListIndexType;
typedef SpatialSearch::DistanceType                           DistanceType;
typedef SpatialSearch::VectorDistanceType                     VectorDistanceType;

/// Pointer definition of BinBasedDEMFluidCoupledMapping
typedef BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle> BinBasedDEMFluidCoupledMapping_TDim_TBaseTypeOfSwimmingParticle;
KRATOS_CLASS_POINTER_DEFINITION(BinBasedDEMFluidCoupledMapping_TDim_TBaseTypeOfSwimmingParticle);

///@}
///@name Life Cycle
///@{

/// Default constructor.
//----------------------------------------------------------------
//                       Key for coupling_type
//----------------------------------------------------------------
//        Averaged variables       |  Fluid Fraction
//   Fluid-to-DEM | DEM-to-fluid   |
//----------------------------------------------------------------
// 0:   Linear         Constant            Constant
// 1:   Linear         Linear              Constant
// 2:   Linear         Linear              Linear
// 3:   Linear         Filtered            Filtered
//----------------------------------------------------------------


BinBasedDEMFluidCoupledMapping(Parameters& rParameters)
                             : mMustCalculateMaxNodalArea(true),
                               mFluidDeltaTime(0.0),
                               mFluidLastCouplingFromDEMTime(0.0),
                               mMaxNodalAreaInv(0.0),
                               mNumberOfDEMSamplesSoFarInTheCurrentFluidStep(0)
{
    Parameters default_parameters( R"(
        {
            "min_fluid_fraction": 0.2,
            "coupling_type": 1,
            "time_averaging_type": 0,
            "viscosity_modification_type" : 0,
            "n_particles_per_depth_distance" : 1,
            "body_force_per_unit_mass_variable_name" : "BODY_FORCE"
        }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mMinFluidFraction = rParameters["min_fluid_fraction"].GetDouble();
    mCouplingType = rParameters["coupling_type"].GetInt();
    mTimeAveragingType = rParameters["time_averaging_type"].GetInt();
    mViscosityModificationType = rParameters["viscosity_modification_type"].GetInt();
    mParticlesPerDepthDistance = rParameters["n_particles_per_depth_distance"].GetInt();
    mpBodyForcePerUnitMassVariable = &( KratosComponents< Variable<array_1d<double,3>> >::Get(rParameters["body_force_per_unit_mass_variable_name"].GetString()) );

    if (TDim == 3){
        mParticlesPerDepthDistance = 1;
    }

    mGravity = ZeroVector(3);
}

/// Destructor.
virtual ~BinBasedDEMFluidCoupledMapping() {}

///@}
///@name Operators
///@{

///@}
///@name Operations
///@{

void AddDEMCouplingVariable(const VariableData& r_variable){
    mDEMCouplingVariables.Add(r_variable);
}

void AddFluidCouplingVariable(const VariableData& r_variable){
    mFluidCouplingVariables.Add(r_variable);
}

void AddDEMVariablesToImpose(const VariableData& r_variable){
    mDEMVariablesToBeImposed.Add(r_variable);
}

void AddFluidVariableToBeTimeFiltered(const VariableData& r_variable, const double time_constant){
    mFluidVariablesToBeTimeFiltered.Add(r_variable);
    mAlphas[r_variable] = time_constant;
    mIsFirstTimeFiltering[r_variable] = true;
}

void InterpolateFromFluidMesh(ModelPart& r_fluid_model_part, ModelPart& r_dem_model_part, Parameters& parameters, BinBasedFastPointLocator<TDim>& bin_of_objects_fluid, const double alpha);
void ImposeFlowOnDEMFromField(FluidFieldUtility& r_flow, ModelPart& r_dem_model_part);
void ImposeVelocityOnDEMFromFieldToSlipVelocity(FluidFieldUtility& r_flow, ModelPart& r_dem_model_part);
void InterpolateVelocityOnSlipVelocity(ModelPart& r_fluid_model_part, ModelPart& r_dem_model_part, BinBasedFastPointLocator<TDim>& bin_of_objects_fluid);
void UpdateOldVelocity(ModelPart& r_dem_model_part);
void UpdateOldAdditionalForce(ModelPart& r_dem_model_part);
void InterpolateFromNewestFluidMesh(ModelPart& r_fluid_model_part, ModelPart& r_dem_model_part, BinBasedFastPointLocator<TDim>& bin_of_objects_fluid);
void InterpolateFromDEMMesh(ModelPart& r_dem_model_part, ModelPart& r_fluid_model_part, BinBasedFastPointLocator<TDim>& bin_of_objects_fluid); // this is a bin of objects which contains the FLUID model part
void VariingRadiusHomogenizeFromDEMMesh(ModelPart& r_dem_model_part, ModelPart& r_fluid_model_part, const double& search_radius, const double& shape_factor, bool must_search = true);
void HomogenizeFromDEMMesh(ModelPart& r_dem_model_part, ModelPart& r_fluid_model_part, const double& search_radius, const double& shape_factor, bool must_search = true);
void ComputePostProcessResults(ModelPart& r_dem_model_part, ModelPart& r_fluid_model_part, ModelPart& rfem_dem_model_part, BinBasedFastPointLocator<TDim>& bin_of_objects_fluid, const ProcessInfo& r_current_process_info);

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

bool mMustCalculateMaxNodalArea;
double mFluidDeltaTime;
double mFluidLastCouplingFromDEMTime;
double mMinFluidFraction;
double mMaxNodalAreaInv;
int mCouplingType;
int mTimeAveragingType;
int mViscosityModificationType;
int mParticlesPerDepthDistance;
int mNumberOfDEMSamplesSoFarInTheCurrentFluidStep;

array_1d<double, 3> mGravity;
VariablesList mDEMCouplingVariables;
VariablesList mFluidCouplingVariables;
VariablesList mDEMVariablesToBeImposed;
VariablesList mFluidVariablesToBeTimeFiltered;
std::map<VariableData, double> mAlphas;
std::map<VariableData, bool> mIsFirstTimeFiltering;
PointPointSearch::Pointer mpPointPointSearch;

FluidFieldUtility mFlowField;

const Variable<array_1d<double,3>>* mpBodyForcePerUnitMassVariable;

// neighbour lists (for mCouplingType = 3)
std::vector<double>  mSearchRadii; // list of nodal search radii (filter radii). It is a vector since spatial search is designed for varying radius
std::vector<SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>* > mSwimmingSphereElementPointers;
VectorResultNodesContainerType mVectorsOfNeighNodes; // list of arrays of pointers to the particle's nodal neighbours
VectorDistanceType mVectorsOfDistances; // list of arrays of distances to the particle's neighbours
VectorDistanceType mVectorsOfRadii;

//***************************************************************************************************************
//***************************************************************************************************************
void ApplyExponentialTimeFiltering(ModelPart& r_model_part, const VariableData *r_current_variable);
void ApplyExponentialTimeFiltering(ModelPart& r_model_part, const Variable<double>& r_current_variable, const Variable<double>& r_previous_averaged_variable = TIME_AVERAGED_DOUBLE);
void ApplyExponentialTimeFiltering(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& r_current_variable, const Variable<array_1d<double, 3> >& r_previous_averaged_variable = TIME_AVERAGED_ARRAY_3);
void CopyValues(ModelPart& r_model_part, const VariableData *r_origin_variable);
void CopyValues(ModelPart& r_model_part, const Variable<double>& r_origin_variable, const Variable<double>& r_destination_variable = TIME_AVERAGED_DOUBLE);
void CopyValues(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& r_origin_variable, const Variable<array_1d<double, 3> >& r_destination_variable = TIME_AVERAGED_ARRAY_3);
void ComputeHomogenizedFluidFraction(ModelPart& r_fluid_model_part, ModelPart& r_dem_model_part);
void InterpolateFluidFraction(ModelPart& r_dem_model_part, ModelPart& r_fluid_model_part, BinBasedFastPointLocator<TDim>& bin_of_objects_fluid); // this is a bin of objects which contains the FLUID model part
void InterpolateOtherFluidVariables(ModelPart& r_dem_model_part, ModelPart& r_fluid_model_part, BinBasedFastPointLocator<TDim>& bin_of_objects_fluid); // this is a bin of objects which contains the FLUID model part
void SearchParticleNodalNeighbours(ModelPart& r_fluid_model_part, ModelPart& r_dem_model_part, const double& search_radius);
void SearchParticleNodalNeighboursFixedRadius(ModelPart& r_fluid_model_part, ModelPart& r_dem_model_part, const double& search_radius);
void RecalculateDistances(ModelPart& r_dem_model_part);
bool IsDEMVariable(const VariableData& var);
bool IsFluidVariable(const VariableData& var);
bool IsFluidVariableToBeTimeFiltered(const VariableData& var);
array_1d<double, 3> CalculateAcceleration(const Geometry<Node<3> >& geom, const Vector& N);
double CalculateNormOfSymmetricGradient(const Geometry<Node<3> >& geom, const int index);
array_1d<double, 3> CalculateVorticity(const Geometry<Node<3> >& geom, const int index);
void Project(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const VariableData *r_destination_variable);
void Project(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const VariableData *r_destination_variable, double alpha);
void DistributeDimensionalContributionToFluidFraction(Element::Pointer p_elem, const Vector& N, ParticleType& particle);
void Distribute(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node,const VariableData *r_destination_variable);
void ComputeHomogenizedNodalVariable(const ParticleType& particle, const ResultNodesContainerType& neighbours, const DistanceType& weights, const VariableData *r_destination_variable);
void CalculateFluidFraction(ModelPart& r_fluid_model_part);
void CalculateFluidMassFraction(ModelPart& r_fluid_model_part);
void Interpolate(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<array_1d<double, 3> >& r_origin_variable, const Variable<array_1d<double, 3> >& r_destination_variable);
void Interpolate(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<array_1d<double, 3> >& r_origin_variable, const Variable<array_1d<double, 3> >& r_destination_variable, double alpha);
void Interpolate(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<double>& r_origin_variable, const Variable<double>& r_destination_variable);
void Interpolate(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<double>& r_origin_variable, const Variable<double>& r_destination_variable, double alpha);
void CalculateVelocityProjectedRate(Node<3>::Pointer p_node);
void InterpolateAcceleration(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<array_1d<double, 3> >& r_destination_variable);
void InterpolateShearRate(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<double>& r_destination_variable);
void InterpolateShearRate(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<double>& r_destination_variable, double alpha);
void InterpolateVorticity(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<array_1d<double, 3> >& r_destination_variable);
void InterpolateVorticity(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<array_1d<double, 3> >& r_destination_variable, double alpha);
void TransferWithConstantWeighing(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node, const Variable<array_1d<double, 3> >& r_destination_variable, const Variable<array_1d<double, 3> >& r_origin_variable);
void TransferWithLinearWeighing(Element::Pointer p_elem, const array_1d<double,TDim + 1>& N, Node<3>::Pointer p_node, const Variable<array_1d<double, 3> >& r_destination_variable, const Variable<array_1d<double, 3> >& r_origin_variable);
void CalculateNodalFluidFractionWithConstantWeighing(Element::Pointer p_elem, const Vector& N, ParticleType& particle);
void CalculateNodalFluidFractionWithLinearWeighing(Element::Pointer p_elem, const Vector& N, ParticleType& particle);
void CalculateNodalFluidFractionByLumpedL2Projection(Element::Pointer p_elem, const Vector& N, Node<3>::Pointer p_node);
//void CalculateFluidFractionGradient(ModelPart& r_model_part);
void TransferByAveraging(const ParticleType& particle, const ResultNodesContainerType& neighbours, const DistanceType& weights, const Variable<array_1d<double, 3> >& r_destination_variable, const Variable<array_1d<double, 3> >& r_origin_variable);
void CalculateNodalFluidFractionByAveraging(ParticleType& particle, const ResultNodesContainerType& neighbours, const DistanceType& weights);
void CalculateNodalSolidFractionByAveraging(const Node<3>::Pointer p_node, const ResultNodesContainerType& neighbours, const DistanceType& weights, const double averaging_volume_inv);
void MultiplyNodalVariableBy(ModelPart& r_model_part, const Variable<double>& r_variable, const double& factor);
void MultiplyNodalVariableBy(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& r_variable, const double& factor);
void ResetDEMVariables(ModelPart& r_dem_model_part);
void ResetFluidVariables(ModelPart& r_fluid_model_part);
void ResetFLuidVelocityRate(const NodeIteratorType& node_it);
void Clear(ModelPart& r_model_part, const VariableData* r_variable);
void Clear(ModelPart& r_model_part, const Variable<double>& r_variable);
inline void ClearVariable(const NodeIteratorType& node_it, const VariableData *var);
void CalculateFluidNodesMaxNodalArea(ModelPart& r_fluid_model_part);
inline void ClearVariable(const NodeIteratorType& node_it, const VariableData& var);
inline unsigned int GetNearestNode(const Vector& N);
void FillVectorOfSwimmingSpheres(ModelPart& r_dem_model_part);
double inline CalculateDistance(Node<3>::Pointer a, SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>* b);
double inline GetAlpha(const VariableData& r_variable);
const Variable<array_1d<double,3>>& GetBodyForcePerUnitMassVariable() const;
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
BinBasedDEMFluidCoupledMapping& operator=(BinBasedDEMFluidCoupledMapping const& rOther);

///@}

}; // Class BinBasedDEMFluidCoupledMapping

/// output stream function
template<std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_BINBASED_DEM_FLUID_COUPLED_MAPPING  defined
