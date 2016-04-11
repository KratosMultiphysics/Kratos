//
//   Project Name:        Kratos
//   Last Modified by:    $Author: gcasas $
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
#include "includes/kratos_flags.h"
#include "geometries/geometry.h"
#include "geometries/triangle_2d_3.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "density_function_polynomial.h"
#include "custom_functions.h"

//Database includes
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/point_point_search.h"

#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"
#include "../../DEM_application/DEM_application.h"
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

void ModifyViscosityLikeEinstein(double & viscosity, const double solid_fraction)
{
    viscosity *= (1.0 + 2.5 * solid_fraction);
}

void ModifyViscosityLikeLiu(double & viscosity, const double solid_fraction)
{
    viscosity *= (1.0 + 1.022 * solid_fraction + 1.358 * solid_fraction * solid_fraction * solid_fraction);
}

//***************************************************************************************************************
//***************************************************************************************************************


template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>

class BinBasedDEMFluidCoupledMapping
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

BinBasedDEMFluidCoupledMapping(double min_fluid_fraction,
                               const int coupling_type,
                               const int viscosity_modification_type,
                               const int n_particles_per_depth_distance = 1)
                             : mMustCalculateMaxNodalArea(true),
                               mMinFluidFraction(min_fluid_fraction),
                               mMaxNodalAreaInv(0.0),
                               mCouplingType(coupling_type),
                               mViscosityModificationType(viscosity_modification_type),
                               mParticlesPerDepthDistance(n_particles_per_depth_distance)


{
    if (TDim == 3){
        mParticlesPerDepthDistance = 1;
    }

}

/// Destructor.
virtual ~BinBasedDEMFluidCoupledMapping() {}

///@}
///@name Operators
///@{

///@}
///@name Operations
///@{


//***************************************************************************************************************
//***************************************************************************************************************

void AddDEMCouplingVariable(const VariableData& r_variable){
    mDEMCouplingVariables.Add(r_variable);
}

//***************************************************************************************************************
//***************************************************************************************************************

void AddFluidCouplingVariable(const VariableData& r_variable){
    mFluidCouplingVariables.Add(r_variable);
}

//***************************************************************************************************************
//***************************************************************************************************************
/// Interpolate fluid data onto the DEM model part.
/**
  * @param r_fluid_model_part: the origin model part from which to project
  * @param r_dem_model_part: the destination model part of which we want to interpolate its nodal values
  * @param bin_of_objects_fluid: pre-assembled bin of objects (elelments of the fluid mesh). It is to be constructed separately
  * @see binbased_nodes_in_element_locator
*/
// data_to_project to DEM mesh = alpha * new_data + (1 - alpha) * old_data

void InterpolateFromFluidMesh(
    ModelPart& r_fluid_model_part,
    ModelPart& r_dem_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
    const double alpha)
{
    KRATOS_TRY

    DerivativeRecovery<TDim> recoberer(r_fluid_model_part);
    recoberer.RecoverGradientOfAScalar(PRESSURE, TORQUE);
    // setting interpolated variables to their default values
    ResetDEMVariables(r_dem_model_part);

    array_1d<double, TDim + 1> N;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

    #pragma omp parallel for firstprivate(results, N)
    for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); i++){
        NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
        Node<3>::Pointer p_particle = *(i_particle.base());

        if (p_particle->IsNot(BLOCKED)){
            Element::Pointer p_element;

            // looking for the fluid element in which the DEM node falls
            bool is_found = bin_of_objects_fluid.FindPointOnMesh(p_particle->Coordinates(), N, p_element, results.begin(), max_results);

            // interpolating the variables

            if (is_found){

                p_particle->Set(INSIDE, true);

                for (unsigned int j = 0; j != mDEMCouplingVariables.size(); ++j){
                    Project(p_element, N, p_particle, mDEMCouplingVariables[j], alpha);
                }
            }

            else {
                p_particle->Set(INSIDE, false);
            }
        }
    }

    if (IsDEMVariable(FLUID_ACCEL_PROJECTED)){

      const double delta_time_inv = 1.0 / r_fluid_model_part.GetProcessInfo().GetValue(DELTA_TIME);
      MultiplyNodalVariableBy(r_dem_model_part, FLUID_ACCEL_PROJECTED, delta_time_inv);

    }

    KRATOS_CATCH("")
}

//***************************************************************************************************************
//***************************************************************************************************************
// data_to_project to DEM mesh = current fluid data

void InterpolateFromNewestFluidMesh(
    ModelPart& r_fluid_model_part,
    ModelPart& r_dem_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid)
{
    KRATOS_TRY

    // setting interpolated variables to their default values
    ResetDEMVariables(r_dem_model_part);

    array_1d<double, TDim + 1> N;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

    #pragma omp parallel for firstprivate(results, N)
    for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); i++){
        NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
        Node<3>::Pointer p_particle = *(i_particle.base());

        if (p_particle->IsNot(BLOCKED)){
            Element::Pointer p_element;

            // looking for the fluid element in which the DEM node falls
            bool is_found = bin_of_objects_fluid.FindPointOnMesh(p_particle->Coordinates(), N, p_element, results.begin(), max_results);

            // interpolating variables
            if (is_found){

                for (unsigned int j = 0; j != mDEMCouplingVariables.size(); ++j){
                    Project(p_element, N, p_particle, mDEMCouplingVariables[j]);
                }
            }
        }
    }

    if (IsDEMVariable(FLUID_ACCEL_PROJECTED)){

        const double delta_time_inv = 1.0 / r_fluid_model_part.GetProcessInfo().GetValue(DELTA_TIME);
        MultiplyNodalVariableBy(r_dem_model_part, FLUID_ACCEL_PROJECTED, delta_time_inv);

    }

    KRATOS_CATCH("")
}

//***************************************************************************************************************
//***************************************************************************************************************

/// Interpolate form the DEM  to the fluid mesh
/**
  * @param r_dem_model_part: the origin model part from which to project
  * @param r_fluid_model_part: the destination model part of which we want to interpolate its nodal values
  * @param bin_of_objects_fluid: pre-assembled bin of objects (elelments of the fluid mesh). It is to be constructed separately
  * @see binbased_nodes_in_element_locator
*/

//***************************************************************************************************************
//***************************************************************************************************************
//  data_to_project to fluid mesh = current DEM data

void InterpolateFromDEMMesh(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid) // this is a bin of objects which contains the FLUID model part
{
    KRATOS_TRY

    // setting interpolated variables to their default values
    ResetFluidVariables(r_fluid_model_part);
    // calculating the fluid fraction
    InterpolateFluidFraction(r_dem_model_part, r_fluid_model_part, bin_of_objects_fluid);
    // calculating the rest of fluid variables (particle-fluid force etc.). The solid fraction must be known at this point as it may be used in this step
    InterpolateOtherFluidVariables(r_dem_model_part, r_fluid_model_part, bin_of_objects_fluid);

    KRATOS_CATCH("")
}

//***************************************************************************************************************
//***************************************************************************************************************
//  data_to_project to fluid mesh = current DEM data

void VariingRadiusHomogenizeFromDEMMesh(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    const double& search_radius,
    const double& shape_factor, // it is the density function's maximum divided by its support's radius
    bool must_search = true)
{
    KRATOS_TRY

    // setting interpolated variables to their default values
    ResetFluidVariables(r_fluid_model_part);

    // searching neighbour nodes to each particle (it will have an influence on them only)

    if (mMustCalculateMaxNodalArea){
        CalculateFluidNodesMaxNodalArea(r_fluid_model_part);
    }

    if (must_search){
        SearchParticleNodalNeighbours(r_fluid_model_part, r_dem_model_part, search_radius);
    }

    FillVectorOfSwimmingSpheres(r_dem_model_part);

    if (!must_search) { // we keep the old neighbours
        RecalculateDistances(r_dem_model_part);
    }

    // calculating weights for each particle's nodal contributions

    DensityFunctionPolynomial<3> weighing_function(search_radius, shape_factor);

    #pragma omp parallel for
    for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
        weighing_function.ComputeWeights(mVectorsOfDistances[i], mVectorsOfRadii[i], mMaxNodalAreaInv, mVectorsOfDistances[i]);
    }

    // transferring fluid fraction information onto the fluid (not a naturally parallel task)

    for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
        NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
        CalculateNodalFluidFractionByAveraging(*(i_particle.base()), mSwimmingSphereElementPointers[i]->mNeighbourNodes, mVectorsOfDistances[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < (int)r_fluid_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;
        double area = i_node->GetSolutionStepValue(NODAL_AREA);
        double& fluid_fraction = i_node->FastGetSolutionStepValue(FLUID_FRACTION);
        fluid_fraction = 1.0 - fluid_fraction/area;
    }

    // transferring the rest of effects onto the fluid (not a naturally parallel task)

    for (unsigned int j = 0; j != mFluidCouplingVariables.size(); ++j){

        for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
            NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;

            ComputeHomogenizedNodalVariable(*(i_particle.base()), mSwimmingSphereElementPointers[i]->mNeighbourNodes, mVectorsOfDistances[i], mFluidCouplingVariables[j]);
        }
    }

    KRATOS_CATCH("")
}

//***************************************************************************************************************
//***************************************************************************************************************
//  data_to_project to fluid mesh = current DEM data

void HomogenizeFromDEMMesh(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    const double& search_radius,
    const double& shape_factor, // it is the density function's maximum divided by its support's radius
    bool must_search = true)
{
    KRATOS_TRY

    // setting interpolated variables to their default values
    ResetFluidVariables(r_fluid_model_part);

    // searching neighbour nodes to each particle (it will have an influence on them only)

    if (must_search){
        SearchParticleNodalNeighbours(r_fluid_model_part, r_dem_model_part, search_radius);
    }

    FillVectorOfSwimmingSpheres(r_dem_model_part);

    if (!must_search) { // we keep the old neighbours
        RecalculateDistances(r_dem_model_part);
    }

    // calculating weights for each particle's nodal contributions

    DensityFunctionPolynomial<3> weighing_function(search_radius, shape_factor);

    #pragma omp parallel for
    for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
        weighing_function.ComputeWeights(mVectorsOfDistances[i], mVectorsOfRadii[i], mVectorsOfDistances[i]);
    }

    // transferring fluid fraction information onto the fluid (not a naturally parallel task)

    for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
        NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
        CalculateNodalFluidFractionByAveraging(*(i_particle.base()), mSwimmingSphereElementPointers[i]->mNeighbourNodes, mVectorsOfDistances[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < (int)r_fluid_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;
        double area = i_node->GetSolutionStepValue(NODAL_AREA);
        double& fluid_fraction = i_node->FastGetSolutionStepValue(FLUID_FRACTION);
        fluid_fraction = 1.0 - fluid_fraction / area;
    }

    // transferring the rest of effects onto the fluid (not a naturally parallel task)

    for (unsigned int j = 0; j != mFluidCouplingVariables.size(); ++j){

        for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
            NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;

            ComputeHomogenizedNodalVariable(*(i_particle.base()), mSwimmingSphereElementPointers[i]->mNeighbourNodes, mVectorsOfDistances[i], mFluidCouplingVariables[j]);
        }
    }

    KRATOS_CATCH("")
}

//***************************************************************************************************************
//***************************************************************************************************************

void ComputePostProcessResults(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    ModelPart& rfem_dem_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
    const ProcessInfo& r_current_process_info)
{

    if (IsFluidVariable(FLUID_FRACTION) && mViscosityModificationType){

        void (*modify_viscosity)(double&, const double);
        modify_viscosity = &ModifyViscosityLikeEinstein;

        if (mViscosityModificationType == 2){
            modify_viscosity = &ModifyViscosityLikeLiu;
        }

        else {
            std::cout << "The viscosity modification type " << mViscosityModificationType << " is not supported";
        }

        //#pragma omp parallel for
        for (int i = 0; i < (int)r_fluid_model_part.Nodes().size(); i++){
            NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;
            const double solid_fraction = 1.0 - i_node->FastGetSolutionStepValue(FLUID_FRACTION);
            double& viscosity           = i_node->FastGetSolutionStepValue(VISCOSITY);
            modify_viscosity(viscosity, solid_fraction);
        }
    }

    if (IsDEMVariable(REYNOLDS_NUMBER)){

        //#pragma omp parallel for
        for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); i++){
            ElementIteratorType ielem = r_dem_model_part.ElementsBegin() + i;
            Geometry<Node<3> >& geom = ielem->GetGeometry();
            double& reynolds_number = geom[0].FastGetSolutionStepValue(REYNOLDS_NUMBER);
            ielem->Calculate(REYNOLDS_NUMBER, reynolds_number, r_current_process_info);
        }
    }

    if (IsFluidVariable(AVERAGED_FLUID_VELOCITY)){

        //#pragma omp parallel for
        for (int i = 0; i < (int)r_fluid_model_part.Nodes().size(); i++){
            NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;
            double fluid_fraction                         = i_node->FastGetSolutionStepValue(FLUID_FRACTION);
            const array_1d<double, 3>& darcy_vel          = i_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3>& space_averaged_fluid_vel = i_node->FastGetSolutionStepValue(AVERAGED_FLUID_VELOCITY);
            space_averaged_fluid_vel                      = darcy_vel / fluid_fraction;
        }
    }

    if (IsFluidVariable(SOLID_FRACTION)){

        //#pragma omp parallel for
        for (int i = 0; i < (int)r_fluid_model_part.Nodes().size(); i++){
            NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;
            double& solid_fraction = i_node->FastGetSolutionStepValue(SOLID_FRACTION);
            solid_fraction = 1.0 - i_node->FastGetSolutionStepValue(FLUID_FRACTION);
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************



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
vector<unsigned int>& GetElementPartition()
{
  return (mElementsPartition);
}

vector<unsigned int>& GetNodePartition()
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

vector<unsigned int> mElementsPartition;
vector<unsigned int> mNodesPartition;

private:

bool mMustCalculateMaxNodalArea;
double mMinFluidFraction;
double mMaxNodalAreaInv;
int mCouplingType;
int mViscosityModificationType;
int mParticlesPerDepthDistance;
VariablesList mDEMCouplingVariables;
VariablesList mFluidCouplingVariables;
PointPointSearch::Pointer mpPointPointSearch;

// neighbour lists (for mCouplingType = 3)
std::vector<double>  mSearchRadii; // list of nodal search radii (filter radii). It is a vector since spatial search is designed for varying radius
std::vector<SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>* > mSwimmingSphereElementPointers;
VectorResultNodesContainerType mVectorsOfNeighNodes; // list of arrays of pointers to the particle's nodal neighbours
VectorDistanceType mVectorsOfDistances; // list of arrays of distances to the particle's neighbours
VectorDistanceType mVectorsOfRadii;

//***************************************************************************************************************
//***************************************************************************************************************

void InterpolateFluidFraction(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid) // this is a bin of objects which contains the FLUID model part
{
    array_1d<double, TDim + 1> N;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

    //#pragma omp parallel for firstprivate(results, N)
//    for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); i++){
//        NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
//        Node<3>::Pointer p_particle = *(i_particle.base());
//        Element::Pointer p_element;

//        // looking for the fluid element in which the DEM node falls
//        bool is_found = bin_of_objects_fluid.FindPointOnMesh(p_particle->Coordinates(), N, p_element, results.begin(), max_results);

//        // interpolating variables

//        if (is_found) {
//            DistributeDimensionalContributionToFluidFraction(p_element, N, p_particle);
//        }
//    }

    for (int i = 0; i < (int)r_dem_model_part.Elements().size(); i++){
        ElementIteratorType i_particle = r_dem_model_part.ElementsBegin() + i;
        if (i_particle->GetGeometry()[0].Is(BLOCKED)) {
            continue;
        }
        ParticleType& particle = dynamic_cast<ParticleType&> (*i_particle);
        Element::Pointer p_element;

        // looking for the fluid element in which the DEM node falls
        bool is_found = bin_of_objects_fluid.FindPointOnMesh(particle.GetGeometry()[0].Coordinates(), N, p_element, results.begin(), max_results);

        // interpolating variables

        if (is_found) {
            DistributeDimensionalContributionToFluidFraction(p_element, N, particle);
        }
    }

    CalculateFluidFraction(r_fluid_model_part);

    if (IsFluidVariable(FLUID_FRACTION_GRADIENT)){
        CalculateFluidFractionGradient(r_fluid_model_part);
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void InterpolateOtherFluidVariables(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid) // this is a bin of objects which contains the FLUID model part
{
    // resetting the variables to be mapped
    array_1d<double, TDim + 1> N;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

    //#pragma omp parallel
    {
        for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); i++){
            NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
            Node<3>::Pointer p_particle = *(i_particle.base());
            Element::Pointer p_element;

            // looking for the fluid element in which the DEM node falls
            bool is_found = bin_of_objects_fluid.FindPointOnMesh(p_particle->Coordinates(), N, p_element, results.begin(), max_results);

            // interpolating variables

            if (is_found) {
                //#pragma omp parallel for firstprivate(N)
                for (int j = 0; j < (int)mFluidCouplingVariables.size(); ++j){
                    Distribute(p_element, N, p_particle, mFluidCouplingVariables[j]);
                }
            }
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void SearchParticleNodalNeighbours(ModelPart& r_fluid_model_part,
                                   ModelPart& r_dem_model_part,
                                   const double& search_radius)
{
    unsigned int n_nodes = r_dem_model_part.Nodes().size();

    mSearchRadii.resize(n_nodes, search_radius);

    if (n_nodes != mVectorsOfNeighNodes.size()){
        mVectorsOfNeighNodes.resize(n_nodes);
        mVectorsOfDistances.resize(n_nodes);
        mVectorsOfRadii.resize(n_nodes);
    }

    for (int i = 0; i != (int)n_nodes; ++i){
        mVectorsOfNeighNodes[i].clear();
        mVectorsOfDistances[i].clear();
        mVectorsOfRadii[i].clear();
    }

    NodesArrayType& p_dem_nodes = r_dem_model_part.GetCommunicator().LocalMesh().Nodes();
    NodesArrayType& p_fluid_nodes = r_fluid_model_part.GetCommunicator().LocalMesh().Nodes();
    mpPointPointSearch->SearchPointsImplementation(p_dem_nodes, p_fluid_nodes, mSearchRadii, mVectorsOfNeighNodes, mVectorsOfDistances);

    for (unsigned int i = 0; i < mVectorsOfNeighNodes.size(); i++){
        unsigned int n_neigh = mVectorsOfNeighNodes[i].size();
        mVectorsOfRadii[i].resize(n_neigh);

        for (unsigned int j = 0; j < n_neigh; j++){
            mVectorsOfRadii[i][j] = mVectorsOfNeighNodes[i][j]->FastGetSolutionStepValue(NODAL_AREA);
        }
    }
    // passing the neighbour's information to the particles

    for (int i = 0; i < (int)n_nodes; i++){
        ElementIteratorType i_particle = r_dem_model_part.ElementsBegin() + i;
        SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>* p_particle = dynamic_cast<SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>* >(&(*i_particle));

        if (mVectorsOfNeighNodes[i].size()){
            p_particle->Set(INSIDE, true);
            p_particle->mNeighbourNodes.clear();
            p_particle->mNeighbourNodesDistances.clear();
            p_particle->mNeighbourNodes.insert((p_particle->mNeighbourNodes).begin(), mVectorsOfNeighNodes[i].begin(), mVectorsOfNeighNodes[i].end());
            p_particle->mNeighbourNodesDistances.insert((p_particle->mNeighbourNodesDistances).begin(), mVectorsOfDistances[i].begin(), mVectorsOfDistances[i].end());
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void SearchParticleNodalNeighboursFixedRadius(ModelPart& r_fluid_model_part,
                                   ModelPart& r_dem_model_part,
                                   const double& search_radius)
{
    unsigned int n_nodes = r_dem_model_part.Nodes().size();
    mSearchRadii.resize(n_nodes, search_radius);

    if (n_nodes != mVectorsOfNeighNodes.size()){
        mVectorsOfNeighNodes.resize(n_nodes);
        mVectorsOfDistances.resize(n_nodes);
    }

    for (int i = 0; i != (int)n_nodes; ++i){
        mVectorsOfNeighNodes[i].clear();
        mVectorsOfDistances[i].clear();
    }

    NodesArrayType& p_dem_nodes = r_dem_model_part.GetCommunicator().LocalMesh().Nodes();
    NodesArrayType& p_fluid_nodes = r_fluid_model_part.GetCommunicator().LocalMesh().Nodes();
    mpPointPointSearch->SearchPointsImplementation(p_dem_nodes, p_fluid_nodes, mSearchRadii, mVectorsOfNeighNodes, mVectorsOfDistances);

    // passing the neighbour's information to the particles

    for (int i = 0; i < (int)n_nodes; i++){
        ElementIteratorType i_particle = r_dem_model_part.ElementsBegin() + i;
        SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>* p_particle = dynamic_cast<SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>* >(&(*i_particle));

        if (mVectorsOfNeighNodes[i].size()){
            p_particle->Set(INSIDE, true);
            p_particle->mNeighbourNodes.clear();
            p_particle->mNeighbourNodesDistances.clear();
            p_particle->mNeighbourNodes.insert((p_particle->mNeighbourNodes).begin(), mVectorsOfNeighNodes[i].begin(), mVectorsOfNeighNodes[i].end());
            p_particle->mNeighbourNodesDistances.insert((p_particle->mNeighbourNodesDistances).begin(), mVectorsOfDistances[i].begin(), mVectorsOfDistances[i].end());
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void RecalculateDistances(ModelPart& r_dem_model_part){

    int n_particles = (int)mSwimmingSphereElementPointers.size();
    mVectorsOfDistances.resize(n_particles);
    mVectorsOfRadii.resize(n_particles);

    for (int i = 0; i != n_particles; ++i){
        SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>* p_particle = mSwimmingSphereElementPointers[i];
        int n_neighbours = (int)p_particle->mNeighbourNodes.size();
        mVectorsOfDistances[i].resize(n_neighbours);
        mVectorsOfRadii[i].resize(n_neighbours);

        for (int j = 0; j != n_neighbours; ++j){

            mVectorsOfDistances[i][j] = CalculateDistance(p_particle->mNeighbourNodes[j], p_particle);
            // spatial search is designed for varying radius
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

bool IsDEMVariable(const VariableData& var)
{
    for (unsigned int i = 0; i != mDEMCouplingVariables.size(); ++i){

        if (*mDEMCouplingVariables[i] == var){
            return true;
        }
    }

    return false;
}

//***************************************************************************************************************
//***************************************************************************************************************

bool IsFluidVariable(const VariableData& var)
{
    for (unsigned int i = 0; i != mFluidCouplingVariables.size(); ++i){

        if (*mFluidCouplingVariables[i] == var){
            return true;
        }
    }

    return false;
}

//***************************************************************************************************************
//***************************************************************************************************************

array_1d<double, 3> CalculateAcceleration(const Geometry<Node<3> >& geom, const array_1d<double, TDim + 1>& N)
{
    array_1d<double, 3> acceleration = ZeroVector(3);

    for (unsigned int i = 0; i < TDim + 1; ++i){
        const array_1d<double, 3>& vel_old = geom[i].FastGetSolutionStepValue(VELOCITY, 1);
        const array_1d<double, 3>& vel_new = geom[i].FastGetSolutionStepValue(VELOCITY);
        acceleration += N[i] * (vel_new - vel_old);
    }

    return(acceleration);
}

//***************************************************************************************************************
//***************************************************************************************************************

double CalculateNormOfSymmetricGradient(const Geometry<Node<3> >& geom, const int index)
{
    Geometry<Node<3> >::ShapeFunctionsGradientsType DN_DX;

    // calculating the gradient of the shape functions on the Gauss points (its ok, since their value is constant over the element)
    geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::GI_GAUSS_1);

    Matrix S = ZeroMatrix(TDim, TDim);
    const unsigned int n_nodes = geom.PointsNumber();

    for (unsigned int n = 0; n < n_nodes; ++n){
        const array_1d<double, 3>& vel = geom[n].FastGetSolutionStepValue(VELOCITY, index);

        for (unsigned int i = 0; i < TDim; ++i){

            for (unsigned int j = 0; j < TDim; ++j){
                S(i, j) += 0.5 * (DN_DX[0](n, j) * vel[i] + DN_DX[0](n, i) * vel[j]);
            }
        }
    }

    // norm of the symetric gradient (shear rate)
    double norm_s = 0.0;

    for (unsigned int i = 0; i < TDim; ++i){

        for (unsigned int j = 0; j < TDim; ++j){
            norm_s += S(i, j) * S(i, j);
        }
    }

    norm_s = sqrt(2.0 * norm_s);

    return(norm_s);
}

//***************************************************************************************************************
//***************************************************************************************************************

array_1d<double, 3> CalculateVorticity(const Geometry<Node<3> >& geom, const int index)
{
    Geometry<Node<3> >::ShapeFunctionsGradientsType DN_DX;

    // calculating the gradient of the shape functions on the Gauss points (its ok, since their value is constant over the element)
    geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::GI_GAUSS_1);

    array_1d<double, 3> vorticity = ZeroVector(3);
    array_1d<double, 3> derivatives = ZeroVector(3);

    const unsigned int n_nodes = geom.PointsNumber();

    for (unsigned int n = 0; n < n_nodes; ++n){

        for (unsigned int i = 0; i < TDim; ++i){
            derivatives[i] = DN_DX[0](n, i);
        }

        const array_1d<double, 3>& vel = geom[n].FastGetSolutionStepValue(VELOCITY, index);
        vorticity += MathUtils<double>::CrossProduct(derivatives, vel);
    }

    return(vorticity);
}

//***************************************************************************************************************
//***************************************************************************************************************

void Project(Element::Pointer p_elem,
             const array_1d<double, TDim + 1>& N,
             Node<3>::Pointer p_node,
             const VariableData *r_destination_variable)
{
    if (*r_destination_variable == FLUID_DENSITY_PROJECTED){
        Interpolate(p_elem, N, p_node, DENSITY, FLUID_DENSITY_PROJECTED);
    }

    else if (*r_destination_variable == FLUID_FRACTION_PROJECTED && IsFluidVariable(FLUID_FRACTION)){
        Interpolate(p_elem, N, p_node, FLUID_FRACTION, FLUID_FRACTION_PROJECTED);
    }

    else if (*r_destination_variable == FLUID_VEL_PROJECTED){
        Interpolate(p_elem, N, p_node, VELOCITY, FLUID_VEL_PROJECTED);
    }

    else if (*r_destination_variable == FLUID_VEL_LAPL_PROJECTED){
        Interpolate(p_elem, N, p_node, VELOCITY_LAPLACIAN, FLUID_VEL_LAPL_PROJECTED);
    }

    else if (*r_destination_variable == FLUID_VEL_LAPL_RATE_PROJECTED){
        Interpolate(p_elem, N, p_node, VELOCITY_LAPLACIAN_RATE, FLUID_VEL_LAPL_RATE_PROJECTED);
    }

    else if (*r_destination_variable == PRESSURE_GRAD_PROJECTED){
        Interpolate(p_elem, N, p_node, PRESSURE_GRADIENT, PRESSURE_GRAD_PROJECTED);
    }

    else if (*r_destination_variable == FLUID_VISCOSITY_PROJECTED){
        Interpolate(p_elem, N, p_node, VISCOSITY, FLUID_VISCOSITY_PROJECTED);
    }

    else if (*r_destination_variable == POWER_LAW_N){
        Interpolate(p_elem, N, p_node, POWER_LAW_N, POWER_LAW_N);
    }

    else if (*r_destination_variable == POWER_LAW_K){
        Interpolate(p_elem, N, p_node, POWER_LAW_K, POWER_LAW_K);
    }

    else if (*r_destination_variable == YIELD_STRESS){
        Interpolate(p_elem, N, p_node, YIELD_STRESS, YIELD_STRESS);
    }

    else if (*r_destination_variable == DISTANCE){
        Interpolate(p_elem, N, p_node, DISTANCE, DISTANCE);
    }

    else if (*r_destination_variable == FLUID_ACCEL_PROJECTED){
        InterpolateAcceleration(p_elem, N, p_node, FLUID_ACCEL_PROJECTED);
    }

    else if (*r_destination_variable == SHEAR_RATE_PROJECTED){
        InterpolateShearRate(p_elem, N, p_node, SHEAR_RATE_PROJECTED);
    }

    else if (*r_destination_variable == FLUID_VORTICITY_PROJECTED){
        InterpolateVorticity(p_elem, N, p_node, FLUID_VORTICITY_PROJECTED);
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void Project(Element::Pointer p_elem,
             const array_1d<double, TDim + 1> N,
             Node<3>::Pointer p_node,
             const VariableData *r_destination_variable,
             double alpha)
{
    if (*r_destination_variable == FLUID_DENSITY_PROJECTED){
        Interpolate(p_elem, N, p_node, DENSITY, FLUID_DENSITY_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_FRACTION_PROJECTED && IsFluidVariable(FLUID_FRACTION)){
        Interpolate(p_elem, N, p_node, FLUID_FRACTION, FLUID_FRACTION_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_VEL_PROJECTED){
        Interpolate(p_elem, N, p_node, VELOCITY, FLUID_VEL_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_VEL_LAPL_PROJECTED){
        Interpolate(p_elem, N, p_node, VELOCITY_LAPLACIAN, FLUID_VEL_LAPL_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_VEL_LAPL_RATE_PROJECTED){
        Interpolate(p_elem, N, p_node, VELOCITY_LAPLACIAN_RATE, FLUID_VEL_LAPL_RATE_PROJECTED, alpha);
    }

    else if (*r_destination_variable == PRESSURE_GRAD_PROJECTED){
        Interpolate(p_elem, N, p_node, PRESSURE_GRADIENT, PRESSURE_GRAD_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_VISCOSITY_PROJECTED){
        Interpolate(p_elem, N, p_node, VISCOSITY, FLUID_VISCOSITY_PROJECTED, alpha);
    }

    else if (*r_destination_variable == POWER_LAW_N){
        Interpolate(p_elem, N, p_node, POWER_LAW_N, POWER_LAW_N, alpha);
    }

    else if (*r_destination_variable == POWER_LAW_K){
        Interpolate(p_elem, N, p_node, POWER_LAW_K, POWER_LAW_K, alpha);
    }

    else if (*r_destination_variable == YIELD_STRESS){
        Interpolate(p_elem, N, p_node, YIELD_STRESS, YIELD_STRESS, alpha);
    }

    else if (*r_destination_variable == DISTANCE){
        Interpolate(p_elem, N, p_node, DISTANCE, DISTANCE, alpha);
    }

    else if (*r_destination_variable == FLUID_ACCEL_PROJECTED){
        InterpolateAcceleration(p_elem, N, p_node, FLUID_ACCEL_PROJECTED);
    }

    else if (*r_destination_variable == SHEAR_RATE_PROJECTED){
        InterpolateShearRate(p_elem, N, p_node, SHEAR_RATE_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_VORTICITY_PROJECTED){
        InterpolateVorticity(p_elem, N, p_node, FLUID_VORTICITY_PROJECTED, alpha);
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

//void DistributeDimensionalContributionToFluidFraction(
//    Element::Pointer p_elem,
//    const array_1d<double, TDim + 1>& N,
//    Node<3>::Pointer p_node)
//{
//    if (mCouplingType == 0 || mCouplingType == 1){
//        CalculateNodalFluidFractionWithConstantWeighing(p_elem, N, p_node);
//    }

//    else if (mCouplingType == 2){
//        CalculateNodalFluidFractionWithLinearWeighing(p_elem, N, p_node);
//    }

//    else if (mCouplingType == 4){
//        CalculateNodalFluidFractionByLumpedL2Projection(p_elem, N, p_node);
//    }
//}

void DistributeDimensionalContributionToFluidFraction(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    ParticleType& particle)
{
    if (mCouplingType == 0 || mCouplingType == 1){
        CalculateNodalFluidFractionWithConstantWeighing(p_elem, N, particle);
    }

    else if (mCouplingType == 2){
        CalculateNodalFluidFractionWithLinearWeighing(p_elem, N, particle);
    }

//    else if (mCouplingType == 4){
//        CalculateNodalFluidFractionByLumpedL2Projection(p_elem, N, p_particle);
//    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void Distribute(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    const VariableData *r_destination_variable)
{
    if (mCouplingType == 0){

        if (*r_destination_variable == BODY_FORCE){
            TransferWithConstantWeighing(p_elem, N, p_node, BODY_FORCE, HYDRODYNAMIC_FORCE);
        }
    }

    else if (mCouplingType == 1){

        if (*r_destination_variable == BODY_FORCE){
            TransferWithLinearWeighing(p_elem, N, p_node, BODY_FORCE, HYDRODYNAMIC_FORCE);
        }
    }

    else if (mCouplingType == 2){

        if (*r_destination_variable == BODY_FORCE){
            TransferWithLinearWeighing(p_elem, N, p_node, BODY_FORCE, HYDRODYNAMIC_FORCE);
        }
    }

    else if (mCouplingType == - 1){

        if (*r_destination_variable == BODY_FORCE){
            TransferWithLinearWeighing(p_elem, N, p_node, BODY_FORCE, HYDRODYNAMIC_FORCE);
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void ComputeHomogenizedNodalVariable(
    const Node<3>::Pointer p_node,
    const ResultNodesContainerType& neighbours,
    const DistanceType& weights,
    const VariableData *r_destination_variable)
{
    if (*r_destination_variable == BODY_FORCE){
        TransferByAveraging(p_node, neighbours, weights, BODY_FORCE, HYDRODYNAMIC_FORCE);
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void  CalculateFluidFraction(ModelPart& r_fluid_model_part)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_fluid_model_part.Nodes().size(), mNodesPartition);

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

        for (NodesArrayType::iterator i_node = this->GetNodePartitionBegin(r_fluid_model_part, k); i_node != this->GetNodePartitionEnd(r_fluid_model_part, k); ++i_node){
            double& fluid_fraction = i_node->FastGetSolutionStepValue(FLUID_FRACTION);

            if (mCouplingType != 4){
                fluid_fraction = 1.0 - fluid_fraction / i_node->FastGetSolutionStepValue(NODAL_AREA);
            }

            else {
                fluid_fraction = 1.0 - fluid_fraction;
            }

            if (fluid_fraction < mMinFluidFraction){
                fluid_fraction = mMinFluidFraction;
            }
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
// project an array1D

void Interpolate(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<array_1d<double, 3> >& r_origin_variable,
    Variable<array_1d<double, 3> >& r_destination_variable)
{
    // geometry of the element of the origin model part   
    Geometry<Node<3> >& geom = p_elem->GetGeometry();
    double N_fast[TDim + 1];
    N_fast[TDim] = 1.0;

    for (unsigned int i = 0; i < TDim; i++){
        N_fast[i] = N[i];
        N_fast[TDim] -= N_fast[i];
    }

    array_1d<double, 3>& first_origin_datum = geom[0].FastGetSolutionStepValue(r_origin_variable);
    double step_datum_fast[TDim];

    for (unsigned int j = 0; j < TDim; j++){
        step_datum_fast[j] = N_fast[0] * first_origin_datum[j];
    }

    for (unsigned int i = 1; i < TDim + 1; i++){
        array_1d<double, 3>& origin_datum = geom[i].FastGetSolutionStepValue(r_origin_variable);

        for (unsigned int j = 0; j < TDim; j++){
            step_datum_fast[j] += N_fast[i] * origin_datum[j];
        }
    }

    array_1d<double, 3>& step_datum = p_node->FastGetSolutionStepValue(r_destination_variable);

    for (unsigned int j = 0; j < TDim; j++){
        step_datum[j] = step_datum_fast[j];
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void Interpolate(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<array_1d<double, 3> >& r_origin_variable,
    Variable<array_1d<double, 3> >& r_destination_variable,
    double alpha)
{
    // Geometry of the element of the origin model part
    Geometry<Node<3> >& geom = p_elem->GetGeometry();
    double N_fast[TDim + 1];
    N_fast[TDim] = 1.0;

    for (unsigned int i = 0; i < TDim; i++){
        N_fast[i] = N[i];
        N_fast[TDim] -= N_fast[i];
    }

    array_1d<double, 3>& first_origin_datum = geom[0].FastGetSolutionStepValue(r_origin_variable);
    array_1d<double, 3>& first_origin_datum_old = geom[0].FastGetSolutionStepValue(r_origin_variable, 1);
    double step_datum_fast[TDim];

    for (unsigned int j = 0; j < TDim; j++){
        step_datum_fast[j] = N_fast[0] * (alpha * first_origin_datum[j] + (1.0 - alpha) * first_origin_datum_old[j]);
    }
    // Destination data

    for (unsigned int i = 1; i < TDim + 1; i++){
        array_1d<double, 3>& origin_datum = geom[i].FastGetSolutionStepValue(r_origin_variable);
        array_1d<double, 3>& origin_datum_old = geom[i].FastGetSolutionStepValue(r_origin_variable, 1);

        for (unsigned int j = 0; j < TDim; j++){
            step_datum_fast[j] += N_fast[i] * (alpha * origin_datum[j] + (1.0 - alpha) * origin_datum_old[j]);
        }
    }

    array_1d<double, 3>& step_datum = p_node->FastGetSolutionStepValue(r_destination_variable);

    for (unsigned int i = 0; i < TDim; i++){
        step_datum[i] = step_datum_fast[i];
    }

 }

//***************************************************************************************************************
//***************************************************************************************************************
// projecting a scalar

void Interpolate(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<double>& r_origin_variable,
    Variable<double>& r_destination_variable)
{
    // Geometry of the element of the origin model part
    Geometry<Node<3> >& geom = p_elem->GetGeometry();

    // Destination data
    double& step_data = p_node->FastGetSolutionStepValue(r_destination_variable);

    step_data = N[0] * geom[0].FastGetSolutionStepValue(r_origin_variable);

    for (unsigned int i = 1; i < TDim + 1; i++){
        step_data += N[i] * geom[i].FastGetSolutionStepValue(r_origin_variable);
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void Interpolate(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<double>& r_origin_variable,
    Variable<double>& r_destination_variable,
    double alpha)
{
    // Geometry of the element of the origin model part
    Geometry<Node<3> >& geom = p_elem->GetGeometry();

    // Destination data
    double& step_data = p_node->FastGetSolutionStepValue(r_destination_variable);
    step_data += N[0] * (alpha * geom[0].FastGetSolutionStepValue(r_origin_variable) + (1 - alpha) * geom[0].FastGetSolutionStepValue(r_origin_variable, 1));

    for (unsigned int i = 1; i < TDim + 1; i++){
      step_data += N[i] * (alpha * geom[i].FastGetSolutionStepValue(r_origin_variable) + (1 - alpha) * geom[i].FastGetSolutionStepValue(r_origin_variable, 1));
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void InterpolateAcceleration(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<array_1d<double, 3> >& r_destination_variable)
{
    p_node->FastGetSolutionStepValue(r_destination_variable) = CalculateAcceleration(p_elem->GetGeometry(), N);
}

//***************************************************************************************************************
//***************************************************************************************************************

void InterpolateShearRate(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<double>& r_destination_variable)
{
    double shear_rate = CalculateNormOfSymmetricGradient(p_elem->GetGeometry(), 0);
    double& step_data = p_node->FastGetSolutionStepValue(r_destination_variable);

    step_data = shear_rate;
}

//***************************************************************************************************************
//***************************************************************************************************************

void InterpolateShearRate(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<double>& r_destination_variable,
    double alpha)
{
    double shear_rate       = CalculateNormOfSymmetricGradient(p_elem->GetGeometry(), 0);
    double prev_shear_rate  = CalculateNormOfSymmetricGradient(p_elem->GetGeometry(), 1);
    double& step_data       = p_node->FastGetSolutionStepValue(r_destination_variable);

    step_data = alpha * shear_rate + (1.0 - alpha) * prev_shear_rate;
}

//***************************************************************************************************************
//***************************************************************************************************************

void InterpolateVorticity(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<array_1d<double, 3> >& r_destination_variable)
{
    array_1d<double, 3> vorticity  = CalculateVorticity(p_elem->GetGeometry(), 0);
    array_1d<double, 3>& step_data = p_node->FastGetSolutionStepValue(r_destination_variable);

    step_data = vorticity;
}

//***************************************************************************************************************
//***************************************************************************************************************

void InterpolateVorticity(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<array_1d<double, 3> >& r_destination_variable,
    double alpha)
{
    array_1d<double, 3> vorticity      = CalculateVorticity(p_elem->GetGeometry(), 0);
    array_1d<double, 3> prev_vorticity = CalculateVorticity(p_elem->GetGeometry(), 1);
    array_1d<double, 3>& step_data     = p_node->FastGetSolutionStepValue(r_destination_variable);

    step_data = alpha * vorticity + (1.0 - alpha) * prev_vorticity;
}

//***************************************************************************************************************
//***************************************************************************************************************

void TransferWithConstantWeighing(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<array_1d<double, 3> >& r_destination_variable,
    Variable<array_1d<double, 3> >& r_origin_variable)
{
    // Geometry of the element of the destination model part
    Geometry<Node<3> >& geom              = p_elem->GetGeometry();
    unsigned int i_nearest_node            = GetNearestNode(N);
    const array_1d<double, 3>& origin_data = p_node->FastGetSolutionStepValue(r_origin_variable);
    array_1d<double, 3>& destination_data  = geom[i_nearest_node].FastGetSolutionStepValue(r_destination_variable);

    if (r_origin_variable == HYDRODYNAMIC_FORCE){
        const double fluid_fraction = geom[i_nearest_node].FastGetSolutionStepValue(FLUID_FRACTION);
        const double nodal_volume   = geom[i_nearest_node].FastGetSolutionStepValue(NODAL_AREA);
        const double density        = geom[i_nearest_node].FastGetSolutionStepValue(DENSITY);
        const double nodal_mass_inv = mParticlesPerDepthDistance / (fluid_fraction * density * nodal_volume);
        destination_data            = - nodal_mass_inv * origin_data;
    }

    else {
        std::cout << "Variable " << r_origin_variable << " is not supported for transference with constant weights";
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void TransferWithLinearWeighing(
    Element::Pointer p_elem,
    const array_1d<double,TDim + 1>& N,
    Node<3>::Pointer p_node,
    Variable<array_1d<double, 3> >& r_destination_variable,
    Variable<array_1d<double, 3> >& r_origin_variable)
{
    // Geometry of the element of the destination model part
    Geometry<Node<3> >& geom = p_elem->GetGeometry();
    const array_1d<double, 3>& origin_data = p_node->FastGetSolutionStepValue(r_origin_variable);

    if (r_origin_variable == HYDRODYNAMIC_FORCE){

        for (unsigned int i = 0; i < TDim + 1; i++){
            array_1d<double, 3>& hydrodynamic_reaction = geom[i].FastGetSolutionStepValue(HYDRODYNAMIC_REACTION);
            array_1d<double, 3>& body_force            = geom[i].FastGetSolutionStepValue(BODY_FORCE);
            const double& fluid_fraction               = geom[i].FastGetSolutionStepValue(FLUID_FRACTION);
            const double& nodal_volume                 = geom[i].FastGetSolutionStepValue(NODAL_AREA);
            const double& density                      = geom[i].FastGetSolutionStepValue(DENSITY);
            hydrodynamic_reaction -= mParticlesPerDepthDistance * N[i] / (fluid_fraction * density * nodal_volume) * origin_data;
            body_force = hydrodynamic_reaction;
            body_force[2] -= 9.81;
        }
    }

    else {
        std::cout << "Variable " << r_origin_variable << " is not supported for transference with linear weights" ;
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateNodalFluidFractionWithConstantWeighing(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    ParticleType& particle)
{
    unsigned int i_nearest_node = GetNearestNode(N);

    // Geometry of the element of the destination model part
    const double particle_volume = particle.CalculateVolume();
    p_elem->GetGeometry()[i_nearest_node].FastGetSolutionStepValue(FLUID_FRACTION) += particle_volume;
}

//***************************************************************************************************************
//***************************************************************************************************************

//void CalculateNodalFluidFractionWithLinearWeighing(
//    Element::Pointer p_elem,
//    const array_1d<double, TDim + 1>& N,
//    Node<3>::Pointer p_node)
//{
//    // Geometry of the element of the origin model part
//    Geometry<Node<3> >& geom = p_elem->GetGeometry();

//    const double& radius         = p_node->FastGetSolutionStepValue(RADIUS);
//    const double particle_volume = 1.33333333333333333333 * KRATOS_M_PI * mParticlesPerDepthDistance * radius * radius * radius;

//    for (unsigned int i = 0; i < TDim + 1; i++){
//        geom[i].FastGetSolutionStepValue(FLUID_FRACTION) += N[i] * particle_volume; // no multiplying by element_volume since we devide by it to get the contributed volume fraction
//    }
//}

void CalculateNodalFluidFractionWithLinearWeighing(
    Element::Pointer p_elem,
    const array_1d<double, TDim + 1>& N,
    ParticleType& particle)
{
    const double particle_volume = particle.CalculateVolume();

    for (unsigned int i = 0; i < TDim + 1; i++){
        p_elem->GetGeometry()[i].FastGetSolutionStepValue(FLUID_FRACTION) += N[i] * particle_volume; // no multiplying by element_volume since we devide by it to get the contributed volume fraction
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateNodalFluidFractionByLumpedL2Projection(
        Element::Pointer p_elem,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer p_node)
    {
        // Geometry of the element of the origin model part
        Geometry<Node<3> >& geom = p_elem->GetGeometry();
        array_1d <double, TDim + 1> Ng;
        boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;
        double elemental_volume;
        GeometryUtils::CalculateGeometryData(geom, DN_DX, Ng, elemental_volume);
        const double& radius         = p_node->FastGetSolutionStepValue(RADIUS);
        const double particle_volume = 4.0 * KRATOS_M_PI_3 * mParticlesPerDepthDistance * radius * radius * radius;

        for (unsigned int i = 0; i < TDim + 1; i++){
            geom[i].FastGetSolutionStepValue(FLUID_FRACTION) += (TDim + 1) * N[i] * particle_volume / elemental_volume;
        }
    }

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateFluidFractionGradient(ModelPart& r_model_part)
{
    for (NodeIteratorType i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); i_node++){
        noalias(i_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT)) = ZeroVector(3);
    }

    array_1d <double, 3> grad = ZeroVector(3); // its dimension is always 3
    array_1d <double, TDim + 1> elemental_fluid_fractions;
    array_1d <double, TDim + 1> N; // shape functions vector
    boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
        // computing the shape function derivatives
        Geometry<Node<3> >& geom = ielem->GetGeometry();
        double volume;
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);

        // getting the fluid fraction gradients;

        for (unsigned int i = 0; i < TDim + 1; ++i){
            elemental_fluid_fractions[i] = geom[i].FastGetSolutionStepValue(FLUID_FRACTION);
        }

        array_1d <double, TDim> grad_aux = prod(trans(DN_DX), elemental_fluid_fractions); // its dimension may be 2

        for (unsigned int i = 0; i < TDim; ++i){
            grad[i] = grad_aux[i];
        }

        double nodal_area = volume / (TDim + 1);
        grad *= nodal_area;

        for (unsigned int i = 0; i < TDim + 1; ++i){
            geom[i].FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT) += grad;
        }
    }

    for (NodeIteratorType i_node = r_model_part.NodesBegin(); i_node != r_model_part.NodesEnd(); i_node++){
        i_node->GetSolutionStepValue(FLUID_FRACTION_GRADIENT) /= i_node->GetSolutionStepValue(NODAL_AREA);
    }
 }

//***************************************************************************************************************
//***************************************************************************************************************

void TransferByAveraging(
    const Node<3>::Pointer p_node,
    const ResultNodesContainerType& neighbours,
    const DistanceType& weights,
    Variable<array_1d<double, 3> >& r_destination_variable,
    Variable<array_1d<double, 3> >& r_origin_variable)
{
    if (p_node->IsNot(INSIDE)){

        return;
    }

    const array_1d<double, 3>& origin_data = p_node->FastGetSolutionStepValue(r_origin_variable);

    if (r_origin_variable == HYDRODYNAMIC_FORCE){

        for (unsigned int i = 0; i != neighbours.size(); ++i){
            const double area = neighbours[i]->FastGetSolutionStepValue(NODAL_AREA);
            const double fluid_density = neighbours[i]->FastGetSolutionStepValue(DENSITY);
            const double fluid_fraction =  neighbours[i]->FastGetSolutionStepValue(FLUID_FRACTION);
            array_1d<double, 3> contribution;
            noalias(contribution) = - weights[i] * origin_data / (area * fluid_density * fluid_fraction);
            //neighbours[i]->FastGetSolutionStepValue(r_destination_variable) += contribution;
            neighbours[i]->FastGetSolutionStepValue(HYDRODYNAMIC_REACTION) += contribution;
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateNodalFluidFractionByAveraging( // it is actually calculating its complementary here; (1 - this value) is performed later
    const Node<3>::Pointer p_node,
    const ResultNodesContainerType& neighbours,
    const DistanceType& weights)
{
    unsigned int vector_size = neighbours.size();
    if (vector_size && p_node->Is(INSIDE)){
        const double& radius = p_node->FastGetSolutionStepValue(RADIUS);
        double solid_volume = 4.0 * KRATOS_M_PI_3 * radius * radius * radius;

        for (unsigned int i = 0; i != vector_size; ++i){
            neighbours[i]->GetSolutionStepValue(FLUID_FRACTION) += weights[i] * solid_volume;
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateNodalSolidFractionByAveraging(
    const Node<3>::Pointer p_node,
    const ResultNodesContainerType& neighbours,
    const DistanceType& weights,
    const double averaging_volume_inv)
{
    unsigned int vector_size = neighbours.size();
    if (vector_size && p_node->Is(INSIDE)){
        const double& radius = p_node->FastGetSolutionStepValue(RADIUS);
        double solid_volume = 4.0 * KRATOS_M_PI_3 * radius * radius * radius;

        for (unsigned int i = 0; i != vector_size; ++i){
            neighbours[i]->GetSolutionStepValue(FLUID_FRACTION) += averaging_volume_inv * weights[i] * solid_volume;
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void MultiplyNodalVariableBy(ModelPart& r_model_part, Variable<double>& r_variable, const double& factor)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); i++){
        NodeIteratorType i_node = r_model_part.NodesBegin() + i;
        Node <3> ::Pointer p_node = *(i_node.base());
        p_node->FastGetSolutionStepValue(r_variable) *= factor;
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void MultiplyNodalVariableBy(ModelPart& r_model_part, Variable<array_1d<double, 3> >& r_variable, const double& factor)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); i++){
        NodeIteratorType i_node = r_model_part.NodesBegin() + i;
        Node<3>::Pointer p_node = *(i_node.base());
        p_node->FastGetSolutionStepValue(r_variable) *= factor;
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void ResetDEMVariables(ModelPart& r_dem_model_part)
{
    for (NodeIteratorType node_it = r_dem_model_part.NodesBegin(); node_it != r_dem_model_part.NodesEnd(); ++node_it){

        for (ListIndexType i = 0; i != mDEMCouplingVariables.size(); ++i){
            ClearVariable(node_it, mDEMCouplingVariables[i]);
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void ResetFluidVariables(ModelPart& r_fluid_model_part)
{
    for (NodeIteratorType node_it = r_fluid_model_part.NodesBegin(); node_it != r_fluid_model_part.NodesEnd(); ++node_it){

        ClearVariable(node_it, FLUID_FRACTION);

        array_1d<double, 3>& body_force                = node_it->FastGetSolutionStepValue(BODY_FORCE);
        array_1d<double, 3>& old_hydrodynamic_reaction = node_it->FastGetSolutionStepValue(HYDRODYNAMIC_REACTION);
        //body_force -= old_hydrodynamic_reaction;
        body_force[0] = 0.0;
        body_force[1] = 0.0;
        body_force[2] = -9.81;

        noalias(old_hydrodynamic_reaction) = ZeroVector(3);
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

inline void ClearVariable(const NodeIteratorType& node_it, const VariableData *var)
{
    var->AssignZero(node_it->SolutionStepData().Data(*var));
}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculateFluidNodesMaxNodalArea(ModelPart& r_fluid_model_part)
{
    double max_nodal_area = 0.0;

    for (int i = 0; i < (int)r_fluid_model_part.Nodes().size(); i++){
        NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;
        Node<3>::Pointer p_node = *(i_node.base());
        double nodal_area = p_node->FastGetSolutionStepValue(NODAL_AREA);
        max_nodal_area = std::max(nodal_area, max_nodal_area);
    }

    mMaxNodalAreaInv = 1.0 / max_nodal_area;

    mMustCalculateMaxNodalArea = false;
}

//***************************************************************************************************************
//***************************************************************************************************************

inline void ClearVariable(const NodeIteratorType& node_it, const VariableData& var)
{
    var.AssignZero(node_it->SolutionStepData().Data(var));
}

//***************************************************************************************************************
//***************************************************************************************************************

inline unsigned int GetNearestNode(const array_1d<double, TDim + 1>& N)
{
    double max = N[0];
    unsigned int i_nearest_node = 0;

    for (unsigned int i = 1; i < TDim + 1; ++i){

        if (N[i] > max){
            max = N[i];
            i_nearest_node = i;
        }
    }

    return i_nearest_node;
}

//***************************************************************************************************************
//***************************************************************************************************************

void FillVectorOfSwimmingSpheres(ModelPart& r_dem_model_part){

    mSwimmingSphereElementPointers.resize(r_dem_model_part.Elements().size());
    unsigned int i = 0;

    for (ElementsArrayType::iterator i_elem = r_dem_model_part.ElementsBegin(); i_elem != r_dem_model_part.ElementsEnd(); ++i_elem){
        mSwimmingSphereElementPointers[i] = &(dynamic_cast<Kratos::SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>&>(*i_elem));
        ++i;
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

double inline CalculateDistance(Node<3>::Pointer a, SphericSwimmingParticle<TBaseTypeOfSwimmingParticle>* b){
    array_1d<double, 3> coor_a = a->Coordinates();   
    array_1d<double, 3> coor_b = b->GetGeometry()[0].Coordinates();
    return sqrt((coor_a[0] - coor_b[0]) * (coor_a[0] - coor_b[0]) + (coor_a[1] - coor_b[1]) * (coor_a[1] - coor_b[1]) + (coor_a[2] - coor_b[2]) * (coor_a[2] - coor_b[2]));
}

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

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

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
