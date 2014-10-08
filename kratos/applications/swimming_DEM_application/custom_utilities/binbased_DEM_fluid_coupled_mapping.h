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

//Database includes
#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/spatial_search.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"
#include "../../DEM_application/DEM_application.h"

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

//class BinBasedDEMFluidCoupledMapping
template<std::size_t TDim>
class BinBasedDEMFluidCoupledMapping
{
public:
    ///@name Type Definitions
    ///@{
    typedef ModelPart::ElementsContainerType             ElementsArrayType;
    typedef ElementsArrayType::ContainerType             ResultElementsContainerType;
    typedef std::vector<ResultElementsContainerType>     VectorResultElementsContainerType;
    typedef ModelPart::ElementsContainerType::iterator   ElementIteratorType;


    typedef ModelPart::NodesContainerType                NodesArrayType;
    typedef NodesArrayType::ContainerType                ResultNodesContainerType;
    typedef std::vector<ResultNodesContainerType>        VectorResultNodesContainerType;
    typedef ModelPart::NodesContainerType::iterator      NodeIteratorType;


    typedef std::size_t                                  ListIndexType;
    typedef SpatialSearch::DistanceType                  DistanceType;
    typedef SpatialSearch::VectorDistanceType            VectorDistanceType;

    /// Pointer definition of BinBasedDEMFluidCoupledMapping
    KRATOS_CLASS_POINTER_DEFINITION(BinBasedDEMFluidCoupledMapping<TDim>);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    //----------------------------------------------------------------
    //                       key for coupling_type
    //----------------------------------------------------------------
    //        Averaged variables       |  Fluid Fraction
    //   Fluid-to-DEM | DEM-to-fluid   |
    //----------------------------------------------------------------
    // 0:   Linear         Constant            Constant
    // 1:   Linear         Linear              Constant
    // 2:   Linear         Linear              Linear
    //----------------------------------------------------------------

    BinBasedDEMFluidCoupledMapping(double min_fluid_fraction, const int coupling_type, typename SpatialSearch::Pointer pSpSearch, const int n_particles_per_depth_distance = 1):
                                   mMinFluidFraction(min_fluid_fraction),
                                   mCouplingType(coupling_type),
                                   mParticlesPerDepthDistance(n_particles_per_depth_distance),
                                   mpSpSearch(pSpSearch)
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

    void AddDEMCouplingVariable(const VariableData& rVariable){
        mDEMCouplingVariables.Add(rVariable);
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void AddFluidCouplingVariable(const VariableData& rVariable){
        mFluidCouplingVariables.Add(rVariable);
    }

    //***************************************************************************************************************
    //***************************************************************************************************************
    /// Interpolate fluid data onto the DEM model part.
    /**
      * @param rfluid_model_part: the origin model part from which to project
      * @param rdem_model_part: the destination model part of which we want to interpolate its nodal values
      * @param bin_of_objects_fluid: pre-assembled bin of objects (elelments of the fluid mesh). It is to be constructed separately
      * @see binbased_nodes_in_element_locator
    */
    // data_to_project to DEM mesh = alpha * new_data + (1 - alpha) * old_data

    void InterpolateFromFluidMesh(
        ModelPart& rfluid_model_part,
        ModelPart& rdem_model_part,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
        const double alpha)
    {
        KRATOS_TRY

        // setting interpolated values to 0
        ResetDEMVariables(rdem_model_part);

        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rdem_model_part.Nodes().size();

        #pragma omp parallel for firstprivate(results, N)
        for (int i = 0; i < nparticles; i++){
            NodeIteratorType iparticle = rdem_model_part.NodesBegin() + i;
            Node < 3 > ::Pointer pparticle = *(iparticle.base());

            if (pparticle->IsNot(BLOCKED)){
                Element::Pointer pelement;

                // looking for the fluid element in which the DEM node falls
                bool is_found = bin_of_objects_fluid.FindPointOnMesh(pparticle->Coordinates(), N, pelement, results.begin(), max_results);

                // interpolating the variables

                if (is_found){

                    for (unsigned int j = 0; j != mDEMCouplingVariables.size(); ++j){
                        Project(pelement, N, pparticle, mDEMCouplingVariables[j], alpha);
                      }

                  }

              }

          }

          if (IsDEMVariable(FLUID_ACCEL_PROJECTED)){

              const double delta_time_inv = 1.0 / rfluid_model_part.GetProcessInfo().GetValue(DELTA_TIME);
              MultiplyNodalVariableBy(rdem_model_part, FLUID_ACCEL_PROJECTED, delta_time_inv);

          }

        KRATOS_CATCH("")
    }

    //***************************************************************************************************************
    //***************************************************************************************************************
    // data_to_project to DEM mesh = current fluid data

    void InterpolateFromNewestFluidMesh(
        ModelPart& rfluid_model_part,
        ModelPart& rdem_model_part,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid)
    {
        KRATOS_TRY

        // setting interpolated values to 0
        ResetDEMVariables(rdem_model_part);

        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rdem_model_part.Nodes().size();

        #pragma omp parallel for firstprivate(results, N)
        for (int i = 0; i < nparticles; i++){
            NodeIteratorType iparticle = rdem_model_part.NodesBegin() + i;
            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            
            if (pparticle->IsNot(BLOCKED)){
                Element::Pointer pelement;

                // looking for the fluid element in which the DEM node falls
                bool is_found = bin_of_objects_fluid.FindPointOnMesh(pparticle->Coordinates(), N, pelement, results.begin(), max_results);

                // interpolating variables

                if (is_found){

                    for (unsigned int j = 0; j != mDEMCouplingVariables.size(); ++j){
                        Project(pelement, N, pparticle, mDEMCouplingVariables[j]);
                      }

                  }

              }

          }

        if (IsDEMVariable(FLUID_ACCEL_PROJECTED)){

            const double delta_time_inv = 1.0 / rfluid_model_part.GetProcessInfo().GetValue(DELTA_TIME);
            MultiplyNodalVariableBy(rdem_model_part, FLUID_ACCEL_PROJECTED, delta_time_inv);

          }

        KRATOS_CATCH("")
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    /// Interpolate form the DEM  to the fluid mesh
    /**
      * @param rdem_model_part: the origin model part from which to project
      * @param rfluid_model_part: the destination model part of which we want to interpolate its nodal values
      * @param bin_of_objects_fluid: pre-assembled bin of objects (elelments of the fluid mesh). It is to be constructed separately
      * @see binbased_nodes_in_element_locator
    */

    //***************************************************************************************************************
    //***************************************************************************************************************
    //  data_to_project to fluid mesh = current DEM data

    void InterpolateFromDEMMesh(
        ModelPart& rdem_model_part,
        ModelPart& rfluid_model_part,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid) // this is a bin of objects which contains the FLUID model part
    {
        KRATOS_TRY

        // setting interpolated values to 0
        ResetFluidVariables(rfluid_model_part);
        // calculating the fluid fraction
        InterpolateFluidFraction(rdem_model_part, rfluid_model_part, bin_of_objects_fluid);
        // calculating the rest of fluid variables (particle-fluid force etc.). The solid fraction must be known at this point as it may be used in this step
        InterpolateOtherFluidVariables(rdem_model_part, rfluid_model_part, bin_of_objects_fluid);

        KRATOS_CATCH("")
    }

    //***************************************************************************************************************
    //***************************************************************************************************************
    //  data_to_project to fluid mesh = current DEM data

    void HomogenizeFromDEMMesh(
        ModelPart& rdem_model_part,
        ModelPart& rfluid_model_part,
        const double& search_radius,
        const double& shape_factor) // it is the density function's maximum over its support's radius
    {
        KRATOS_TRY

        // setting interpolated values to 0
        ResetFluidVariables(rfluid_model_part);

        // searching neighbours
        std::vector<double> search_radii; // list of nodal search radii (filter radii). It is a vector since spatial search is designed for varying radius
        VectorResultNodesContainerType vectors_of_neighbouring_balls; // list of nodal arrays of pointers to he node's neighbours
        VectorDistanceType vectors_of_distances; // list of nodal arrays of distances to the node's neighbours

        SearchNodalNeighbours(rfluid_model_part, rdem_model_part, search_radius, search_radii, vectors_of_neighbouring_balls, vectors_of_distances);

        DensityFunctionPolynomial<3> weighing_function(search_radius, shape_factor);

        for (unsigned int i = 0; i < rfluid_model_part.Nodes().size(); ++i){
            weighing_function.ComputeWeights(vectors_of_distances[i], vectors_of_distances[i]);
          }

        for (unsigned int i = 0; i < rfluid_model_part.Nodes().size(); ++i){
            NodeIteratorType inode = rfluid_model_part.NodesBegin() + i;            

            CalculateNodalFluidFractionByAveraging(*(inode.base()), vectors_of_neighbouring_balls[i], vectors_of_distances[i]);

            for (unsigned int j = 0; j != mFluidCouplingVariables.size(); ++j){
                ComputeHomogenizedNodalVariable(*(inode.base()), vectors_of_neighbouring_balls[i], vectors_of_distances[i], mFluidCouplingVariables[j]);
              }

          }

        KRATOS_CATCH("")
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void ComputePostProcessResults(
        ModelPart& rdem_model_part,
        ModelPart& rfluid_model_part,
        ModelPart& rfem_dem_model_part,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
        const ProcessInfo& r_current_process_info)
    {
        const int n_dem_elements = rdem_model_part.Elements().size();
        const int n_fluid_nodes = rfluid_model_part.Nodes().size();

        if (IsDEMVariable(REYNOLDS_NUMBER)){

            #pragma omp parallel for
            for (int i = 0; i < n_dem_elements; i++){
                ElementIteratorType ielem = rdem_model_part.ElementsBegin() + i;
                Geometry< Node<3> >& geom = ielem->GetGeometry();
                double& reynolds_number = geom[0].FastGetSolutionStepValue(REYNOLDS_NUMBER);
                ielem->Calculate(REYNOLDS_NUMBER, reynolds_number, r_current_process_info);
              }
          }
       
        if (IsFluidVariable(MESH_VELOCITY1)){

            #pragma omp parallel for
            for (int i = 0; i < n_fluid_nodes; i++){
                NodeIteratorType inode = rfluid_model_part.NodesBegin() + i;
                double fluid_fraction                         = inode->FastGetSolutionStepValue(FLUID_FRACTION);
                const array_1d<double, 3>& darcy_vel          = inode->FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3>& space_averaged_fluid_vel = inode->FastGetSolutionStepValue(MESH_VELOCITY1);
                space_averaged_fluid_vel                      = darcy_vel / fluid_fraction;
              }

          }

        if (IsFluidVariable(SOLID_FRACTION)){

            #pragma omp parallel for
            for (int i = 0; i < n_fluid_nodes; i++){
                NodeIteratorType inode = rfluid_model_part.NodesBegin() + i;
                double& solid_fraction = inode->FastGetSolutionStepValue(SOLID_FRACTION);
                solid_fraction = 1.0 - inode->FastGetSolutionStepValue(FLUID_FRACTION);
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

    ///@name Protected static Member rVariables
    ///@{

    ///@}
    ///@name Protected member rVariables
    ///@{ template<class T, std::size_t dim>

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{
   vector<unsigned int> mElementsPartition;
   vector<unsigned int> mNodesPartition;
    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member rVariables
    ///@{

    ///@}
    ///@name Member rVariables
    ///@{

    double mMinFluidFraction;
    int mCouplingType;
    int mParticlesPerDepthDistance;
    VariablesList mDEMCouplingVariables;
    VariablesList mFluidCouplingVariables;
    typename SpatialSearch::Pointer mpSpSearch; // it is not be used for some interpolation options

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateFluidFraction(
        ModelPart& rdem_model_part,
        ModelPart& rfluid_model_part,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid) // this is a bin of objects which contains the FLUID model part
    {
        KRATOS_TRY

        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rdem_model_part.Nodes().size();

        #pragma omp parallel for firstprivate(results, N)
        for (int i = 0; i < nparticles; i++){
            NodeIteratorType iparticle = rdem_model_part.NodesBegin() + i;
            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            Element::Pointer pelement;

            // looking for the fluid element in which the DEM node falls
            bool is_found = bin_of_objects_fluid.FindPointOnMesh(pparticle->Coordinates(), N, pelement, results.begin(), max_results);

            // interpolating variables

            if (is_found) {
                DistributeDimensionalContributionToFluidFraction(pelement, N, pparticle);
              }

          }

        CalculateFluidFraction(rfluid_model_part);

        if (IsFluidVariable(FLUID_FRACTION_GRADIENT)){
            CalculateFluidFractionGradient(rfluid_model_part);
          }

        KRATOS_CATCH("")
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateOtherFluidVariables(
        ModelPart& rdem_model_part,
        ModelPart& rfluid_model_part,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid) // this is a bin of objects which contains the FLUID model part
    {
        KRATOS_TRY

        // resetting the variables to be mapped
        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rdem_model_part.Nodes().size();

        #pragma omp parallel for firstprivate(results, N)
        for (int i = 0; i < nparticles; i++){
            NodeIteratorType iparticle = rdem_model_part.NodesBegin() + i;
            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            Element::Pointer pelement;

            // looking for the fluid element in which the DEM node falls
            bool is_found = bin_of_objects_fluid.FindPointOnMesh(pparticle->Coordinates(), N, pelement, results.begin(), max_results);

            // interpolating variables

            if (is_found) {

                for (unsigned int j = 0; j != mFluidCouplingVariables.size(); ++j){
                    Distribute(pelement, N, pparticle, mFluidCouplingVariables[j]);
                  }

              }

          }

        KRATOS_CATCH("")
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void SearchNodalNeighbours(ModelPart& rfluid_model_part,
                               ModelPart& rdem_model_part,
                               const double& search_radius,
                               std::vector<double>& search_radii,
                               VectorResultNodesContainerType& vectors_of_neighbouring_balls,
                               VectorDistanceType& vectors_of_distances)
    {
      KRATOS_TRY

      search_radii.resize(rfluid_model_part.Nodes().size());
      vectors_of_neighbouring_balls.resize(rfluid_model_part.Nodes().size());
      vectors_of_distances.resize(rfluid_model_part.Nodes().size());

      for (unsigned int i = 0; i != rfluid_model_part.Nodes().size(); ++i){
          search_radii[i] = search_radius; // spatial search is designed for varying radius
        }

      mpSpSearch->SearchNodesInRadiusExclusive(rdem_model_part, rfluid_model_part.NodesArray(), search_radii, vectors_of_neighbouring_balls, vectors_of_distances);

      KRATOS_CATCH("")
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

    array_1d<double, 3> CalculateAcceleration(const Geometry< Node < 3 > >& geom, const array_1d<double, TDim + 1 >& N)
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

    double CalculateNormOfSymmetricGradient(const Geometry< Node < 3 > >& geom, const int index)
    {
        Geometry< Node < 3 > >::ShapeFunctionsGradientsType DN_DX;

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

    array_1d<double, 3> CalculateVorticity(const Geometry< Node < 3 > >& geom, const int index)
    {
        Geometry< Node < 3 > >::ShapeFunctionsGradientsType DN_DX;

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

    void Project(Element::Pointer el_it,
                 const array_1d<double, TDim + 1>& N,
                 Node<3>::Pointer pnode,
                 const VariableData *r_destination_variable)
    {

        if (*r_destination_variable == FLUID_DENSITY_PROJECTED){
            Interpolate(el_it, N, pnode, DENSITY, FLUID_DENSITY_PROJECTED);
          }

        else if (*r_destination_variable == FLUID_FRACTION_PROJECTED && IsFluidVariable(FLUID_FRACTION)){
            Interpolate(el_it, N, pnode, FLUID_FRACTION, FLUID_FRACTION_PROJECTED);
          }

        else if (*r_destination_variable == FLUID_VEL_PROJECTED){
            Interpolate(el_it, N, pnode, VELOCITY, FLUID_VEL_PROJECTED);
          }

        else if (*r_destination_variable == PRESSURE_GRAD_PROJECTED){
            Interpolate(el_it, N, pnode, PRESSURE_GRADIENT, PRESSURE_GRAD_PROJECTED);
          }

        else if (*r_destination_variable == FLUID_VISCOSITY_PROJECTED){
            Interpolate(el_it, N, pnode, VISCOSITY, FLUID_VISCOSITY_PROJECTED);
          }

        else if (*r_destination_variable == POWER_LAW_N){
            Interpolate(el_it, N, pnode, POWER_LAW_N, POWER_LAW_N);
          }

        else if (*r_destination_variable == POWER_LAW_K){
            Interpolate(el_it, N, pnode, POWER_LAW_K, POWER_LAW_K);
          }

        else if (*r_destination_variable == GEL_STRENGTH){
            Interpolate(el_it, N, pnode, GEL_STRENGTH, GEL_STRENGTH);
          }

        else if (*r_destination_variable == DISTANCE){
            Interpolate(el_it, N, pnode, DISTANCE, DISTANCE);
          }

        else if (*r_destination_variable == FLUID_ACCEL_PROJECTED){
            InterpolateAcceleration(el_it, N, pnode, FLUID_ACCEL_PROJECTED);
          }

        else if (*r_destination_variable == SHEAR_RATE_PROJECTED){
            InterpolateShearRate(el_it, N, pnode, SHEAR_RATE_PROJECTED);
          }

        else if (*r_destination_variable == FLUID_VORTICITY_PROJECTED){
            InterpolateVorticity(el_it, N, pnode, FLUID_VORTICITY_PROJECTED);
          }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void Project(Element::Pointer el_it,
                 const array_1d<double, TDim + 1> N,
                 Node<3>::Pointer pnode,
                 const VariableData *r_destination_variable,
                 double alpha)
    {

        if (*r_destination_variable == FLUID_DENSITY_PROJECTED){
            Interpolate(el_it, N, pnode, DENSITY, FLUID_DENSITY_PROJECTED, alpha);
          }

        else if (*r_destination_variable == FLUID_FRACTION_PROJECTED && IsFluidVariable(FLUID_FRACTION)){
            Interpolate(el_it, N, pnode, FLUID_FRACTION, FLUID_FRACTION_PROJECTED, alpha);
          }

        else if (*r_destination_variable == FLUID_VEL_PROJECTED){
            Interpolate(el_it, N, pnode, VELOCITY, FLUID_VEL_PROJECTED, alpha);
          }

        else if (*r_destination_variable == PRESSURE_GRAD_PROJECTED){
            Interpolate(el_it, N, pnode, PRESSURE_GRADIENT, PRESSURE_GRAD_PROJECTED, alpha);
          }

        else if (*r_destination_variable == FLUID_VISCOSITY_PROJECTED){
            Interpolate(el_it, N, pnode, VISCOSITY, FLUID_VISCOSITY_PROJECTED, alpha);
          }

        else if (*r_destination_variable == POWER_LAW_N){
            Interpolate(el_it, N, pnode, POWER_LAW_N, POWER_LAW_N, alpha);
          }

        else if (*r_destination_variable == POWER_LAW_K){
            Interpolate(el_it, N, pnode, POWER_LAW_K, POWER_LAW_K, alpha);
          }

        else if (*r_destination_variable == GEL_STRENGTH){
            Interpolate(el_it, N, pnode, GEL_STRENGTH, GEL_STRENGTH, alpha);
          }

        else if (*r_destination_variable == DISTANCE){
            Interpolate(el_it, N, pnode, DISTANCE, DISTANCE, alpha);
          }

        else if (*r_destination_variable == FLUID_ACCEL_PROJECTED){
            InterpolateAcceleration(el_it, N, pnode, FLUID_ACCEL_PROJECTED);
          }

        else if (*r_destination_variable == SHEAR_RATE_PROJECTED){
            InterpolateShearRate(el_it, N, pnode, SHEAR_RATE_PROJECTED, alpha);
          }

        else if (*r_destination_variable == FLUID_VORTICITY_PROJECTED){
            InterpolateVorticity(el_it, N, pnode, FLUID_VORTICITY_PROJECTED, alpha);
          }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void DistributeDimensionalContributionToFluidFraction(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode)
    {

        if (mCouplingType == 0 || mCouplingType == 1){
            CalculateNodalFluidFractionWithConstantWeighing(el_it, N, pnode);
          }

        else if (mCouplingType == 2){
            CalculateNodalFluidFractionWithLinearWeighing(el_it, N, pnode);
          }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void Distribute(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        const VariableData *r_destination_variable)
    {

        if (mCouplingType == 0){

            if (*r_destination_variable == BODY_FORCE){
                TransferWithConstantWeighing(el_it, N, pnode, BODY_FORCE, HYDRODYNAMIC_FORCE);
              }

          }

        else if (mCouplingType == 1){

            if (*r_destination_variable == BODY_FORCE){
                TransferWithLinearWeighing(el_it, N, pnode, BODY_FORCE, HYDRODYNAMIC_FORCE);
              }

          }

        else if (mCouplingType == 2){
            
            if (*r_destination_variable == BODY_FORCE){
                TransferWithLinearWeighing(el_it, N, pnode, BODY_FORCE, HYDRODYNAMIC_FORCE);
              }

          }

        else if (mCouplingType == - 1){

             if (*r_destination_variable == BODY_FORCE){
                 TransferWithLinearWeighing(el_it, N, pnode, BODY_FORCE, HYDRODYNAMIC_FORCE);
               }

          }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void ComputeHomogenizedNodalVariable(
        const Node<3>::Pointer pnode,
        const ResultNodesContainerType& neighbours,
        const DistanceType& weights,
        const VariableData *r_destination_variable)
    {

        if (mCouplingType < 0 || mCouplingType > - 1){

            if (*r_destination_variable == BODY_FORCE){
                TransferByAveraging(pnode, neighbours, weights, BODY_FORCE, HYDRODYNAMIC_FORCE);
              }

          }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void  CalculateFluidFraction(ModelPart& rfluid_model_part)
    {

        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), rfluid_model_part.Nodes().size(), mNodesPartition);

        #pragma omp parallel for
        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

            for (NodesArrayType::iterator inode = this->GetNodePartitionBegin(rfluid_model_part, k); inode != this->GetNodePartitionEnd(rfluid_model_part, k); ++inode){
                double& fluid_fraction          = inode->FastGetSolutionStepValue(FLUID_FRACTION);
                fluid_fraction                  = 1.0 - fluid_fraction / inode->FastGetSolutionStepValue(NODAL_AREA);

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
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_origin_variable,
        Variable<array_1d<double, 3> >& r_destination_variable)
    {

        // Geometry of the element of the origin model part
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        // Destination data
        array_1d<double, 3>& step_data = pnode->FastGetSolutionStepValue(r_destination_variable);

        step_data = N[0] * geom[0].FastGetSolutionStepValue(r_origin_variable);

        for (unsigned int i = 1; i < TDim + 1; i++){
            step_data += N[i] * geom[i].FastGetSolutionStepValue(r_origin_variable);
          }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_origin_variable,
        Variable<array_1d<double, 3> >& r_destination_variable,
        double alpha)
    {

        // Geometry of the element of the origin model part
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        // Destination data
        array_1d<double, 3>& step_data = pnode->FastGetSolutionStepValue(r_destination_variable);
        step_data += N[0] * (alpha * geom[0].FastGetSolutionStepValue(r_origin_variable) + (1 - alpha) * geom[0].FastGetSolutionStepValue(r_origin_variable, 1));

        for (unsigned int i = 1; i < TDim + 1; i++){
            step_data += N[i] * (alpha * geom[i].FastGetSolutionStepValue(r_origin_variable) + (1 - alpha) * geom[i].FastGetSolutionStepValue(r_origin_variable, 1));
          }

     }

    //***************************************************************************************************************
    //***************************************************************************************************************
    // projecting a scalar

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<double>& r_origin_variable,
        Variable<double>& r_destination_variable)
    {

        // Geometry of the element of the origin model part
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        // Destination data
        double& step_data = pnode->FastGetSolutionStepValue(r_destination_variable);

        step_data = N[0] * geom[0].FastGetSolutionStepValue(r_origin_variable);

        for (unsigned int i = 1; i < TDim + 1; i++){
            step_data += N[i] * geom[i].FastGetSolutionStepValue(r_origin_variable);
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<double>& r_origin_variable,
        Variable<double>& r_destination_variable,
        double alpha)
    {

      // Geometry of the element of the origin model part
      Geometry< Node<3> >& geom = el_it->GetGeometry();

      // Destination data
      double& step_data = pnode->FastGetSolutionStepValue(r_destination_variable);
      step_data += N[0] * (alpha * geom[0].FastGetSolutionStepValue(r_origin_variable) + (1 - alpha) * geom[0].FastGetSolutionStepValue(r_origin_variable, 1));

      for (unsigned int i = 1; i < TDim + 1; i++){
          step_data += N[i] * (alpha * geom[i].FastGetSolutionStepValue(r_origin_variable) + (1 - alpha) * geom[i].FastGetSolutionStepValue(r_origin_variable, 1));
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateAcceleration(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_destination_variable)
    {
        pnode->FastGetSolutionStepValue(r_destination_variable) = CalculateAcceleration(el_it->GetGeometry(), N);
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateShearRate(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<double>& r_destination_variable)
    {
        double shear_rate = CalculateNormOfSymmetricGradient(el_it->GetGeometry(), 0);
        double& step_data = pnode->FastGetSolutionStepValue(r_destination_variable);

        step_data = shear_rate;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateShearRate(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<double>& r_destination_variable,
        double alpha)
    {
        double shear_rate       = CalculateNormOfSymmetricGradient(el_it->GetGeometry(), 0);
        double prev_shear_rate  = CalculateNormOfSymmetricGradient(el_it->GetGeometry(), 1);
        double& step_data       = pnode->FastGetSolutionStepValue(r_destination_variable);

        step_data = alpha * shear_rate + (1.0 - alpha) * prev_shear_rate;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateVorticity(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_destination_variable)
    {
        array_1d<double, 3> vorticity  = CalculateVorticity(el_it->GetGeometry(), 0);
        array_1d<double, 3>& step_data = pnode->FastGetSolutionStepValue(r_destination_variable);

        step_data = vorticity;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateVorticity(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_destination_variable,
        double alpha)
    {
        array_1d<double, 3> vorticity      = CalculateVorticity(el_it->GetGeometry(), 0);
        array_1d<double, 3> prev_vorticity = CalculateVorticity(el_it->GetGeometry(), 1);
        array_1d<double, 3>& step_data     = pnode->FastGetSolutionStepValue(r_destination_variable);

        step_data = alpha * vorticity + (1.0 - alpha) * prev_vorticity;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void TransferWithConstantWeighing(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_destination_variable,
        Variable<array_1d<double, 3> >& r_origin_variable)
    {
        // Geometry of the element of the destination model part
        Geometry< Node<3> >& geom              = el_it->GetGeometry();
        unsigned int i_nearest_node            = GetNearestNode(N);
        const array_1d<double, 3>& origin_data = pnode->FastGetSolutionStepValue(r_origin_variable);
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
        Element::Pointer el_it,
        const array_1d<double,TDim + 1>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_destination_variable,
        Variable<array_1d<double, 3> >& r_origin_variable)
    {
        // Geometry of the element of the destination model part
        Geometry< Node<3> >& geom = el_it->GetGeometry();
        const array_1d<double, 3>& origin_data = pnode->FastGetSolutionStepValue(r_origin_variable);

        if (r_origin_variable == HYDRODYNAMIC_FORCE){
            array_1d<double, 3>& hydrodynamic_reaction = geom[0].FastGetSolutionStepValue(HYDRODYNAMIC_REACTION);
            const double& fluid_fraction               = geom[0].FastGetSolutionStepValue(FLUID_FRACTION);
            const double& nodal_volume                 = geom[0].FastGetSolutionStepValue(NODAL_AREA);
            const double& density                      = geom[0].FastGetSolutionStepValue(DENSITY);
            hydrodynamic_reaction = - (mParticlesPerDepthDistance * N[0] / (fluid_fraction * density * nodal_volume)) * origin_data;

            for (unsigned int i = 1; i < TDim + 1; i++){
                array_1d<double, 3>& hydrodynamic_reaction = geom[i].FastGetSolutionStepValue(HYDRODYNAMIC_REACTION);
                const double& fluid_fraction               = geom[i].FastGetSolutionStepValue(FLUID_FRACTION);
                const double& nodal_volume                 = geom[i].FastGetSolutionStepValue(NODAL_AREA);
                const double& density                      = geom[i].FastGetSolutionStepValue(DENSITY);
                hydrodynamic_reaction -= (mParticlesPerDepthDistance * N[i] / (fluid_fraction * density * nodal_volume)) * origin_data;
              }

          }

        else {
            std::cout << "Variable " << r_origin_variable << " is not supported for transference with linear weights" ;
          }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void CalculateNodalFluidFractionWithConstantWeighing(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer pnode)
    {
        unsigned int i_nearest_node = GetNearestNode(N);

        // Geometry of the element of the destination model part
        const double& radius         = pnode->FastGetSolutionStepValue(RADIUS);
        const double particle_volume = 1.33333333333333333333 * KRATOS_M_PI * mParticlesPerDepthDistance * radius * radius * radius;

        (el_it->GetGeometry())[i_nearest_node].FastGetSolutionStepValue(FLUID_FRACTION) += particle_volume;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void CalculateNodalFluidFractionWithLinearWeighing(
        Element::Pointer el_it,
        const array_1d<double, TDim + 1>& N,
        Node<3>::Pointer inode)
    {
        // Geometry of the element of the origin model part
        Geometry< Node<3> >& geom = el_it->GetGeometry();
        array_1d<double,TDim + 1> N_2; // a dummy since we are not interested in its value at the Gauss points
        boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX; // a dummy
        double element_volume;
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N_2, element_volume);

        // Destination data
        const double& radius         = (inode)->FastGetSolutionStepValue(RADIUS);
        const double particle_volume = 1.33333333333333333333 * KRATOS_M_PI * mParticlesPerDepthDistance * radius * radius * radius;

        for (unsigned int i = 0; i < TDim + 1; i++){
            (el_it->GetGeometry())[i].FastGetSolutionStepValue(FLUID_FRACTION) += N[i] * element_volume * particle_volume;
          }

    }
    //***************************************************************************************************************
    //***************************************************************************************************************

    void CalculateFluidFractionGradient(ModelPart& r_model_part)
    {

        for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
            noalias(inode->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT)) = ZeroVector(3);
        }

        array_1d <double, TDim + 1 > elemental_fluid_fractions;
        array_1d <double, TDim> grad;
        array_1d <double, TDim + 1 > N; // shape functions vector
        boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;

        for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
            // computing the shape function derivatives
            Geometry< Node < 3 > >& geom = ielem->GetGeometry();
            double volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);

            // getting the fluid fraction gradients;

            for (unsigned int i = 0; i < TDim + 1; ++i){
                elemental_fluid_fractions[i] = geom[i].FastGetSolutionStepValue(FLUID_FRACTION);
              }

            noalias(grad) = prod(trans(DN_DX), elemental_fluid_fractions);
            double nodal_area = volume / (TDim + 1);
            grad *= nodal_area;

            for (unsigned int i = 0; i < TDim + 1; ++i){
                geom[i].FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT) += grad;
              }

        }

        for (NodeIteratorType inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
            inode->GetSolutionStepValue(FLUID_FRACTION_GRADIENT) /= inode->GetSolutionStepValue(NODAL_AREA);
          }

     }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void TransferByAveraging(
        const Node<3>::Pointer pnode,
        const ResultNodesContainerType& neighbours,
        const DistanceType& weights,
        Variable<array_1d<double, 3> >& r_destination_variable,
        Variable<array_1d<double, 3> >& r_origin_variable)
    {
        array_1d<double, 3>& destination_data = pnode->FastGetSolutionStepValue(r_destination_variable);
        const double& fluid_density  = pnode->FastGetSolutionStepValue(DENSITY);
        const double& fluid_fraction = pnode->FastGetSolutionStepValue(FLUID_FRACTION);
        array_1d<double, 3> neighbours_contribution = ZeroVector(3);

        for (unsigned int i = 0; i < neighbours.size(); ++i){
            const array_1d<double, 3>& origin_data = neighbours[i]->FastGetSolutionStepValue(r_origin_variable);
            neighbours_contribution += origin_data * weights[i];
          }

        if (r_origin_variable == HYDRODYNAMIC_FORCE){
            neighbours_contribution /= (fluid_density * fluid_fraction);
          }

        destination_data += neighbours_contribution;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void CalculateNodalFluidFractionByAveraging(
        const Node<3>::Pointer pnode,
        const ResultNodesContainerType& neighbours,
        const DistanceType& weights)
    {
        double& destination_data = pnode->FastGetSolutionStepValue(FLUID_FRACTION);

        for (unsigned int i = 0; i != neighbours.size(); ++i){
            const double& radius = neighbours[i]->FastGetSolutionStepValue(RADIUS);
            destination_data += weights[i] * radius * radius * radius;
          }

        destination_data = 1 -  4 / 3 * KRATOS_M_PI * destination_data;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void MultiplyNodalVariableBy(ModelPart& r_model_part, Variable<double>& r_variable, const double& factor){

        #pragma omp parallel for
        for (unsigned int i = 0; i < r_model_part.Nodes().size(); i++){
            NodeIteratorType inode = r_model_part.NodesBegin() + i;
            Node < 3 > ::Pointer pnode = *(inode.base());
            pnode->FastGetSolutionStepValue(r_variable) *= factor;
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void MultiplyNodalVariableBy(ModelPart& r_model_part, Variable<array_1d<double, 3> >& r_variable, const double& factor){

        #pragma omp parallel for
        for (unsigned int i = 0; i < r_model_part.Nodes().size(); i++){
            NodeIteratorType inode = r_model_part.NodesBegin() + i;
            Node < 3 > ::Pointer pnode = *(inode.base());
            pnode->FastGetSolutionStepValue(r_variable) *= factor;
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************


    void ResetDEMVariables(ModelPart& rdem_model_part)
    {

        for (NodeIteratorType node_it = rdem_model_part.NodesBegin(); node_it != rdem_model_part.NodesEnd(); ++node_it){

            for (ListIndexType i = 0; i != mDEMCouplingVariables.size(); ++i){
                ClearVariable(node_it, mDEMCouplingVariables[i]);
              }

          }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void ResetFluidVariables(ModelPart& rfluid_model_part)
    {

        for (NodeIteratorType node_it = rfluid_model_part.NodesBegin(); node_it != rfluid_model_part.NodesEnd(); ++node_it){

            ClearVariable(node_it, FLUID_FRACTION);

            array_1d<double, 3>& body_force                = node_it->FastGetSolutionStepValue(BODY_FORCE);
            array_1d<double, 3>& old_hydrodynamic_reaction = node_it->FastGetSolutionStepValue(HYDRODYNAMIC_REACTION);
            body_force -= old_hydrodynamic_reaction;

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

    inline void ClearVariable(const NodeIteratorType& node_it, const VariableData& var)
    {
        var.AssignZero(node_it->SolutionStepData().Data(var));
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline double CalculateVol(const double x0, const double y0,
                               const double x1, const double y1,
                               const double x2, const double y2)
    {
        return 0.5 * ((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0));
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline double CalculateVol(const double x0, const double y0, const double z0,
                               const double x1, const double y1, const double z1,
                               const double x2, const double y2, const double z2,
                               const double x3, const double y3, const double z3)
    {
        double x10 = x1 - x0;
        double y10 = y1 - y0;
        double z10 = z1 - z0;

        double x20 = x2 - x0;
        double y20 = y2 - y0;
        double z20 = z2 - z0;

        double x30 = x3 - x0;
        double y30 = y3 - y0;
        double z30 = z3 - z0;

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 +
                      y10 * z20 * x30 - y10 * x20 * z30 +
                      z10 * x20 * y30 - z10 * y20 * x30;

        return  detJ * 0.1666666666666666666667;

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline unsigned int GetNearestNode(const array_1d<double, TDim + 1>& N)

    {
      unsigned int i_nearest_node = 0;
      double max                  = N[0];

      for (unsigned int inode = 1; inode < TDim + 1; ++inode){

          if (N[inode] > max) {
              max = N[inode];
              i_nearest_node = inode;
            }

        }

      return i_nearest_node;
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
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BinBasedDEMFluidCoupledMapping<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_BINBASED_DEM_FLUID_COUPLED_MAPPING  defined


