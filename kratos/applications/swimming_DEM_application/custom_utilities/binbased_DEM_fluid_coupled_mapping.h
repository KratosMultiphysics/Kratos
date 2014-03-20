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

//Database includes
#include "spatial_containers/spatial_containers.h"
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
    typedef ModelPart::NodesContainerType                NodesArrayType;
    typedef ModelPart::NodesContainerType::ContainerType NodesContainerType;
    typedef ModelPart::NodesContainerType::iterator      NodeIteratorType;
    typedef ModelPart::ElementsContainerType::iterator   ElementIteratorType;
    typedef std::size_t                                  ListIndexType;


    /// Pointer definition of BinBasedDEMFluidCoupledMapping
    KRATOS_CLASS_POINTER_DEFINITION(BinBasedDEMFluidCoupledMapping<TDim>);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    //----------------------------------------------------------------
    //                       key for coupling_type
    //----------------------------------------------------------------
    //        Averaged variables       |  Solid Fraction
    //   Fluid-to-DEM | DEM-to-fluid   |
    //----------------------------------------------------------------
    // 0:   Linear         Constant            Constant
    // 1:   Linear         Linear              Constant
    // 2:   Linear         Linear              Linear
    //----------------------------------------------------------------

    BinBasedDEMFluidCoupledMapping(double max_solid_fraction, const int coupling_type, const int n_particles_per_depth_distance = 1):
                                   mMaxSolidFraction(max_solid_fraction),
                                   mCouplingType(coupling_type),
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

    /// Interpolate fluid data onto the DEM model part.
    /**
      * @param rfluid_model_part: the origin model part from which to project
      * @param rdem_model_part: the destination model part of which we want to interpolate its nodal values
      * @param bin_of_objects_fluid: pre-assembled bin of objects (elelments of the fluid mesh). It is to be constructed separately
      * @see binbased_nodes_in_element_locator
    */

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
    // data_to_project to DEM mesh = alpha * new_data + (1 - alpha) * old_data

    void InterpolateFromFluidMesh(
        ModelPart& rfluid_model_part,
        ModelPart& rdem_model_part,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
        const double alpha)
    {
        KRATOS_TRY

        //Clear all the variables to be mapped

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

        // resetting all variables to be interpolated
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

        // resetting the variables to be mapped
        ResetFluidVariables(rfluid_model_part);

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

        if (IsFluidVariable(SOLID_FRACTION)){
            DivideSolidFractionByNodalVolume(rfluid_model_part);
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

        if (IsDEMVariable(REYNOLDS_NUMBER)){

            #pragma omp parallel
            for (int i = 0; i < n_dem_elements; i++){
                ElementIteratorType ielem = rdem_model_part.ElementsBegin() + i;
                Geometry< Node<3> >& geom = ielem->GetGeometry();
                double& reynolds_number = geom[0].FastGetSolutionStepValue(REYNOLDS_NUMBER, 0);
                ielem->Calculate(REYNOLDS_NUMBER, reynolds_number, r_current_process_info);
            }

        }

        const int n_fluid_nodes = rfluid_model_part.Nodes().size();

        if (IsFluidVariable(MESH_VELOCITY1)){

            #pragma omp parallel
            for (int i = 0; i < n_fluid_nodes; i++){
                NodeIteratorType inode = rfluid_model_part.NodesBegin() + i;
                double fluid_fraction                         = 1 - inode->FastGetSolutionStepValue(SOLID_FRACTION, 0);
                const array_1d<double, 3>& darcy_vel          = inode->FastGetSolutionStepValue(VELOCITY, 0);
                array_1d<double, 3>& space_averaged_fluid_vel = inode->FastGetSolutionStepValue(MESH_VELOCITY1, 0);
                space_averaged_fluid_vel                      = darcy_vel / fluid_fraction;
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

    double mMaxSolidFraction;
    int mCouplingType;
    int mParticlesPerDepthDistance;
    VariablesList mDEMCouplingVariables;
    VariablesList mFluidCouplingVariables;

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
        array_1d<double, 3> derivatives;

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
                 const array_1d<double,4>& N,
                 Node<3>::Pointer pnode,
                 const VariableData *r_destination_variable)
    {

        if (*r_destination_variable == FLUID_DENSITY_PROJECTED){
            Interpolate(el_it, N, pnode, DENSITY, FLUID_DENSITY_PROJECTED);
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
                 const array_1d<double,4>& N,
                 Node<3>::Pointer pnode,
                 const VariableData *r_destination_variable,
                 double alpha)
    {

        if (*r_destination_variable == FLUID_DENSITY_PROJECTED){
            Interpolate(el_it, N, pnode, DENSITY, FLUID_DENSITY_PROJECTED, alpha);
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

        else if (*r_destination_variable == SHEAR_RATE_PROJECTED){
            InterpolateShearRate(el_it, N, pnode, SHEAR_RATE_PROJECTED, alpha);
        }

        else if (*r_destination_variable == FLUID_VORTICITY_PROJECTED){
            InterpolateVorticity(el_it, N, pnode, FLUID_VORTICITY_PROJECTED, alpha);
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void Distribute(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        const VariableData *r_destination_variable)
    {

        if (mCouplingType == 0){

            if (*r_destination_variable == BODY_FORCE){
                TransferWithConstantWeighing(el_it, N, pnode, BODY_FORCE, HYDRODYNAMIC_FORCE);
            }

            else if (*r_destination_variable == SOLID_FRACTION){
                CalculateNodalSolidFractionWithConstantWeighing(el_it, N, pnode);
            }

        }

        else if (mCouplingType == 1){

            if (*r_destination_variable == BODY_FORCE){
                TransferWithLinearWeighing(el_it, N, pnode, BODY_FORCE, HYDRODYNAMIC_FORCE);
            }

            else if (*r_destination_variable == SOLID_FRACTION){
                CalculateNodalSolidFractionWithConstantWeighing(el_it, N, pnode);
            }

        }

        else if (mCouplingType == 2){

            if (*r_destination_variable == BODY_FORCE){
                TransferWithLinearWeighing(el_it, N, pnode, BODY_FORCE, HYDRODYNAMIC_FORCE);
            }

            else if (*r_destination_variable == SOLID_FRACTION && IsFluidVariable(SOLID_FRACTION_GRADIENT)){
                CalculateNodalSolidFractionWithLinearWeighing(el_it, N, pnode);
            }

            else if (*r_destination_variable == SOLID_FRACTION){
                CalculateNodalSolidFractionWithLinearWeighingNoGradient(el_it, N, pnode);
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

    void  DivideSolidFractionByNodalVolume(ModelPart& rfluid_model_part)
    {

        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), rfluid_model_part.Nodes().size(), mNodesPartition);

        #pragma omp parallel for
        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

            for (NodesArrayType::iterator inode = this->GetNodePartitionBegin(rfluid_model_part, k); inode != this->GetNodePartitionEnd(rfluid_model_part, k); ++inode){
                double& solid_fraction = inode->FastGetSolutionStepValue(SOLID_FRACTION, 0);
                solid_fraction /= inode->FastGetSolutionStepValue(NODAL_AREA, 0);

                if (solid_fraction > mMaxSolidFraction){
                    solid_fraction = mMaxSolidFraction;
                  }
              }

        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************
    // project an array1D (2Dversion)

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double, 3>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_origin_variable,
        Variable<array_1d<double, 3> >& r_destination_variable)
    {

        // Geometry element of the origin model part
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        // getting the data of the solution step
        array_1d<double, 3>& step_data = (pnode)->FastGetSolutionStepValue(r_destination_variable, 0);
        const array_1d<double, 3>& velocity = (pnode)->FastGetSolutionStepValue(VELOCITY, 0);
        const array_1d<double, 3>& node0_data = geom[0].FastGetSolutionStepValue(r_origin_variable, 0);
        const array_1d<double, 3>& node1_data = geom[1].FastGetSolutionStepValue(r_origin_variable, 0);
        const array_1d<double, 3>& node2_data = geom[2].FastGetSolutionStepValue(r_origin_variable, 0);

        // copying this data in the position of the vector we are interested in

        for (unsigned int j= 0; j< TDim; j++){
            step_data[j] = N[0] * node0_data[j] + N[1] * node1_data[j] + N[2] * node2_data[j];
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************
    // projecting an array1D 3Dversion

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_origin_variable,
        Variable<array_1d<double, 3> >& r_destination_variable)
    {
        // Geometry element of the origin model part
        Geometry< Node < 3 > >& geom = el_it->GetGeometry();

        // getting the data of the solution step
        array_1d<double, 3 > & step_data = (pnode)->FastGetSolutionStepValue(r_destination_variable);

        const array_1d<double, 3 > & node0_data = geom[0].FastGetSolutionStepValue(r_origin_variable);
        const array_1d<double, 3 > & node1_data = geom[1].FastGetSolutionStepValue(r_origin_variable);
        const array_1d<double, 3 > & node2_data = geom[2].FastGetSolutionStepValue(r_origin_variable);
        const array_1d<double, 3 > & node3_data = geom[3].FastGetSolutionStepValue(r_origin_variable);

        for (unsigned int j = 0; j < TDim; j++) {
            step_data[j] = N[0] * node0_data[j] + N[1] * node1_data[j] + N[2] * node2_data[j] + N[3] * node3_data[j];
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_origin_variable,
        Variable<array_1d<double, 3> >& r_destination_variable,
        double alpha)
    {

        // Geometry element of the origin model part
        Geometry< Node < 3 > >& geom = el_it->GetGeometry();

        // getting the data of the solution step
        array_1d<double, 3 > & step_data = (pnode)->FastGetSolutionStepValue(r_destination_variable, 0);

        const array_1d<double, 3 > & node0_data = geom[0].FastGetSolutionStepValue(r_origin_variable, 0);
        const array_1d<double, 3 > & node1_data = geom[1].FastGetSolutionStepValue(r_origin_variable, 0);
        const array_1d<double, 3 > & node2_data = geom[2].FastGetSolutionStepValue(r_origin_variable, 0);
        const array_1d<double, 3 > & node3_data = geom[3].FastGetSolutionStepValue(r_origin_variable, 0);

        const array_1d<double, 3 > & node0_data_prev = geom[0].FastGetSolutionStepValue(r_origin_variable, 1);
        const array_1d<double, 3 > & node1_data_prev = geom[1].FastGetSolutionStepValue(r_origin_variable, 1);
        const array_1d<double, 3 > & node2_data_prev = geom[2].FastGetSolutionStepValue(r_origin_variable, 1);
        const array_1d<double, 3 > & node3_data_prev = geom[3].FastGetSolutionStepValue(r_origin_variable, 1);

        for (unsigned int j = 0; j < TDim; j++) {
            step_data[j] = alpha * (N[0] * node0_data[j]      + N[1] * node1_data[j]      + N[2] * node2_data[j]      + N[3] * node3_data[j]) +
                   (1.0 - alpha) * (N[0] * node0_data_prev[j] + N[1] * node1_data_prev[j] + N[2] * node2_data_prev[j] + N[3] * node3_data_prev[j]);
        }

     }

    //***************************************************************************************************************
    //**************************************************************************************************************
    // projecting a scalar 2Dversion

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double, 3>& N,
        Node<3>::Pointer pnode,
        Variable<double>& r_origin_variable,
        Variable<double>& r_destination_variable)
    {

        Geometry< Node < 3 > >& geom = el_it->GetGeometry();

        // getting the data of the solution step
        double& step_data = (pnode)->FastGetSolutionStepValue(r_destination_variable, 0);
        const double node0_data = geom[0].FastGetSolutionStepValue(r_origin_variable, 0);
        const double node1_data = geom[1].FastGetSolutionStepValue(r_origin_variable, 0);
        const double node2_data = geom[2].FastGetSolutionStepValue(r_origin_variable, 0);

        // copying this data in the position of the vector we are interested in
        step_data = N[0] * node0_data + N[1] * node1_data + N[2] * node2_data;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************
    // projecting a scalar 3Dversion

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<double>& r_origin_variable,
        Variable<double>& r_destination_variable)
    {
        // Geometry element of the origin model part
        Geometry< Node < 3 > >& geom = el_it->GetGeometry();

        // getting the data of the solution step
        double& step_data = (pnode)->FastGetSolutionStepValue(r_destination_variable, 0);
        const double node0_data = geom[0].FastGetSolutionStepValue(r_origin_variable, 0);
        const double node1_data = geom[1].FastGetSolutionStepValue(r_origin_variable, 0);
        const double node2_data = geom[2].FastGetSolutionStepValue(r_origin_variable, 0);
        const double node3_data = geom[3].FastGetSolutionStepValue(r_origin_variable, 0);
        step_data = N[0] * node0_data + N[1] * node1_data + N[2] * node2_data + N[3] * node3_data;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<double>& r_origin_variable,
        Variable<double>& r_destination_variable,
        double alpha)
    {
        // Geometry element of the origin model part
        Geometry< Node < 3 > >& geom = el_it->GetGeometry();

        // getting the data of the solution step
        double& step_data = (pnode)->FastGetSolutionStepValue(r_destination_variable, 0);

        const double node0_data      = geom[0].FastGetSolutionStepValue(r_origin_variable, 0);
        const double node1_data      = geom[1].FastGetSolutionStepValue(r_origin_variable, 0);
        const double node2_data      = geom[2].FastGetSolutionStepValue(r_origin_variable, 0);
        const double node3_data      = geom[3].FastGetSolutionStepValue(r_origin_variable, 0);

        const double node0_data_prev = geom[0].FastGetSolutionStepValue(r_origin_variable, 1);
        const double node1_data_prev = geom[1].FastGetSolutionStepValue(r_origin_variable, 1);
        const double node2_data_prev = geom[2].FastGetSolutionStepValue(r_origin_variable, 1);
        const double node3_data_prev = geom[3].FastGetSolutionStepValue(r_origin_variable, 1);

        step_data = alpha * (N[0] * node0_data + N[1] * node1_data + N[2] * node2_data + N[3] * node3_data) +
            (1.0 - alpha) * (N[0] * node0_data_prev + N[1] * node1_data_prev + N[2] * node2_data_prev + N[3] * node3_data_prev);
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateShearRate(
        Element::Pointer el_it,
        const array_1d<double, 3>& N,
        Node<3>::Pointer pnode,
        Variable<double>& r_destination_variable)
    {
        double shear_rate = CalculateNormOfSymmetricGradient(el_it->GetGeometry(), 0);
        double& step_data = (pnode)->FastGetSolutionStepValue(r_destination_variable, 0);

        step_data = shear_rate;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateShearRate(
        Element::Pointer el_it,
        const array_1d<double, 3>& N,
        Node<3>::Pointer pnode,
        Variable<double>& r_destination_variable,
        double alpha)
    {
        double shear_rate       = CalculateNormOfSymmetricGradient(el_it->GetGeometry(), 0);
        double prev_shear_rate  = CalculateNormOfSymmetricGradient(el_it->GetGeometry(), 1);
        double& step_data       = (pnode)->FastGetSolutionStepValue(r_destination_variable, 0);

        step_data = alpha * shear_rate + (1.0 - alpha) * prev_shear_rate;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateVorticity(
        Element::Pointer el_it,
        const array_1d<double, 3>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_destination_variable)
    {
        array_1d<double, 3> vorticity  = CalculateVorticity(el_it->GetGeometry(), 0);
        array_1d<double, 3>& step_data = (pnode)->FastGetSolutionStepValue(r_destination_variable, 0);

        step_data = vorticity;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void InterpolateVorticity(
        Element::Pointer el_it,
        const array_1d<double, 3>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_destination_variable,
        double alpha)
    {
        array_1d<double, 3> vorticity      = CalculateVorticity(el_it->GetGeometry(), 0);
        array_1d<double, 3> prev_vorticity = CalculateVorticity(el_it->GetGeometry(), 1);
        array_1d<double, 3>& step_data     = (pnode)->FastGetSolutionStepValue(r_destination_variable, 0);

        step_data = alpha * vorticity + (1.0 - alpha) * prev_vorticity;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************
    // 2D version

    void TransferWithConstantWeighing(
        Element::Pointer el_it,
        const array_1d<double, 3>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_destination_variable,
        Variable<array_1d<double, 3> >& r_origin_variable)
    {
        // Geometry element of the origin model part
        Geometry< Node<3> >& geom              = el_it->GetGeometry();
        unsigned int i_nearest_node            = GetNearestNode(N);
        const array_1d<double, 3>& origin_data = (pnode)->FastGetSolutionStepValue(r_origin_variable, 0);
        array_1d<double, 3>& destination_data  = geom[i_nearest_node].FastGetSolutionStepValue(r_destination_variable, 0);

        if (r_origin_variable == HYDRODYNAMIC_FORCE){
            const double density                       = geom[i_nearest_node].FastGetSolutionStepValue(DENSITY, 0);
            const double nodal_volume                  = geom[i_nearest_node].FastGetSolutionStepValue(NODAL_AREA, 0);
            array_1d<double, 3>& old_drag_contribution = geom[i_nearest_node].FastGetSolutionStepValue(HYDRODYNAMIC_REACTION, 0);
            const double nodal_mass_inv                = mParticlesPerDepthDistance / (density * nodal_volume);

            for (unsigned int j= 0; j< TDim; j++){
                array_1d<double, 3> origin_nodal_contribution = - nodal_mass_inv * origin_data;
                destination_data += origin_nodal_contribution;
                old_drag_contribution += origin_nodal_contribution;
            }

        }

        else {
            std::cout << "Variable " << r_origin_variable << " is not supported for transference with constant weights";
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void TransferWithLinearWeighing(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double, 3> >& r_destination_variable,
        Variable<array_1d<double, 3> >& r_origin_variable)
    {
        // Geometry element of the origin model part
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        const array_1d<double, 3>& origin_data = (pnode)->FastGetSolutionStepValue(r_origin_variable, 0);
        array_1d<double, 3>& node0_data     = geom[0].FastGetSolutionStepValue(r_destination_variable, 0);
        array_1d<double, 3>& node1_data     = geom[1].FastGetSolutionStepValue(r_destination_variable, 0);
        array_1d<double, 3>& node2_data     = geom[2].FastGetSolutionStepValue(r_destination_variable, 0);
        array_1d<double, 3>& node3_data     = geom[3].FastGetSolutionStepValue(r_destination_variable, 0);

        if (r_origin_variable == HYDRODYNAMIC_FORCE){
            array_1d<double, 3>& node0_drag = geom[0].FastGetSolutionStepValue(HYDRODYNAMIC_REACTION, 0);
            array_1d<double, 3>& node1_drag = geom[1].FastGetSolutionStepValue(HYDRODYNAMIC_REACTION, 0);
            array_1d<double, 3>& node2_drag = geom[2].FastGetSolutionStepValue(HYDRODYNAMIC_REACTION, 0);
            array_1d<double, 3>& node3_drag = geom[3].FastGetSolutionStepValue(HYDRODYNAMIC_REACTION, 0);

            const double fluid_fraction0    = 1 - geom[0].FastGetSolutionStepValue(SOLID_FRACTION, 0);
            const double fluid_fraction1    = 1 - geom[1].FastGetSolutionStepValue(SOLID_FRACTION, 0);
            const double fluid_fraction2    = 1 - geom[2].FastGetSolutionStepValue(SOLID_FRACTION, 0);
            const double fluid_fraction3    = 1 - geom[3].FastGetSolutionStepValue(SOLID_FRACTION, 0);

            const double& node0_volume      = geom[0].FastGetSolutionStepValue(NODAL_AREA, 0);
            const double& node1_volume      = geom[1].FastGetSolutionStepValue(NODAL_AREA, 0);
            const double& node2_volume      = geom[2].FastGetSolutionStepValue(NODAL_AREA, 0);
            const double& node3_volume      = geom[3].FastGetSolutionStepValue(NODAL_AREA, 0);

            const double& node0_density     = geom[0].FastGetSolutionStepValue(DENSITY, 0);
            const double& node1_density     = geom[1].FastGetSolutionStepValue(DENSITY, 0);
            const double& node2_density     = geom[2].FastGetSolutionStepValue(DENSITY, 0);
            const double& node3_density     = geom[3].FastGetSolutionStepValue(DENSITY, 0);

            const double node0_mass_inv     = mParticlesPerDepthDistance / (fluid_fraction0 * node0_volume * node0_density);
            const double node1_mass_inv     = mParticlesPerDepthDistance / (fluid_fraction1 * node1_volume * node1_density);
            const double node2_mass_inv     = mParticlesPerDepthDistance / (fluid_fraction2 * node2_volume * node2_density);
            const double node3_mass_inv     = mParticlesPerDepthDistance / (fluid_fraction3 * node3_volume * node3_density);

            for (unsigned int j= 0; j< TDim; j++){
                double data   = origin_data[j];
                double data_0 = -N[0] * data * node0_mass_inv;
                double data_1 = -N[1] * data * node1_mass_inv;
                double data_2 = -N[2] * data * node2_mass_inv;
                double data_3 = -N[3] * data * node3_mass_inv;
                node0_data[j] += data_0;
                node1_data[j] += data_1;
                node2_data[j] += data_2;
                node3_data[j] += data_3;
                node0_drag[j] += data_0;
                node1_drag[j] += data_1;
                node2_drag[j] += data_2;
                node3_drag[j] += data_3;
            }

        }

        else {
            std::cout << "Variable " << r_origin_variable << " is not supported for transference with constant weights" ;
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void CalculateNodalSolidFractionWithConstantWeighing(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode)
    {
        unsigned int i_nearest_node;
        i_nearest_node = GetNearestNode(N);

        // getting the data of the solution step
        const double& radius         = (pnode)->FastGetSolutionStepValue(RADIUS, 0);
        const double particle_volume = 1.33333333333333333333 * M_PI * mParticlesPerDepthDistance * radius * radius * radius;

        (el_it->GetGeometry())[i_nearest_node].FastGetSolutionStepValue(SOLID_FRACTION, 0) += particle_volume;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void CalculateNodalSolidFractionWithLinearWeighing(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode)
    {
        // Geometry element of the origin model part

        Geometry< Node<3> >& geom = el_it->GetGeometry();
        array_1d<double,4> N_2; // a dummy since we are not interested in its value at the Gauss points
        boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX; // its value is constant over the element so its value on the Gauss point will do
        double element_volume; // a dummy
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N_2, element_volume);

        // getting the data of the solution step
        const double& radius         = (pnode)->FastGetSolutionStepValue(RADIUS, 0);
        const double particle_volume = 1.33333333333333333333 * M_PI * mParticlesPerDepthDistance * radius * radius * radius;

        for (unsigned int inode = 0; inode < N.size(); inode++){
            geom[inode].FastGetSolutionStepValue(SOLID_FRACTION, 0) += N[inode] * particle_volume;
            geom[inode].FastGetSolutionStepValue(SOLID_FRACTION_GRADIENT, 0)[0] += DN_DX(inode, 0) * particle_volume;
            geom[inode].FastGetSolutionStepValue(SOLID_FRACTION_GRADIENT, 0)[1] += DN_DX(inode, 1) * particle_volume;
            geom[inode].FastGetSolutionStepValue(SOLID_FRACTION_GRADIENT, 0)[2] += DN_DX(inode, 2) * particle_volume;
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void CalculateNodalSolidFractionWithLinearWeighingNoGradient(
        Element::Pointer el_it,
        const array_1d<double,4>& NN,
        Node<3>::Pointer pnode)
    {
        // getting the data of the solution step
        const double& radius         = (pnode)->FastGetSolutionStepValue(RADIUS, 0);
        const double particle_volume = 1.33333333333333333333 * M_PI * mParticlesPerDepthDistance * radius * radius * radius;

        for (unsigned int inode = 0; inode < NN.size(); inode++){
            (el_it->GetGeometry())[inode].FastGetSolutionStepValue(SOLID_FRACTION, 0) += NN[inode] * particle_volume;
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

              if (mCouplingType != - 1){
                   ClearVariable(node_it, SOLID_FRACTION);
              }

              array_1d<double, 3>& body_force                = node_it->FastGetSolutionStepValue(BODY_FORCE, 0);
              array_1d<double, 3>& old_hydrodynamic_reaction = node_it->FastGetSolutionStepValue(HYDRODYNAMIC_REACTION, 0);
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

    inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
                                               double& xc, double& yc, double& zc, double& R,
                                               array_1d<double, 3>& N)
    {
        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double x2 = geom[2].X();
        double y2 = geom[2].Y();

        xc = 0.3333333333333333333 * (x0 + x1 + x2);
        yc = 0.3333333333333333333 * (y0 + y1 + y2);
        zc = 0.0;

        double R1 = (xc - x0) * (xc - x0) + (yc - y0) * (yc - y0);
        double R2 = (xc - x1) * (xc - x1) + (yc - y1) * (yc - y1);
        double R3 = (xc - x2) * (xc - x2) + (yc - y2) * (yc - y2);

        R = R1;
        if (R2 > R) R = R2;
        if (R3 > R) R = R3;

        R = 1.01 * sqrt(R);
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
                                               double& xc, double& yc, double& zc, double& R,
                                               array_1d<double,4>& N)
    {
        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double z0 = geom[0].Z();

        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double z1 = geom[1].Z();

        double x2 = geom[2].X();
        double y2 = geom[2].Y();
        double z2 = geom[2].Z();

        double x3 = geom[3].X();
        double y3 = geom[3].Y();
        double z3 = geom[3].Z();

        xc = 0.25 * (x0 + x1 + x2 + x3);
        yc = 0.25 * (y0 + y1 + y2 + y3);
        zc = 0.25 * (z0 + z1 + z2 + z3);

        double R1 = (xc - x0) * (xc - x0) + (yc - y0) * (yc - y0) + (zc - z0) * (zc - z0);
        double R2 = (xc - x1) * (xc - x1) + (yc - y1) * (yc - y1) + (zc - z1) * (zc - z1);
        double R3 = (xc - x2) * (xc - x2) + (yc - y2) * (yc - y2) + (zc - z2) * (zc - z2);
        double R4 = (xc - x3) * (xc - x3) + (yc - y3) * (yc - y3) + (zc - z3) * (zc - z3);

        R = R1;

        if (R2 > R) R = R2;
        if (R3 > R) R = R3;
        if (R4 > R) R = R4;

        R = sqrt(R);
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline double CalculateVol(	const double x0, const double y0,
                                const double x1, const double y1,
                                const double x2, const double y2)
    {
        return 0.5 * ((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0));
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline double CalculateVol(	const double x0, const double y0, const double z0,
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

    inline bool CalculatePosition(Geometry<Node<3> >&geom,
                                  const double xc, const double yc, const double zc,
                                  array_1d<double, 3>& N)
    {
        double x0 = geom[0].X();
        double y0 = geom[0].Y();

        double x1 = geom[1].X();
        double y1 = geom[1].Y();

        double x2 = geom[2].X();
        double y2 = geom[2].Y();

        double area = CalculateVol(x0, y0, x1, y1, x2, y2);
        double inv_area = 0.0;

        if (area == 0.0){
            return false;
        }

        else {
            inv_area = 1.0 / area;
        }

        N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
        N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
        N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;

        // if the xc yc is inside the triangle return true

        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0){
            return true;
        }

        return false;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline bool CalculatePosition(Geometry<Node<3> >&geom,
                                  const double xc, const double yc, const double zc,
                                  array_1d<double, 4>& N)
    {
        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double z0 = geom[0].Z();

        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double z1 = geom[1].Z();

        double x2 = geom[2].X();
        double y2 = geom[2].Y();
        double z2 = geom[2].Z();

        double x3 = geom[3].X();
        double y3 = geom[3].Y();
        double z3 = geom[3].Z();

        double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

        double inv_vol = 0.0;

        if (vol < 0.0000000000001){
            return false;
        }

        else {
            inv_vol = 1.0 / vol;
        }

        N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
        N[1] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;
        N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
        N[3] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;

        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
            N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <=1.0){
            return true;
        }

        return false;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline unsigned int GetNearestNode(const array_1d<double, 4>& N)

    {
      unsigned int i_nearest_node = 0;
      double max                  = N[0];

      for (unsigned int inode = 1; inode < N.size(); ++inode){

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


