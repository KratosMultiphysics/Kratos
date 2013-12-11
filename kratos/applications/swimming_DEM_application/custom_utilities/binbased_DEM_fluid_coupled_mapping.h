//
//   Project Name:        Kratos
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2012-03-08 08:56:42 $
//
//

#if !defined(KRATOS_BINBASED_DEM_FLUID_COUPLED_MAPPING )
#define  KRATOS_BINBASED_DEM_FLUID_COUPLED_MAPPING

// /* External includes */
// #include "boost/smart_ptr.hpp"

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "utilities/timer.h"

// #include "geometries/tetrahedra_3d_4.h"

//#include "PFEM_application.h"

//Database includes
#include "spatial_containers/spatial_containers.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"

//iNCLUDE THE DRAG UTILITIED TO CALCULATE THE SEEPAGE DRAG
// #include "custom_utilities/drag_utilities.h"

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

/// This class allows the interpolation between non-matching meshes in 2D and 3D. it is designed for DEM-CFD coupling problems
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
      * @param rOrigin_ModelPart: the model part  all the variable should be taken from
      * @param rDestination_ModelPart: the destination model part where we want to know the values of the variables
      * @param bin_of_objects_fluid: precomputed bin of objects (elelments of the fluid mesh). It is to be constructed separately
      * @see binbased_nodes_in_element_locator
    */

    //***************************************************************************************************************
    //***************************************************************************************************************

    // The fluid data are a weighted: data_to_project = alpha * new_data + (1 - alpha) * old_data

    void InterpolateFromFluidMesh(
        ModelPart& rfluid_model_part,
        ModelPart& rdem_model_part,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
        const double alpha)
    {
        KRATOS_TRY

        //Clear all the variables to be mapped
        for (ModelPart::NodesContainerType::iterator node_it = rdem_model_part.NodesBegin();
                node_it != rdem_model_part.NodesEnd(); ++node_it){

            ClearVariables(node_it, FLUID_VEL_PROJECTED);
            ClearVariables(node_it, PRESSURE_GRAD_PROJECTED);
            ClearVariables(node_it, FLUID_DENSITY_PROJECTED);
            ClearVariables(node_it, FLUID_VISCOSITY_PROJECTED);
        }

        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rdem_model_part.Nodes().size();

        #pragma omp parallel for firstprivate(results, N)
        for (int i = 0; i < nparticles; i++){
            ModelPart::NodesContainerType::iterator iparticle = rdem_model_part.NodesBegin() + i;
            Node < 3 > ::Pointer pparticle = *(iparticle.base());

            if (!pparticle->IsFixed(VELOCITY_X)){
                typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
                Element::Pointer pelement;

                // looking for the fluid element in which the DEM node falls
                bool is_found = bin_of_objects_fluid.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);

                // interpolating the variables

                if (is_found){
                    //Interpolate(el_it,  N, *it_found , rOriginVariable , rDestinationVariable, Present/(Present + Past)Coeff);
                    Interpolate(pelement, N, pparticle, DENSITY, FLUID_DENSITY_PROJECTED, alpha);
                    Interpolate(pelement, N, pparticle, VELOCITY, FLUID_VEL_PROJECTED, alpha);
                    Interpolate(pelement, N, pparticle, PRESSURE_GRADIENT, PRESSURE_GRAD_PROJECTED, alpha);
                    Interpolate(pelement, N, pparticle, VISCOSITY, FLUID_VISCOSITY_PROJECTED, alpha);
                }
            }
        }

        KRATOS_CATCH("")
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    // The current fluid data are projected to the DEM mesh
    void InterpolateFromNewestFluidMesh(
        ModelPart& rfluid_model_part,
        ModelPart& rdem_model_part,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid)
    {
        KRATOS_TRY

        // resetting all variables to be interpolated

        for (ModelPart::NodesContainerType::iterator node_it = rdem_model_part.NodesBegin();
             node_it != rdem_model_part.NodesEnd(); ++node_it){
            ClearVariables(node_it, FLUID_VEL_PROJECTED);
            ClearVariables(node_it, PRESSURE_GRAD_PROJECTED);
            ClearVariables(node_it, FLUID_DENSITY_PROJECTED);
            ClearVariables(node_it, FLUID_VISCOSITY_PROJECTED);
            ClearVariables(node_it, SOLID_FRACTION_PROJECTED);
        }

        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rdem_model_part.Nodes().size();

        #pragma omp parallel for firstprivate(results, N)
        for (int i = 0; i < nparticles; i++){
            ModelPart::NodesContainerType::iterator iparticle = rdem_model_part.NodesBegin() + i;
            Node < 3 > ::Pointer pparticle = *(iparticle.base());

            if (!pparticle->IsFixed(VELOCITY_X)){
                typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
                Element::Pointer pelement;

                // looking for the fluid element in which the DEM node falls
                bool is_found = bin_of_objects_fluid.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);

                // interpolating variables

                if (is_found){
                    //Interpolate(el_it,  N, *it_found , rOriginVariable , rDestinationVariable);
                    Interpolate(pelement, N, pparticle, DENSITY, FLUID_DENSITY_PROJECTED);
                    Interpolate(pelement, N, pparticle, VELOCITY, FLUID_VEL_PROJECTED);
                    Interpolate(pelement, N, pparticle, PRESSURE_GRADIENT, PRESSURE_GRAD_PROJECTED);
                    Interpolate(pelement, N, pparticle, VISCOSITY, FLUID_VISCOSITY_PROJECTED);
                    Interpolate(pelement, N, pparticle, SOLID_FRACTION, SOLID_FRACTION_PROJECTED);
                }

            }

        }

        KRATOS_CATCH("")
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    /// Interpolate form the DEM  to the fluid mesh
    /**
      * @param rOrigin_ModelPart: the model part  all the variable should be taken from
      * @param rDestination_ModelPart: the destination model part where we want to know the values of the variables
      * @param bin_of_nodes_fluid: precomputed bin of nodes of the fluid mesh. It is to be constructed separately @see binbased_nodes_in_element_locator
    */

    //***************************************************************************************************************
    //***************************************************************************************************************

    // The current DEM data are projected to the fluid mesh
    void InterpolateFromDEMMesh(
        ModelPart& rdem_model_part,
        ModelPart& rfluid_model_part,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid) //this is a bin of objects which contains the FLUID model part
     {
        KRATOS_TRY
        const int n_fluid_nodes = rfluid_model_part.Nodes().size();

        // resetting the variables to be mapped

        for (int i = 0; i < n_fluid_nodes; i++){
            ModelPart::NodesContainerType::iterator inode = rfluid_model_part.NodesBegin() + i;
            ClearVariables(inode, SOLID_FRACTION);
            array_1d<double,3>& body_force              = inode->FastGetSolutionStepValue(BODY_FORCE, 0);
            const array_1d<double,3>& old_drag_reaction = inode->FastGetSolutionStepValue(DRAG_REACTION, 0);
            body_force -= old_drag_reaction;
            ClearVariables(inode, DRAG_REACTION);
        }

        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rdem_model_part.Nodes().size();

        #pragma omp parallel for firstprivate(results, N)
        for (int i = 0; i < nparticles; i++){
            ModelPart::NodesContainerType::iterator iparticle = rdem_model_part.NodesBegin() + i;
            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
            Element::Pointer pelement;

            // looking for the fluid element in which the DEM node falls
            bool is_found = bin_of_objects_fluid.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);

            // interpolating variables

            if (is_found) {

                if (mCouplingType == 0){
                    TransferWithConstantWeighing(pelement, N, pparticle, BODY_FORCE, DRAG_FORCE);
                    CalculateNodalSolidFractionWithConstantWeighing(pelement, N, pparticle);
                }

                else if (mCouplingType == 1){
                    TransferWithLinearWeighing(pelement, N, pparticle, BODY_FORCE, DRAG_FORCE);
                    CalculateNodalSolidFractionWithConstantWeighing(pelement, N, pparticle);
                }

                else if (mCouplingType == 2){
                    TransferWithLinearWeighing(pelement, N, pparticle, BODY_FORCE, DRAG_FORCE);
                    CalculateNodalSolidFractionWithLinearWeighing(pelement, N, pparticle);
                }

            }

        }

        for (int i = 0; i < n_fluid_nodes; i++){
            ModelPart::NodesContainerType::iterator inode = rfluid_model_part.NodesBegin() + i;
            double& solid_fraction = inode->FastGetSolutionStepValue(SOLID_FRACTION, 0);
            solid_fraction /= inode->FastGetSolutionStepValue(NODAL_AREA, 0);

            if (solid_fraction > mMaxSolidFraction){
                solid_fraction = mMaxSolidFraction;
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

      for (int i = 0; i < n_dem_elements; i++){
          ModelPart::ElementsContainerType::iterator ielem = rdem_model_part.ElementsBegin() + i;
          Geometry< Node<3> >& geom = ielem->GetGeometry();
          double& reynolds_number = geom[0].FastGetSolutionStepValue(REYNOLDS_NUMBER, 0);
          ielem->Calculate(REYNOLDS_NUMBER, reynolds_number, r_current_process_info);
      }

      for (int i = 0; i < n_fluid_nodes; i++){
          ModelPart::NodesContainerType::iterator inode = rfluid_model_part.NodesBegin() + i;
          double fluid_fraction                         = 1 - inode->FastGetSolutionStepValue(SOLID_FRACTION, 0);
          const array_1d<double,3>& darcy_vel           = inode->FastGetSolutionStepValue(VELOCITY, 0);
          array_1d<double,3>& space_averaged_fluid_vel  = inode->FastGetSolutionStepValue(MESH_VELOCITY1, 0);
          space_averaged_fluid_vel                      = darcy_vel / fluid_fraction;
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

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
                                               double& xc, double& yc, double& zc, double& R, array_1d<double,3>& N)
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
                                               double& xc, double& yc, double& zc, double& R, array_1d<double,4>& N)
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
                                  const double xc, const double yc, const double zc, array_1d<double,3>& N)
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
                                  const double xc, const double yc, const double zc, array_1d<double,4>& N)
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

    inline unsigned int GetNearestNode(const array_1d<double,4>& N)

    {
      unsigned int i_nearest_node = 0;
      double max                  = N[0];

      for (unsigned int inode = 1; inode < N.size(); inode++){

          if (N[inode] > max) {
              max = N[inode];
              i_nearest_node = inode;
          }

      }

      return i_nearest_node;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    // project an array1D (2Dversion)
    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double,3>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double,3> >& rOriginVariable,
        Variable<array_1d<double,3> >& rDestinationVariable)
    {

        //Geometry element of the rOrigin_ModelPart
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        //getting the data of the solution step
        array_1d<double,3>& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable, 0);
        const array_1d<double,3>& velocity = (pnode)->FastGetSolutionStepValue(VELOCITY, 0);
        const array_1d<double,3>& node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable, 0);
        const array_1d<double,3>& node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable, 0);
        const array_1d<double,3>& node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable, 0);

        //copying this data in the position of the vector we are interested in
        for (unsigned int j= 0; j< TDim; j++){
            step_data[j] = N[0] * node0_data[j] + N[1] * node1_data[j] + N[2] * node2_data[j];
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    //projecting an array1D 3Dversion

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double,3> >& rOriginVariable,
        Variable<array_1d<double,3> >& rDestinationVariable)
    {
          //Geometry element of the rOrigin_ModelPart
          Geometry< Node < 3 > >& geom = el_it->GetGeometry();

          //getting the data of the solution step
          array_1d<double, 3 > & step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable);

          const array_1d<double, 3 > & node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable);
          const array_1d<double, 3 > & node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable);
          const array_1d<double, 3 > & node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable);
          const array_1d<double, 3 > & node3_data = geom[3].FastGetSolutionStepValue(rOriginVariable);

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
        Variable<array_1d<double,3> >& rOriginVariable,
        Variable<array_1d<double,3> >& rDestinationVariable,
        double alpha)
    {

        //Geometry element of the rOrigin_ModelPart
        Geometry< Node < 3 > >& geom = el_it->GetGeometry();

        //getting the data of the solution step
        array_1d<double, 3 > & step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable, 0);

        const array_1d<double, 3 > & node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable, 0);
        const array_1d<double, 3 > & node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable, 0);
        const array_1d<double, 3 > & node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable, 0);
        const array_1d<double, 3 > & node3_data = geom[3].FastGetSolutionStepValue(rOriginVariable, 0);

        const array_1d<double, 3 > & node0_data_prev = geom[0].FastGetSolutionStepValue(rOriginVariable, 1);
        const array_1d<double, 3 > & node1_data_prev = geom[1].FastGetSolutionStepValue(rOriginVariable, 1);
        const array_1d<double, 3 > & node2_data_prev = geom[2].FastGetSolutionStepValue(rOriginVariable, 1);
        const array_1d<double, 3 > & node3_data_prev = geom[3].FastGetSolutionStepValue(rOriginVariable, 1);

        for (unsigned int j = 0; j < TDim; j++) {
            step_data[j] = alpha * (N[0] * node0_data[j]      + N[1] * node1_data[j]      + N[2] * node2_data[j]      + N[3] * node3_data[j]) +
                   (1.0 - alpha) * (N[0] * node0_data_prev[j] + N[1] * node1_data_prev[j] + N[2] * node2_data_prev[j] + N[3] * node3_data_prev[j]);
        }

     }

    //***************************************************************************************************************
    //***************************************************************************************************************

    //projecting a scalar 2Dversion
    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double,3>& N,
        Node<3>::Pointer pnode,
        Variable<double>& rOriginVariable,
        Variable<double>& rDestinationVariable)
    {

        Geometry< Node < 3 > >& geom = el_it->GetGeometry();

        //getting the data of the solution step
        double& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable, 0);
        const double node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable, 0);
        const double node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable, 0);
        const double node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable, 0);

        //copying this data in the position of the vector we are interested in
        step_data = N[0] * node0_data + N[1] * node1_data + N[2] * node2_data;
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    //projecting a scalar 3Dversion
    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<double>& rOriginVariable,
        Variable<double>& rDestinationVariable)
    {
        //Geometry element of the rOrigin_ModelPart
        Geometry< Node < 3 > >& geom = el_it->GetGeometry();

        //getting the data of the solution step
        double& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable, 0);
        const double node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable, 0);
        const double node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable, 0);
        const double node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable, 0);
        const double node3_data = geom[3].FastGetSolutionStepValue(rOriginVariable, 0);
        step_data = N[0] * node0_data + N[1] * node1_data + N[2] * node2_data + N[3] * node3_data;

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void Interpolate(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<double>& rOriginVariable,
        Variable<double>& rDestinationVariable,
        double alpha)
    {
        // Geometry element of the rOrigin_ModelPart
        Geometry< Node < 3 > >& geom = el_it->GetGeometry();

        // getting the data of the solution step
        double& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable, 0);

        const double node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable, 0);
        const double node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable, 0);
        const double node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable, 0);
        const double node3_data = geom[3].FastGetSolutionStepValue(rOriginVariable, 0);

        const double node0_data_prev = geom[0].FastGetSolutionStepValue(rOriginVariable, 1);
        const double node1_data_prev = geom[1].FastGetSolutionStepValue(rOriginVariable, 1);
        const double node2_data_prev = geom[2].FastGetSolutionStepValue(rOriginVariable, 1);
        const double node3_data_prev = geom[3].FastGetSolutionStepValue(rOriginVariable, 1);

        step_data = alpha * (N[0] * node0_data + N[1] * node1_data + N[2] * node2_data + N[3] * node3_data) +
        (1.0 - alpha) * (N[0] * node0_data_prev + N[1] * node1_data_prev + N[2] * node2_data_prev + N[3] * node3_data_prev);

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    //2Dversion
    void TransferWithConstantWeighing(
        Element::Pointer el_it,
        const array_1d<double,3>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double,3> >& rDestinationVariable,
        Variable<array_1d<double,3> >& rOriginVariable
        )

    {
        //Geometry element of the rOrigin_ModelPart
        Geometry< Node<3> >& geom             = el_it->GetGeometry();
        unsigned int i_nearest_node           = GetNearestNode(N);
        const array_1d<double,3>& origin_data = (pnode)->FastGetSolutionStepValue(rOriginVariable, 0);
        array_1d<double,3>& destination_data  = geom[i_nearest_node].FastGetSolutionStepValue(rDestinationVariable, 0);

        if (rOriginVariable == DRAG_FORCE){
            const double density                      = geom[i_nearest_node].FastGetSolutionStepValue(DENSITY, 0);
            const double nodal_volume                 = geom[i_nearest_node].FastGetSolutionStepValue(NODAL_AREA, 0);
            array_1d<double,3>& old_drag_contribution = geom[i_nearest_node].FastGetSolutionStepValue(DRAG_REACTION, 0);
            const double nodal_mass_inv               = mParticlesPerDepthDistance / (density * nodal_volume);

            for (unsigned int j= 0; j< TDim; j++){
                array_1d<double,3> origin_nodal_contribution = - nodal_mass_inv * origin_data;
                destination_data += origin_nodal_contribution;
                old_drag_contribution += origin_nodal_contribution;
            }

        }

        else {
            std::cout << "Variable " << rOriginVariable << " is not supported for transference with constant weights";
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void TransferWithLinearWeighing(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double,3> >& rDestinationVariable,
        Variable<array_1d<double,3> >& rOriginVariable)
        
    {
        //Geometry element of the rOrigin_ModelPart
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        const array_1d<double,3>& origin_data = (pnode)->FastGetSolutionStepValue(rOriginVariable, 0);
        array_1d<double,3>& node0_data = geom[0].FastGetSolutionStepValue(rDestinationVariable, 0);
        array_1d<double,3>& node1_data = geom[1].FastGetSolutionStepValue(rDestinationVariable, 0);
        array_1d<double,3>& node2_data = geom[2].FastGetSolutionStepValue(rDestinationVariable, 0);
        array_1d<double,3>& node3_data = geom[3].FastGetSolutionStepValue(rDestinationVariable, 0);

        if (rOriginVariable == DRAG_FORCE){
            array_1d<double,3>& node0_drag = geom[0].FastGetSolutionStepValue(DRAG_REACTION, 0);
            array_1d<double,3>& node1_drag = geom[1].FastGetSolutionStepValue(DRAG_REACTION, 0);
            array_1d<double,3>& node2_drag = geom[2].FastGetSolutionStepValue(DRAG_REACTION, 0);
            array_1d<double,3>& node3_drag = geom[3].FastGetSolutionStepValue(DRAG_REACTION, 0);
            const double fluid_fraction0   = 1 - geom[0].FastGetSolutionStepValue(SOLID_FRACTION, 0);
            const double fluid_fraction1   = 1 - geom[1].FastGetSolutionStepValue(SOLID_FRACTION, 0);
            const double fluid_fraction2   = 1 - geom[2].FastGetSolutionStepValue(SOLID_FRACTION, 0);
            const double fluid_fraction3   = 1 - geom[3].FastGetSolutionStepValue(SOLID_FRACTION, 0);
            const double& node0_volume     = geom[0].FastGetSolutionStepValue(NODAL_AREA, 0);
            const double& node1_volume     = geom[1].FastGetSolutionStepValue(NODAL_AREA, 0);
            const double& node2_volume     = geom[2].FastGetSolutionStepValue(NODAL_AREA, 0);
            const double& node3_volume     = geom[3].FastGetSolutionStepValue(NODAL_AREA, 0);
            const double& node0_density    = geom[0].FastGetSolutionStepValue(DENSITY, 0);
            const double& node1_density    = geom[1].FastGetSolutionStepValue(DENSITY, 0);
            const double& node2_density    = geom[2].FastGetSolutionStepValue(DENSITY, 0);
            const double& node3_density    = geom[3].FastGetSolutionStepValue(DENSITY, 0);
            const double node0_mass_inv    = mParticlesPerDepthDistance / (fluid_fraction0 * node0_volume * node0_density);
            const double node1_mass_inv    = mParticlesPerDepthDistance / (fluid_fraction1 * node1_volume * node1_density);
            const double node2_mass_inv    = mParticlesPerDepthDistance / (fluid_fraction2 * node2_volume * node2_density);
            const double node3_mass_inv    = mParticlesPerDepthDistance / (fluid_fraction3 * node3_volume * node3_density);

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
            std::cout << "Variable " << rOriginVariable << " is not supported for transference with constant weights" ;
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void CalculateNodalSolidFractionWithConstantWeighing(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode)

    {
        //Geometry element of the rOrigin_ModelPart
        Geometry< Node<3> >& geom = el_it->GetGeometry();
        unsigned int i_nearest_node;
        i_nearest_node = GetNearestNode(N);

        //getting the data of the solution step
        const double& radius         = (pnode)->FastGetSolutionStepValue(RADIUS, 0);
        const double particle_volume = 1.33333333333333333333 * M_PI * mParticlesPerDepthDistance * radius * radius * radius;

        geom[i_nearest_node].FastGetSolutionStepValue(SOLID_FRACTION, 0) += particle_volume;

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    void CalculateNodalSolidFractionWithLinearWeighing(
        Element::Pointer el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode)

    {
        //Geometry element of the rOrigin_ModelPart

        Geometry< Node<3> >& geom = el_it->GetGeometry();
        array_1d<double,4> N_2; // a dummy since we are nbot interested in its value at the Gauss points
        boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX; // its value is constant over the element so its value on the Gauss point will do
        double element_volume; // a dummy
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N_2, element_volume);

        //getting the data of the solution step
        const double& radius         = (pnode)->FastGetSolutionStepValue(RADIUS, 0);
        const double particle_volume = 1.33333333333333333333 * M_PI * mParticlesPerDepthDistance * radius * radius * radius;

        for (unsigned int inode = 0; inode < N.size(); inode++){
            geom[inode].FastGetSolutionStepValue(SOLID_FRACTION, 0) += N[inode] * particle_volume;
//            geom[inode].FastGetSolutionStepValue(SOLID_FRACTION_GRADIENT, 0)[0] += DN_DX(inode, 0) * particle_volume;
//            geom[inode].FastGetSolutionStepValue(SOLID_FRACTION_GRADIENT, 0)[1] += DN_DX(inode, 1) * particle_volume;
//            geom[inode].FastGetSolutionStepValue(SOLID_FRACTION_GRADIENT, 0)[2] += DN_DX(inode, 2) * particle_volume;
        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline void Clear(ModelPart::NodesContainerType::iterator node_it,  int step_data_size )
    {
        unsigned int buffer_size = node_it->GetBufferSize();

        for (unsigned int step = 0; step<buffer_size; step++){
            //getting the data of the solution step
            double* step_data = (node_it)->SolutionStepData().Data(step);

            //copying this data in the position of the vector we are interested in

            for (int j= 0; j< step_data_size; j++){
                step_data[j] = 0.0;
            }

        }

    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it , Variable<array_1d<double,3> >& rVariable)
    {
        array_1d<double, 3>& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);        
        noalias(Aux_var) = ZeroVector(3);
    }

    //***************************************************************************************************************
    //***************************************************************************************************************

    inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it,  Variable<double>& rVariable)
    {
        double& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
        Aux_var = 0.0;
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


