
/*
 * Author: Miguel Angel Celigueta, Akhil
 *
 *  maceli@cimne.upc.edu
 */

#ifndef KRATOS_DEMFEM_VOLUME_COUPLING_UTILITIES_H
#define KRATOS_DEMFEM_VOLUME_COUPLING_UTILITIES_H
// /* External includes */

// System includes

// Project includes
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"

// for variables - application includes
#include "structural_mechanics_application_variables.h"

#include "DEMFEM_volume_coupling_application.h"


namespace Kratos
{
class DEMFEMVolumeCouplingUtilities
{
public:
typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

KRATOS_CLASS_POINTER_DEFINITION(DEMFEMVolumeCouplingUtilities);

/// Default constructor
DEMFEMVolumeCouplingUtilities(){}

/// Destructor
virtual ~DEMFEMVolumeCouplingUtilities(){}

//*************************************************************************************************************** FEM RELATED FUNCTIONS- START
//***************************************************************************************************************


void AssignPointLoads(ModelPart& rFEMModelPart, std::vector<int>& node_ids, double& pointload)
{
    for (ModelPart::NodesContainerType::iterator node_it = rFEMModelPart.NodesBegin(); node_it != rFEMModelPart.NodesEnd(); ++node_it)
    {
        if (std::find(node_ids.begin(), node_ids.end(), node_it->Id()) != node_ids.end()) // assigning point loads - only if force is exerted on FEM.
        {
            array_1d<double, 3> point_load = ZeroVector(3);
            point_load[1] = pointload;
            node_it->FastGetSolutionStepValue(POINT_LOAD) = point_load;
        }
    }
}



void SetNodalCouplingWeightsOnFEMLinearly(ModelPart& rFEMModelPart,double& y_fem_boundary,double& y_dem_boundary, double& tolerance,double& weight_fem_boundary,double& weight_dem_boundary) 
{
//loop through all nodes in modelpart and assign weights to the nodes in the hybrid region
for (ModelPart::NodesContainerType::iterator node_it = rFEMModelPart.NodesBegin(); node_it != rFEMModelPart.NodesEnd(); ++node_it)
{
    if (node_it->Y() >= y_fem_boundary - tolerance && node_it->Y() <= y_dem_boundary + tolerance) // assigning weights to the nodes in the hybrid region
    {
        double weight= weight_fem_boundary + (weight_dem_boundary-weight_fem_boundary)*(node_it->Y()-y_fem_boundary)/(y_dem_boundary-y_fem_boundary);
        node_it->FastGetSolutionStepValue(NODAL_COUPLING_WEIGHT)=weight;
    }
    else if (node_it->Y() <= y_fem_boundary + tolerance && node_it->Y() >= y_dem_boundary - tolerance) // assigning weights to the nodes in the hybrid region
    {
        double weight= weight_dem_boundary + (weight_fem_boundary-weight_dem_boundary)*(node_it->Y()-y_dem_boundary)/(y_fem_boundary-y_dem_boundary);
        node_it->FastGetSolutionStepValue(NODAL_COUPLING_WEIGHT)=weight;
    }
    else
    {
        node_it->FastGetSolutionStepValue(NODAL_COUPLING_WEIGHT)=0.0; // assigning weights to the nodes in the non-hybrid region
    }
}
}

// void SetNodalCouplingWeightsOnFEMLinearly(ModelPart& rFEMModelPart, double& y_fem_boundary, double& y_dem_boundary, double& tolerance, double& weight_fem_boundary, double& weight_dem_boundary)
// {
//     // loop through all nodes in modelpart and assign weights to the nodes in the hybrid region
//     for (ModelPart::NodesContainerType::iterator node_it = rFEMModelPart.NodesBegin(); node_it != rFEMModelPart.NodesEnd(); ++node_it)
//     {
//         double y_pos = node_it->Y();

//         if (y_pos >= y_fem_boundary - tolerance && y_pos <= y_dem_boundary + tolerance) // assigning weights to the nodes in the hybrid region
//         {
//             double normalized_y = (y_pos - y_fem_boundary) / (y_dem_boundary - y_fem_boundary); // normalized position between boundaries
//             double weight = weight_fem_boundary + (weight_dem_boundary - weight_fem_boundary) * (3 * std::pow(normalized_y, 2) - 2 * std::pow(normalized_y, 3)); // cubic interpolation
//             node_it->FastGetSolutionStepValue(NODAL_COUPLING_WEIGHT) = weight;
//         }
//         else if (y_pos <= y_fem_boundary + tolerance && y_pos >= y_dem_boundary - tolerance) // assigning weights to the nodes in the hybrid region
//         {
//             double normalized_y = (y_pos - y_dem_boundary) / (y_fem_boundary - y_dem_boundary); // normalized position between boundaries
//             double weight = weight_dem_boundary + (weight_fem_boundary - weight_dem_boundary) * (3 * std::pow(normalized_y, 2) - 2 * std::pow(normalized_y, 3)); // cubic interpolation
//             node_it->FastGetSolutionStepValue(NODAL_COUPLING_WEIGHT) = weight;
//         }
//         else
//         {
//             node_it->FastGetSolutionStepValue(NODAL_COUPLING_WEIGHT) = 0.0; // assigning weights to the nodes in the non-hybrid region
//         }
//     }
// }





void CalculateDisplacementDifference(ModelPart& rFEMModelPart, double& dt)
{
    for (ModelPart::NodesContainerType::iterator node_it = rFEMModelPart.NodesBegin(); node_it != rFEMModelPart.NodesEnd(); ++node_it)
    {
        double total_mass = node_it->FastGetSolutionStepValue(NODAL_MAUX);
        if (total_mass != 0) // check for hybrid region
        {
            array_1d<double, 3> zero_array = ZeroVector(3);
            node_it->FastGetSolutionStepValue(POINT_LOAD) = zero_array; // nodal coupling forces set to zero before every iteration
            node_it->FastGetSolutionStepValue(DEMFEM_VOLUME_COUPLING_FORCE) = zero_array; // force to map to dem (set to zero before every iteration)
            // array_1d<double, 3> displacement_dem = node_it->FastGetSolutionStepValue(DISPLACEMENT_MULTIPLIED_MASS) / total_mass; // updated lagrange method-> calculating homogenised displacement
            // node_it->FastGetSolutionStepValue(PENALIZE_DISPLACEMENT) += (displacement_dem - node_it->FastGetSolutionStepValue(DISPLACEMENT)) ; // updated lagrange method : calcualting displacement difference
            // // condition checking if the above displacement difference (PENALIZE_DISPLACEMENT) is greater than 0.1*element size then set the displacement difference to 0.1*element size
            // if (norm_2(node_it->FastGetSolutionStepValue(PENALIZE_DISPLACEMENT)) > 0.001 * 0.02) 
            // {
            //     node_it->FastGetSolutionStepValue(PENALIZE_DISPLACEMENT) = 0.001 * 0.02 * node_it->FastGetSolutionStepValue(PENALIZE_DISPLACEMENT) / norm_2(node_it->FastGetSolutionStepValue(PENALIZE_DISPLACEMENT));
            // }
            //print value of nodal_H
            //std::cout<<"NODAL_H"<<node_it->FastGetSolutionStepValue(NODAL_H)<<std::endl;


            array_1d<double, 3> velocity_dem = node_it->FastGetSolutionStepValue(VELOCITY_MULTIPLIED_MASS) / total_mass; // updated lagrange method-> calculating homogenised velocity
            node_it->FastGetSolutionStepValue(PENALIZE_DISPLACEMENT) += (velocity_dem - node_it->FastGetSolutionStepValue(VELOCITY))*dt ; // updated lagrange method : calcualting displacement difference in velocity method
       }
    }

}
        




// void CalculateNodalCouplingForces(ModelPart& rFEMModelPart, double& penalty_max)
// {
//     for (ModelPart::ElementsContainerType::iterator elem_it = rFEMModelPart.ElementsBegin(); elem_it != rFEMModelPart.ElementsEnd(); ++elem_it)
//     {   
//         //std::vector<double> V = elem_it->CalculateOnIntegrationPoints(INTEGRATION_WEIGHT,rFEMModelPart.ProcessInfo);// V is an array of integration weights 
//         // to check if any of the nodes in the element is in the hybrid region by checking elem_it->GetGeometry()[n].FastGetSolutionStepValue(NODAL_MAUX) is 0 or not
//         if(elem_it->GetGeometry()[0].FastGetSolutionStepValue(NODAL_MAUX) != 0)
//         {
//             for (unsigned int i = 0; i < elem_it->GetGeometry().IntegrationPointsNumber(); i++)
//             {
//                 double J = (elem_it->GetGeometry().DeterminantOfJacobian(i));
//                 // double w = V[i] / J;
//                 double w = 1.0;
//                 const Kratos::Matrix shape_functions = elem_it->GetGeometry().ShapeFunctionsValues();
//                 for (unsigned int n = 0; n < elem_it->GetGeometry().size(); n++)
//                 {
//                     for (unsigned int m = 0; m < elem_it->GetGeometry().size(); m++)
//                     {
//                         double vol = penalty_max * w * J * shape_functions(i,n) * shape_functions(i,m);
//                         elem_it->GetGeometry()[n].FastGetSolutionStepValue(POINT_LOAD) += (vol * elem_it->GetGeometry()[n].FastGetSolutionStepValue(PENALIZE_DISPLACEMENT) );
//                     }
//                 }
//             }
//         }
//     }
// }

void CalculateNodalCouplingForces(ModelPart& rFEMModelPart, double& penalty_max)
{
    for (ModelPart::ElementsContainerType::iterator elem_it = rFEMModelPart.ElementsBegin(); elem_it != rFEMModelPart.ElementsEnd(); ++elem_it)
    {   
        // Check if any of the nodes in the element is in the hybrid region
        if (elem_it->GetGeometry()[0].FastGetSolutionStepValue(NODAL_MAUX) != 0)
        {
            for (unsigned int i = 0; i < elem_it->GetGeometry().IntegrationPointsNumber(); i++)
            {
                double J = elem_it->GetGeometry().DeterminantOfJacobian(i);
                // Adjust this weight if needed. Assuming V[i] / J is omitted here
                double w = 1.0;
                const Kratos::Matrix shape_functions = elem_it->GetGeometry().ShapeFunctionsValues();
                
                for (unsigned int n = 0; n < elem_it->GetGeometry().size(); n++)
                {
                    for (unsigned int m = 0; m < elem_it->GetGeometry().size(); m++)
                    {
                        // Define penalties for each direction
                        std::array<double, 3> penalties = {penalty_max, penalty_max, penalty_max};

                        // Loop over directions and apply specific penalty for each
                        for (size_t dir = 0; dir < 3; ++dir)
                        {
                            double vol = penalties[dir] * w * J * shape_functions(i, n) * shape_functions(i, m);
                            elem_it->GetGeometry()[n].FastGetSolutionStepValue(POINT_LOAD)[dir] += vol * elem_it->GetGeometry()[n].FastGetSolutionStepValue(PENALIZE_DISPLACEMENT)[dir];
                        }
                    }
                }
            }
        }
    }
}


void CalculateNodalDEMCouplingForces(ModelPart& rFEMModelPart)
{
    for (ModelPart::NodesContainerType::iterator node_it = rFEMModelPart.NodesBegin(); node_it != rFEMModelPart.NodesEnd(); ++node_it)
    {
        double total_mass = node_it->FastGetSolutionStepValue(NODAL_MAUX);
        if (total_mass != 0) // check for hybrid region
        {
             node_it->FastGetSolutionStepValue(DEMFEM_VOLUME_COUPLING_FORCE) = -1 * node_it->FastGetSolutionStepValue(POINT_LOAD) / total_mass; // force to map to dem
        }
    }
}

//*************************************************************************************************************** FEM RELATED FUNCTIONS- END
//*************************************************************************************************************** DEM RELATED FUNCTIONS- START


void CalculateMomentum(ModelPart& rDEMModelPart)
{
    for (ModelPart::NodesContainerType::iterator node_it = rDEMModelPart.NodesBegin(); node_it != rDEMModelPart.NodesEnd(); ++node_it)
    {
        node_it->FastGetSolutionStepValue(DISPLACEMENT_MULTIPLIED_MASS) = node_it->FastGetSolutionStepValue(NODAL_MASS) * node_it->FastGetSolutionStepValue(DISPLACEMENT); //FOR SENDING DISPLACEMENT
        node_it->FastGetSolutionStepValue(VELOCITY_MULTIPLIED_MASS) = node_it->FastGetSolutionStepValue(NODAL_MASS) * node_it->FastGetSolutionStepValue(VELOCITY); // FOR SENDING VELOCITY
    }
}


void CalculateDEMForces(ModelPart& rDEMModelPart)
{
    for (ModelPart::NodesContainerType::iterator node_it = rDEMModelPart.NodesBegin(); node_it != rDEMModelPart.NodesEnd(); ++node_it)
    {
        double particle_weight = node_it->FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT);
        if (particle_weight == 0)
        {
            node_it->FastGetSolutionStepValue(PARTICLE_COUPLING_WEIGHT) = 1;
        }
        else
        {
            node_it->FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) = node_it->FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) * node_it->FastGetSolutionStepValue(NODAL_MASS); // if nodal mass needs to be multiplied with external applied force
        }
    }
}
//*************************************************************************************************************** DEM RELATED FUNCTIONS- END
//***************************************************************************************************************




/// Turn back information as a stemplate<class T, std::size_t dim> tring.

virtual std::string Info() const
{
    return "";
}

/// Print information about this object.

virtual void PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.

virtual void PrintData(std::ostream& rOStream) const
{
}

protected:

private:

/// Assignment operator
DEMFEMVolumeCouplingUtilities & operator=(DEMFEMVolumeCouplingUtilities const& rOther);


///@}

}; // Class DEMFEMVolumeCouplingUtilities

}  // namespace Python.

#endif // DEMFEM_VOLUME_COUPLING_UTILITIES_H