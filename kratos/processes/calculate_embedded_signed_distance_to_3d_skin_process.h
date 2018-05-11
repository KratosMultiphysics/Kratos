//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Ruben Zorrilla
//                   Daniel Baumgaertner
//                   Johannes Wolf
//

#if !defined(KRATOS_CALCULATE_EMBEDDED_SIGNED_DISTANCE_TO_3D_SKIN_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_EMBEDDED_SIGNED_DISTANCE_TO_3D_SKIN_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "includes/kratos_flags.h"

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"
#include "utilities/openmp_utils.h"

namespace Kratos {

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

/// Short class definition.
class CalculateEmbeddedSignedDistanceTo3DSkinProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of CalculateEmbeddedSignedDistanceTo3DSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateEmbeddedSignedDistanceTo3DSkinProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    CalculateEmbeddedSignedDistanceTo3DSkinProcess(ModelPart& rThisModelPartStruc, ModelPart& rThisModelPartFluid, bool DiscontinuousDistance = false)
        : mrSkinModelPart(rThisModelPartStruc), mrFluidModelPart(rThisModelPartFluid), mDiscontinuousDistance(DiscontinuousDistance)
    {
    }

    /// Destructor.
    ~CalculateEmbeddedSignedDistanceTo3DSkinProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        // Create a pointer to the discontinuous or continuos distance calculation process
        CalculateDiscontinuousDistanceToSkinProcess<3>::Pointer pdistance_calculator;
        if(mDiscontinuousDistance)
        {
            pdistance_calculator = CalculateDiscontinuousDistanceToSkinProcess<3>::Pointer(
                new CalculateDiscontinuousDistanceToSkinProcess<3>(mrFluidModelPart, mrSkinModelPart));
        }
        else
        {
            pdistance_calculator = CalculateDiscontinuousDistanceToSkinProcess<3>::Pointer(
                new CalculateDistanceToSkinProcess<3>(mrFluidModelPart, mrSkinModelPart));
        }

        // Call the distance calculator methods
        pdistance_calculator->Initialize();
        pdistance_calculator->FindIntersections();
        pdistance_calculator->CalculateDistances(pdistance_calculator->GetIntersections());

        // TODO: Raycasting

        // Distance positive and negative peak values correction
        this->PeakValuesCorrection(); //TODO: Check the correct behaviour of this method once the raycasting has been implemented

        // Compute the embedded velocity
        this->CalculateEmbeddedVelocity(pdistance_calculator->GetIntersections());

        // Call the distance calculation Clear() to delete the intersection data
        pdistance_calculator->Clear();
    }

    virtual void Clear()
    {
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CalculateEmbeddedSignedDistanceTo3DSkinProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateEmbeddedSignedDistanceTo3DSkinProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

    protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void CalculateEmbeddedVelocity(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
    {
        const array_1d<double, 3> aux_zero = ZeroVector(3);

        // #pragma omp parallel for firstprivate(aux_zero)
        for (int k = 0; k < static_cast<int>(mrFluidModelPart.NumberOfElements()); ++k)
        {
            ModelPart::ElementsContainerType::iterator itFluidElement = mrFluidModelPart.ElementsBegin() + k;
            const PointerVector<GeometricalObject>& intersected_skin_elems = rIntersectedObjects[k];

            // Initialize the element EMBEDDED_VELOCITY
            itFluidElement->SetValue(EMBEDDED_VELOCITY, aux_zero);

            // Accumulate the VELOCITY from all the structure conditions that intersect the element
            unsigned int intersection_counter = 0;

            for(auto itSkinElement : intersected_skin_elems)
            {
                array_1d<double,3> emb_vel = (itSkinElement.GetGeometry()[0]).GetSolutionStepValue(VELOCITY);
                emb_vel += (itSkinElement.GetGeometry()[1]).GetSolutionStepValue(VELOCITY);
                emb_vel += (itSkinElement.GetGeometry()[2]).GetSolutionStepValue(VELOCITY);

                itFluidElement->GetValue(EMBEDDED_VELOCITY) += emb_vel/3;
                intersection_counter++;
            }

            // Set the EMBEDDED_VELOCITY as the average of the accumulated values
            if (intersection_counter!=0)
            {
                itFluidElement->GetValue(EMBEDDED_VELOCITY) /= intersection_counter;
            }
        }
    }

    void PeakValuesCorrection()
    {
        // Obtain the maximum and minimum distance values to be set
        double max_distance, min_distance;
        this->SetMaximumAndMinimumDistanceValues(max_distance, min_distance);

        // Bound the distance value in the non splitted nodes
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(mrFluidModelPart.NumberOfNodes()); ++k)
        {
            ModelPart::NodesContainerType::iterator itFluidNode = mrFluidModelPart.NodesBegin() + k;
            if(itFluidNode->IsNot(TO_SPLIT))
            {
                double& rnode_distance = itFluidNode->FastGetSolutionStepValue(DISTANCE);
                rnode_distance = (rnode_distance > 0.0) ? max_distance : min_distance;
            }
        }
    }

    void SetMaximumAndMinimumDistanceValues(double& max_distance, double& min_distance)
    {
        // Flag the nodes belonging to the splitted elements
        for (int k = 0; k < static_cast<int>(mrFluidModelPart.NumberOfElements()); ++k)
        {
            ModelPart::ElementsContainerType::iterator itFluidElement = mrFluidModelPart.ElementsBegin() + k;

            if(itFluidElement->Is(TO_SPLIT))
            {
                Geometry<Node<3>>& rGeom = itFluidElement->GetGeometry();
                for (unsigned int i=0; i<rGeom.size(); ++i)
                {
                    rGeom[i].Set(TO_SPLIT, true);
                }
            }
        }

        // Obtain the maximum and minimum nodal distance values of the nodes flagged as TO_SPLIT
        const unsigned int num_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector nodes_partition;
        OpenMPUtils::DivideInPartitions(mrFluidModelPart.NumberOfNodes(), num_threads, nodes_partition);

        std::vector<double> max_distance_vect(num_threads, 1.0);
        std::vector<double> min_distance_vect(num_threads, 1.0);

        #pragma omp parallel shared(max_distance_vect, min_distance_vect)
        {
            const int k = OpenMPUtils::ThisThread();
            ModelPart::NodeIterator nodes_begin = mrFluidModelPart.NodesBegin() + nodes_partition[k];
            ModelPart::NodeIterator nodes_end   = mrFluidModelPart.NodesBegin() + nodes_partition[k+1];

            double max_local_distance = 1.0;
            double min_local_distance = 1.0;

            for( ModelPart::NodeIterator itFluidNode = nodes_begin; itFluidNode != nodes_end; ++itFluidNode)
            {
                if(itFluidNode->Is(TO_SPLIT))
                {
                    const double node_distance = itFluidNode->FastGetSolutionStepValue(DISTANCE);
                    max_local_distance = (node_distance>max_local_distance) ? node_distance : max_local_distance;
                    min_local_distance = (node_distance<min_local_distance) ? node_distance : min_local_distance;
                }
            }

            max_distance_vect[k] = max_local_distance;
            min_distance_vect[k] = min_local_distance;
        }

        // Reduce to maximum and minimum the thread results
        // Note that MSVC14 does not support max reductions, which are part of OpenMP 3.1
        max_distance = max_distance_vect[0];
        min_distance = min_distance_vect[0];
        for (unsigned int k = 1; k < num_threads; k++)
        {
             max_distance = (max_distance > max_distance_vect[k]) ?  max_distance : max_distance_vect[k];
             min_distance = (min_distance < min_distance_vect[k]) ?  min_distance : min_distance_vect[k];
        }
    }

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

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrSkinModelPart;
    ModelPart& mrFluidModelPart;

    bool mDiscontinuousDistance;

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
    CalculateEmbeddedSignedDistanceTo3DSkinProcess& operator=(CalculateEmbeddedSignedDistanceTo3DSkinProcess const& rOther);

    /// Copy constructor.
    //CalculateEmbeddedSignedDistanceTo3DSkinProcess(CalculateEmbeddedSignedDistanceTo3DSkinProcess const& rOther);

    ///@}

}; // Class CalculateEmbeddedSignedDistanceTo3DSkinProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateEmbeddedSignedDistanceTo3DSkinProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateEmbeddedSignedDistanceTo3DSkinProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_EMBEDDED_SIGNED_DISTANCE_TO_3D_SKIN_PROCESS_H_INCLUDED  defined
