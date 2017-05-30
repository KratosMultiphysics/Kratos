//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Ruben Zorrilla
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
    virtual ~CalculateEmbeddedSignedDistanceTo3DSkinProcess()
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

    virtual void Execute() override
    {
        // Create a pointer to the discontinuous or continuos distance calculation process
        CalculateDiscontinuousDistanceToSkinProcess::Pointer pdistance_calculator;
        if(mDiscontinuousDistance)
        {
            pdistance_calculator = CalculateDiscontinuousDistanceToSkinProcess::Pointer(
                new CalculateDiscontinuousDistanceToSkinProcess(mrFluidModelPart, mrSkinModelPart));
        }
        else
        {
            pdistance_calculator = CalculateDiscontinuousDistanceToSkinProcess::Pointer(
                new CalculateDistanceToSkinProcess(mrFluidModelPart, mrSkinModelPart));
        }

        // Call the distance calculator methods
        pdistance_calculator->Initialize();
        pdistance_calculator->FindIntersections();
        pdistance_calculator->ComputeDistances(pdistance_calculator->GetIntersections());

        // Compute the embedded velocity
        this->ComputeEmbeddedVelocity(pdistance_calculator->GetIntersections());

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
    virtual std::string Info() const
    {
        return "CalculateEmbeddedSignedDistanceTo3DSkinProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "CalculateEmbeddedSignedDistanceTo3DSkinProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    void ComputeEmbeddedVelocity(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
    {
        #pragma omp for
        for (int k = 0; k < static_cast<int>(mrFluidModelPart.NumberOfElements()); k++)
        {
            ModelPart::ElementsContainerType::iterator itFluidElement = mrFluidModelPart.ElementsBegin() + k;
            PointerVector<GeometricalObject> intersected_skin_elems = rIntersectedObjects[itFluidElement->Id()];

            unsigned int intersection_counter = 0;

            // Accumulate the VELOCITY from all the structure conditions that intersect the element
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
