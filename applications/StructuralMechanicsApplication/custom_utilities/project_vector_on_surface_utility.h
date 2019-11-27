//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus Sautter
//                   Philipp Bucher
//
//

#if !defined(KRATOS_PROJECT_VECTOR_ON_SURFACE_UTILITY_H )
#define  KRATOS_PROJECT_VECTOR_ON_SURFACE_UTILITY_H

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ProjectVectorOnSurfaceUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ProjectVectorOnSurfaceUtility
    KRATOS_CLASS_POINTER_DEFINITION(ProjectVectorOnSurfaceUtility);

    typedef array_1d<double, 3> Vector3;
    typedef Variable< array_1d< double, 3> > ArrayVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    //ProjectVectorOnSurfaceUtility();

    /// Destructor.
    virtual ~ProjectVectorOnSurfaceUtility() = default;


    ///@}
    ///@name Operations
    ///@{

    void Execute(ModelPart& rModelPart, const Parameters MethodParameters);


    ///@}
    ///@name Input and output
    ///@{

    ///@}

    void CheckAndReadVectors(Parameters ThisParameters, const std::string KeyName, Vector3 &rVector);

    void PlanarProjection(ModelPart& rModelPart, const Parameters& MethodParameters);
    void RadialProjection(ModelPart& rModelPart, const Parameters& MethodParameters);
    void SphericalProjection(ModelPart& rModelPart, const Parameters& MethodParameters);

private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    ///@}

}; // Class ProjectVectorOnSurfaceUtility

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ASSIGN_MATERIAL_ORIENTATION_UTILITY_H  defined
