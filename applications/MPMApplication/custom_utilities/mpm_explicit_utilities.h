//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

#if !defined(KRATOS_MPM_EXPLICIT_UTILITIES)
#define KRATOS_MPM_EXPLICIT_UTILITIES

// Project includes
#include "includes/model_part.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "mpm_application_variables.h"
#include "custom_utilities/mpm_math_utilities.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{
    /**
     * @namespace MPMExplicitUtilities
     * @ingroup MPMApplication
     * @brief This namespace includes several utilities necessaries for the computation of the explicit integration
     * @author Peter Wilson
     */
    namespace MPMExplicitUtilities
    {
        /// The size type definition
        typedef std::size_t SizeType;

        /// The index type definition
        typedef std::size_t IndexType;

        /// The arrays of elements and nodes
        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ModelPart::NodesContainerType NodesArrayType;
        typedef Node NodeType;
        typedef Geometry<NodeType> GeometryType;

        void KRATOS_API(MPM_APPLICATION) CalculateAndAddExplicitInternalForce(const ProcessInfo& rProcessInfo,
            Element& rElement, const Vector& rMPStress, const double rMPVolume,
            const SizeType StrainSize, Vector& rRightHandSideVector);

        void KRATOS_API(MPM_APPLICATION) UpdateGaussPointExplicit(const ProcessInfo& rCurrentProcessInfo,
            Element& rElement);

        void KRATOS_API(MPM_APPLICATION) CalculateMUSLGridVelocity(const ProcessInfo& rCurrentProcessInfo,
            Element& rElement);

        void KRATOS_API(MPM_APPLICATION) CalculateExplicitKinematics(const ProcessInfo& rCurrentProcessInfo,
            Element& rElement, Vector& rMPStrain, Matrix& rDeformationGradient,
            const SizeType StrainSize);

        inline void GetCartesianDerivatives(std::vector<Matrix>& rDN_DXVec, GeometryType& rGeom);
    }; // namespace ExplicitIntegrationUtilities
}  // namespace Kratos
#endif /* KRATOS_MPM_EXPLICIT_UTILITIES defined */