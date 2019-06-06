//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//


#if !defined(KRATOS_SUPPORT_PENALTY_CURVE_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_SUPPORT_PENALTY_CURVE_DISCRETE_CONDITION_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "custom_utilities/iga_flags.h"

#include "iga_application_variables.h"
#include "custom_conditions/base_discrete_condition.h"

#include "custom_utilities/geometry_utilities/iga_surface_utilities.h"
#include "custom_utilities/geometry_utilities/iga_curve_on_surface_utilities.h"

namespace Kratos
{
    class SupportPenaltyCurveDiscreteCondition
        : public BaseDiscreteCondition
    {
    public:
        /// Counted pointer of SupportPenaltyCurveDiscreteCondition
        KRATOS_CLASS_POINTER_DEFINITION(SupportPenaltyCurveDiscreteCondition);

        /// Default constructor.
        SupportPenaltyCurveDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
            : BaseDiscreteCondition(NewId, pGeometry)
        {};

        // with properties
        SupportPenaltyCurveDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
            : BaseDiscreteCondition(NewId, pGeometry, pProperties)
        {};

        SupportPenaltyCurveDiscreteCondition() : BaseDiscreteCondition()
        {};

        /// Destructor.
        virtual ~SupportPenaltyCurveDiscreteCondition() override
        {};

        /**
        * @brief Creates a new element
        * @param NewId The Id of the new created element
        * @param pGeom The pointer to the geometry of the element
        * @param pProperties The pointer to property
        * @return The pointer to the created element
        */
        Condition::Pointer Create(
            IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties
        ) const override
        {
            return Kratos::make_intrusive<SupportPenaltyCurveDiscreteCondition>(
                NewId, pGeom, pProperties);
        };

        Condition::Pointer Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties
        ) const override
        {
            return Kratos::make_intrusive< SupportPenaltyCurveDiscreteCondition >(
                NewId, GetGeometry().Create(ThisNodes), pProperties);
        };

        /**
        * @brief Sets on rResult the ID's of the element degrees of freedom
        * @param rResult The vector containing the equation id
        * @param rCurrentProcessInfo The current process info instance
        */
        void EquationIdVector(
            EquationIdVectorType& rResult,
            ProcessInfo& rCurrentProcessInfo
        ) override;

        /**
        * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
        * @param rElementalDofList The vector containing the dof of the element
        * @param rCurrentProcessInfo The current process info instance
        */
        void GetDofList(
            DofsVectorType& rElementalDofList,
            ProcessInfo& rCurrentProcessInfo
        ) override;

        /**
        * This functions calculates both the RHS and the LHS
        * @param rLeftHandSideMatrix: The LHS
        * @param rRightHandSideVector: The RHS
        * @param rCurrentProcessInfo: The current process info instance
        * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
        * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
        */
        void CalculateAll(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo,
            const bool CalculateStiffnessMatrixFlag,
            const bool CalculateResidualVectorFlag
        );

        /**
        * @brief Called to initialize the element. Sets the base vectors of the 
        * integration point on the underlying surface.
        */
        void Initialize() override;

        /// Turn back information as a string.
        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "\"SupportPenaltyCurveDiscreteCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"SupportPenaltyCurveDiscreteCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const {
            pGetGeometry()->PrintData(rOStream);
        }

    private:

        array_1d<double, 3> mG10;
        array_1d<double, 3> mG20;
        array_1d<double, 3> mG30;

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseDiscreteCondition);
            rSerializer.save("mG10", mG10);
            rSerializer.save("mG20", mG20);
            rSerializer.save("mG30", mG30);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseDiscreteCondition);
            rSerializer.load("mG10", mG10);
            rSerializer.load("mG20", mG20);
            rSerializer.load("mG30", mG30);
        }

    }; // Class SupportPenaltyCurveDiscreteCondition

}  // namespace Kratos.

#endif // KRATOS_SUPPORT_PENALTY_CURVE_DISCRETE_CONDITION_H_INCLUDED  defined 