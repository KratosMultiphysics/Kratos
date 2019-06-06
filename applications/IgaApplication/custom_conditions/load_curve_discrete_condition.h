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

#if !defined(KRATOS_LOAD_CURVE_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_LOAD_CURVE_DISCRETE_CONDITION_H_INCLUDED


// System includes

// External includes

// Project includes
#include "iga_application_variables.h"

#include "custom_conditions/base_discrete_condition.h"

#include "custom_utilities/geometry_utilities/iga_curve_on_surface_utilities.h"
#include "custom_utilities/geometry_utilities/iga_surface_utilities.h"

#include "includes/define.h"
#include "includes/condition.h"

namespace Kratos
{
    class LoadCurveDiscreteCondition
        : public BaseDiscreteCondition
    {
    public:

        /// Counted pointer of LoadCurveDiscreteCondition
        KRATOS_CLASS_POINTER_DEFINITION(LoadCurveDiscreteCondition);

        /// Default constructor
        LoadCurveDiscreteCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : BaseDiscreteCondition(NewId, pGeometry)
        {};

        LoadCurveDiscreteCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties)
            : BaseDiscreteCondition(NewId, pGeometry, pProperties)
        {};

        LoadCurveDiscreteCondition() : BaseDiscreteCondition()
        {};

        /// Destructor
        virtual ~LoadCurveDiscreteCondition() override
        {};

        Condition::Pointer Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties
        ) const override
        {
            return Kratos::make_intrusive<LoadCurveDiscreteCondition>(
                NewId, 
                GetGeometry().Create(ThisNodes), 
                pProperties);
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

        /// Turn back information as a string.
        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "\"LoadCurveDiscreteCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"LoadCurveDiscreteCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const {
            pGetGeometry()->PrintData(rOStream);
        }

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        }
        ///@}

    }; // Class LoadCurveDiscreteCondition

}  // namespace Kratos.

#endif // KRATOS_LOAD_CURVE_DISCRETE_CONDITION_H_INCLUDED  defined