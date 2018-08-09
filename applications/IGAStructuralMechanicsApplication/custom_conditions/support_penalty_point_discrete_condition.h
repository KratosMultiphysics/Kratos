#if !defined(KRATOS_SUPPORT_PENALTY_POINT_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_SUPPORT_PENALTY_POINT_DISCRETE_CONDITION_H_INCLUDED



// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "iga_structural_mechanics_application_variables.h"
#include "custom_conditions/point_base_discrete_condition.h"


namespace Kratos
{
    class SupportPenaltyPointDiscreteCondition
        : public PointBaseDiscreteCondition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer of SupportPenaltyPointDiscreteCondition
        KRATOS_CLASS_POINTER_DEFINITION(SupportPenaltyPointDiscreteCondition);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        // Constructor using an array of nodes
        SupportPenaltyPointDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
            : PointBaseDiscreteCondition(NewId, pGeometry)
        {};
        // Constructor using an array of nodes with properties
        SupportPenaltyPointDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
            : PointBaseDiscreteCondition(NewId, pGeometry, pProperties)
        {};

        SupportPenaltyPointDiscreteCondition() : PointBaseDiscreteCondition()
        {};

        /// Destructor.
        virtual ~SupportPenaltyPointDiscreteCondition() override
        {};

        Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
        {
            return Kratos::make_shared< SupportPenaltyPointDiscreteCondition >(NewId, GetGeometry().Create(ThisNodes), pProperties);
        };


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

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
            buffer << "\"SupportPenaltyPointDiscreteCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"SupportPenaltyPointDiscreteCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const {
            pGetGeometry()->PrintData(rOStream);
        }

        ///@}

    protected:

    private:

        ///@name Static Member Variables
        ///@{
        
        ///@}
        ///@name Member Variables
        ///@{

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PointBaseDiscreteCondition);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PointBaseDiscreteCondition);
        }
        ///@}

    }; // Class SupportPenaltyPointDiscreteCondition


}  // namespace Kratos.

#endif // KRATOS_SUPPORT_PENALTY_POINT_DISCRETE_CONDITION_H_INCLUDED  defined 