#if !defined(KRATOS_SUPPORT_PENALTY_CURVE_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_SUPPORT_PENALTY_CURVE_DISCRETE_CONDITION_H_INCLUDED



// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "iga_application_variables.h"
#include "custom_conditions/curve_base_discrete_condition.h"


namespace Kratos
{
    class SupportPenaltyCurveDiscreteCondition
        : public CurveBaseDiscreteCondition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer of SupportPenaltyCurveDiscreteCondition
        KRATOS_CLASS_POINTER_DEFINITION(SupportPenaltyCurveDiscreteCondition);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        // Constructor using an array of nodes
        SupportPenaltyCurveDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
            : CurveBaseDiscreteCondition(NewId, pGeometry)
        {};
        // Constructor using an array of nodes with properties
        SupportPenaltyCurveDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
            : CurveBaseDiscreteCondition(NewId, pGeometry, pProperties)
        {};

        SupportPenaltyCurveDiscreteCondition() : CurveBaseDiscreteCondition()
        {};

        /// Destructor.
        virtual ~SupportPenaltyCurveDiscreteCondition() override
        {};

        Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
        {
            return Kratos::make_shared< SupportPenaltyCurveDiscreteCondition >(NewId, GetGeometry().Create(ThisNodes), pProperties);
        };


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{
        /**
        * Called to initialize the element.
        * Must be called before any calculation is done
        */
        void Initialize() override;

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

        ///@}

    protected:

    private:

        Vector m_g10;
        Vector m_g20;
        Vector m_g30;

        void CalculateRotation(const Matrix& ShapeFunctionDerivatives,
            Vector& Phi_r, Matrix& Phi_rs, array_1d<double, 2>& Phi);
        ///@name Static Member Variables
        ///@{
        
        ///@}
        ///@name Member Variables
        ///@{

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, CurveBaseDiscreteCondition);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, CurveBaseDiscreteCondition);
        }
        ///@}

    }; // Class SupportPenaltyCurveDiscreteCondition


}  // namespace Kratos.

#endif // KRATOS_SUPPORT_PENALTY_CURVE_DISCRETE_CONDITION_H_INCLUDED  defined 