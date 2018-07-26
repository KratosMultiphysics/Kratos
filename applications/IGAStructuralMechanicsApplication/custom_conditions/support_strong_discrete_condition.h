#if !defined(KRATOS_SUPPORT_STRONG_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_SUPPORT_STRONG_DISCRETE_CONDITION_H_INCLUDED



// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "iga_structural_mechanics_application_variables.h"
#include "custom_conditions/curve_base_discrete_condition.h"


namespace Kratos
{
    class SupportStrongDiscreteCondition
        : public BaseDiscreteCondition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer of SupportStrongDiscreteCondition
        KRATOS_CLASS_POINTER_DEFINITION(SupportStrongDiscreteCondition);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        // Constructor using an array of nodes
        SupportStrongDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
            : BaseDiscreteCondition(NewId, pGeometry)
        {};
        // Constructor using an array of nodes with properties
        SupportStrongDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
            : BaseDiscreteCondition(NewId, pGeometry, pProperties)
        {};

        SupportStrongDiscreteCondition() : BaseDiscreteCondition()
        {};

        /// Destructor.
        virtual ~SupportStrongDiscreteCondition() override
        {};

        Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
        {
            return Kratos::make_shared< SupportStrongDiscreteCondition >(NewId, GetGeometry().Create(ThisNodes), pProperties);
        };


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /**
        * this is called in the beginning of each solution step
        */
        void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);


        /// Turn back information as a string.
        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "\"SupportStrongDiscreteCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"SupportStrongDiscreteCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const {
            pGetGeometry()->PrintData(rOStream);
        }

        ///@}

    protected:

    private:

        void CalculateRotation(const Matrix& ShapeFunctionDerivatives,
            Vector& Phi_r, array_1d<double, 2>& Phi);

        void GetBaseVectorsSurface(
            const Matrix& DN_De,
            Vector& g1,
            Vector& g2,
            Vector& g3);

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
        ) override;

        ///@name Static Member Variables
        ///@{
        
        ///@}
        ///@name Member Variables
        ///@{

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseDiscreteCondition);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseDiscreteCondition);
        }
        ///@}

    }; // Class SupportStrongDiscreteCondition


}  // namespace Kratos.

#endif // KRATOS_SUPPORT_STRONG_DISCRETE_CONDITION_H_INCLUDED  defined 