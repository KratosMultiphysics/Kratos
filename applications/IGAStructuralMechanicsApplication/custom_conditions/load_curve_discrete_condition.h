#if !defined(KRATOS_LOAD_CURVE_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_LOAD_CURVE_DISCRETE_CONDITION_H_INCLUDED



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
    class LoadCurveDiscreteCondition
        : public CurveBaseDiscreteCondition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer of LoadCurveDiscreteCondition
        KRATOS_CLASS_POINTER_DEFINITION(LoadCurveDiscreteCondition);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        // Constructor using an array of nodes
        LoadCurveDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
            : CurveBaseDiscreteCondition(NewId, pGeometry)
        {};
        // Constructor using an array of nodes with properties
        LoadCurveDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
            : CurveBaseDiscreteCondition(NewId, pGeometry, pProperties)
        {};

        LoadCurveDiscreteCondition() : CurveBaseDiscreteCondition()
        {};

        /// Destructor.
        virtual ~LoadCurveDiscreteCondition() override
        {};

        Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
        {
            return Kratos::make_shared< LoadCurveDiscreteCondition >(NewId, GetGeometry().Create(ThisNodes), pProperties);
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

        void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

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

        ///@}

    protected:

    private:

        void get_var_of_small_rotation_wrt_disp_global(
            Matrix& phi_r);
        ///@name Static Member Variables
        ///@{
        Vector m_g10;
        Vector m_g20;
        Vector m_g30;
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

    }; // Class LoadCurveDiscreteCondition


}  // namespace Kratos.

#endif // KRATOS_LOAD_CURVE_DISCRETE_CONDITION_H_INCLUDED  defined 