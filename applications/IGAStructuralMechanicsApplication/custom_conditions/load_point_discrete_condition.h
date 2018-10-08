#if !defined(KRATOS_LOAD_POINT_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_LOAD_POINT_DISCRETE_CONDITION_H_INCLUDED



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
    class LoadPointDiscreteCondition
        : public PointBaseDiscreteCondition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer of LoadPointDiscreteCondition
        KRATOS_CLASS_POINTER_DEFINITION(LoadPointDiscreteCondition);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        // Constructor using an array of nodes
        LoadPointDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
            : PointBaseDiscreteCondition(NewId, pGeometry)
        {};
        // Constructor using an array of nodes with properties
        LoadPointDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
            : PointBaseDiscreteCondition(NewId, pGeometry, pProperties)
        {};

        LoadPointDiscreteCondition() : PointBaseDiscreteCondition()
        {};

        /// Destructor.
        virtual ~LoadPointDiscreteCondition() override
        {};

        Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
        {
            return Kratos::make_shared< LoadPointDiscreteCondition >(NewId, GetGeometry().Create(ThisNodes), pProperties);
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
            buffer << "\"LoadPointDiscreteCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"LoadPointDiscreteCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const {
            pGetGeometry()->PrintData(rOStream);
        }

        /**
        * this is called at the end of each solution step
        */
        void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

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
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        }
        ///@}

    }; // Class LoadPointDiscreteCondition


}  // namespace Kratos.

#endif // KRATOS_LOAD_SURFACE_DISCRETE_CONDITION_H_INCLUDED  defined 