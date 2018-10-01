// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//                   Armin Geiser, https://github.com/armingeiser
//

// System includes
#if !defined(ADJOINT_SEMI_ANALYTIC_POINT_LOAD_CONDITION )
#define  ADJOINT_SEMI_ANALYTIC_POINT_LOAD_CONDITION

// System includes

// External includes

// Project includes
#include "adjoint_semi_analytic_base_condition.h"

namespace Kratos
{

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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  AdjointSemiAnalyticPointLoadCondition
    : public AdjointSemiAnalyticBaseCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AdjointSemiAnalyticPointLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( AdjointSemiAnalyticPointLoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    AdjointSemiAnalyticPointLoadCondition(): AdjointSemiAnalyticBaseCondition(){};

    AdjointSemiAnalyticPointLoadCondition(Condition::Pointer pPrimalCondition);

    /// Destructor.
    ~AdjointSemiAnalyticPointLoadCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo ) override;

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo ) override;

    void GetValuesVector(Vector& rValues, int Step = 0 ) override;

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, AdjointSemiAnalyticBaseCondition );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, AdjointSemiAnalyticBaseCondition );
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class AdjointSemiAnalyticPointLoadCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // ADJOINT_SEMI_ANALYTIC_POINT_LOAD_CONDITION  defined


