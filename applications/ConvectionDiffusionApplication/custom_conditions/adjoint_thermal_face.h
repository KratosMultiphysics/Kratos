// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Jordi Cotela
//

#if !defined(KRATOS_ADJOINT_THERMAL_FACE_H_INCLUDED )
#define  KRATOS_ADJOINT_THERMAL_FACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/condition.h"

namespace Kratos
{

///@addtogroup ConvectionDiffusionApplication
///@{

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

/// Heat flux Neumann condition for the ajdoint thermal diffusion problem.
template< class PrimalCondition >
class AdjointThermalFace: public PrimalCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AdjointThermalFace
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointThermalFace);

    using IndexType = typename PrimalCondition::IndexType;
    using GeometryType = typename PrimalCondition::GeometryType;
    using NodesArrayType = typename PrimalCondition::NodesArrayType;
    using MatrixType = typename PrimalCondition::MatrixType;
    using VectorType = typename PrimalCondition::VectorType;
    using EquationIdVectorType = typename PrimalCondition::EquationIdVectorType;
    using DofsVectorType = typename PrimalCondition::DofsVectorType;

    constexpr static double StefanBoltzmann = PrimalCondition::StefanBoltzmann;

    ///@}
    ///@name Life Cycle
    ///@{

    AdjointThermalFace(IndexType NewId, typename GeometryType::Pointer pGeometry);

    AdjointThermalFace(IndexType NewId, typename GeometryType::Pointer pGeometry,  Properties::Pointer pProperties);

    /// Destructor.
    ~AdjointThermalFace() override;

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  Properties::Pointer pProperties) const override;

    Condition::Pointer Create(IndexType NewId, typename GeometryType::Pointer pGeom,  Properties::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step = 0) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    AdjointThermalFace() : PrimalCondition()
    {
    }

    void save(Serializer& rSerializer) const override
    {
        using BaseType = PrimalCondition;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    void load(Serializer& rSerializer) override
    {
        using BaseType = PrimalCondition;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // Note: the Jacobian should come from the geometry, but non-square ones are often wrongly implemented at the moment :(
    MatrixType GetJacobian(GeometryData::IntegrationMethod QuadratureOrder, unsigned int IntegrationPointIndex) const;

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    AdjointThermalFace& operator=(const AdjointThermalFace& rOther) = delete;

    /// Copy constructor.
    AdjointThermalFace(const AdjointThermalFace& rOther) = delete;

    ///@}

}; // Class AdjointThermalFace

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class PrimalCondition >
inline std::istream& operator >> (std::istream& rIStream, AdjointThermalFace<PrimalCondition>& rThis)
{
    return rIStream;
}

/// output stream function
template< class PrimalCondition >
inline std::ostream& operator << (std::ostream& rOStream, const AdjointThermalFace<PrimalCondition>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@}

}  // namespace Kratos.

#endif // KRATOS_ADJOINT_THERMAL_FACE_H_INCLUDED  defined


