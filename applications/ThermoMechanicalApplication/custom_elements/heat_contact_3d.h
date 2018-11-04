//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Kazem Kamran
//                   Riccardo Rossi
//

#if !defined(KRATOS_HEAT_CONTACT_3D_CONDITION_H_INCLUDED )
#define  KRATOS_HEAT_CONTACT_3D_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "includes/serializer.h"

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

/// Short class definition.
/** Detail class definition.
*/
class HeatContact3D
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Monolithic2DNeumann
    KRATOS_CLASS_POINTER_DEFINITION(HeatContact3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HeatContact3D() {};
    HeatContact3D(IndexType NewId, GeometryType::Pointer pGeometry);
    HeatContact3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~HeatContact3D();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;
    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo) override;
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo) override;
//        void GetFirstDerivativesVector(Vector& values, int Step = 0);
//        void GetSecondDerivativesVector(Vector& values, int Step = 0);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
//      virtual String Info() const;

    /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


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
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
//       Monolithic2DNeumann() : Condition()
//       {
//       }

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }


    ///@}
    ///@name Private Operators
    ///@{
    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                      ProcessInfo& rCurrentProcessInfo,
                      bool CalculateStiffnessMatrixFlag,
                      bool CalculateResidualVectorFlag);


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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //Monolithic2DNeumann& operator=(const Monolithic2DNeumann& rOther);

    /// Copy constructor.
    //Monolithic2DNeumann(const Monolithic2DNeumann& rOther);


    ///@}

}; // Class Monolithic2DNeumann

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    Monolithic2DNeumann& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const Monolithic2DNeumann& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_HEAT_CONTACT_3D_CONDITION_H_INCLUDED   defined


