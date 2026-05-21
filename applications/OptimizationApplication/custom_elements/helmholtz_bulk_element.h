//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//

#if !defined(KRATOS_HELMHOLTZ_BULK_ELEMENT_H_INCLUDED )
#define  KRATOS_HELMHOLTZ_BULK_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "optimization_application_variables.h"


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
/**
 * @class HelmholtzBulkElement
 * @ingroup OptimizationApplication
 * @brief Helmholtz filtering element for 3D geometries.
 * @details Implements Sobolev/Helmholtz filter for optimization and regularization. This works for arbitrary geometries in 3D
 * @author Reza Najian Asl
 */
class KRATOS_API(OPTIMIZATION_APPLICATION) HelmholtzBulkElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of HelmholtzBulkElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(HelmholtzBulkElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HelmholtzBulkElement(IndexType NewId, GeometryType::Pointer pGeometry);
    HelmholtzBulkElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~HelmholtzBulkElement();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    void GetValuesVector(VectorType &rValues, int Step = 0) const override;

    void Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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

    // Protected default constructor necessary for serialization
    HelmholtzBulkElement() : Element()
    {
    }

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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
    void CalculateBulkMassMatrix(MatrixType& rMassMatrix,const ProcessInfo& rCurrentProcessInfo) const;
    void CalculateBulkStiffnessMatrix(MatrixType& rStiffnessMatrix,const ProcessInfo& rCurrentProcessInfo) const;

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
    //HelmholtzBulkElement& operator=(const HelmholtzBulkElement& rOther);

    /// Copy constructor.
    //HelmholtzBulkElement(const HelmholtzBulkElement& rOther);


    ///@}

}; // Class HelmholtzBulkElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    HelmholtzBulkElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const HelmholtzBulkElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_HELMHOLTZ_BULK_ELEMENT_H_INCLUDED  defined


