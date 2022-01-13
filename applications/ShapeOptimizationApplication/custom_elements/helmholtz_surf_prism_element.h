//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Reza Najian Asl
//

#if !defined(KRATOS_HELMHOLTZ_SURF_PRISM_ELEMENT_H_INCLUDED )
#define  KRATOS_HELMHOLTZ_SURF_PRISM_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "geometries/prism_3d_6.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "shape_optimization_application_variables.h"


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
class HelmholtzSurfPrismElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3>                                PointType;
    typedef Node<3>::Pointer                    PointPtrType;
    typedef Prism3D6<PointType>            PrismGeometryType;
    typedef Tetrahedra3D4<PointType>            TetrahedraGeometryType;    

    /// Counted pointer of HelmholtzSurfPrismElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(HelmholtzSurfPrismElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HelmholtzSurfPrismElement(IndexType NewId, GeometryType::Pointer pGeometry);
    HelmholtzSurfPrismElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~HelmholtzSurfPrismElement();


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
    HelmholtzSurfPrismElement() : Element()
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

    void CalculateSurfaceMassMatrix(MatrixType& rMassMatrix,const ProcessInfo& rCurrentProcessInfo) const;
    void CalculateSurfaceStiffnessMatrix(MatrixType& rStiffnessMatrix,const ProcessInfo& rCurrentProcessInfo) const;
    void CalculatePrismDN_DXMatrix(MatrixType& rDN_DX, const IntegrationMethod& rIntegrationMethod, const IndexType PointNumber, const ProcessInfo& rCurrentProcessInfo) const;
    void GetPrismShapeFunctionsValues(MatrixType& rNMatrix,const IntegrationMethod& rIntegrationMethod, const ProcessInfo& rCurrentProcessInfo) const;
    void CalculateTetrahedraDN_DXMatrix(MatrixType& rDN_DX, const IntegrationMethod& rIntegrationMethod, const IndexType PointNumber, const ProcessInfo& rCurrentProcessInfo) const;
    void GetTetrahedraShapeFunctionsValues(MatrixType& rNMatrix,const IntegrationMethod& rIntegrationMethod, const ProcessInfo& rCurrentProcessInfo) const;    
    void CalculateRotationMatrix(MatrixType& rRotMatrix,const ProcessInfo& rCurrentProcessInfo) const;
    void CalculateNormal(VectorType & r_n) const;
    void CalculateCMatrix(MatrixType& rCMatrix, const IntegrationMethod& rIntegrationMethod, const IndexType PointNumber) const;
    void CalculateBMatrix(MatrixType& rBMatrix, const MatrixType& rDN_DX_tMatrix, const IntegrationMethod& rIntegrationMethod, const IndexType PointNumber) const;

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
    //HelmholtzSurfPrismElement& operator=(const HelmholtzSurfPrismElement& rOther);

    /// Copy constructor.
    //HelmholtzSurfPrismElement(const HelmholtzSurfPrismElement& rOther);


    ///@}

}; // Class HelmholtzSurfPrismElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    HelmholtzSurfPrismElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const HelmholtzSurfPrismElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_HELMHOLTZ_SURF_PRISM_ELEMENT_H_INCLUDED  defined


