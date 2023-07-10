//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

#if !defined(KRATOS_HELMHOLTZ_SURF_SHAPE_ELEMENT_H_INCLUDED )
#define  KRATOS_HELMHOLTZ_SURF_SHAPE_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/pyramid_3d_5.h"
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
/** Detail class definition.
*/
class HelmholtzSurfShapeElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node                             PointType;
    typedef Node::Pointer                    PointPtrType;
    typedef Geometry<PointType>                 GeometryType;
    typedef Pyramid3D5<PointType>               PyramidGeometryType;
    typedef Tetrahedra3D4<PointType>            TetrahedraGeometryType;    

    /// Counted pointer of HelmholtzSurfShapeElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(HelmholtzSurfShapeElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HelmholtzSurfShapeElement(IndexType NewId, GeometryType::Pointer pGeometry);
    HelmholtzSurfShapeElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~HelmholtzSurfShapeElement();


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
    HelmholtzSurfShapeElement() : Element()
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
    void GetPseudoBulkSurfaceShapeFunctionsValues(MatrixType& rNMatrix,const IntegrationMethod& rIntegrationMethod, const ProcessInfo& rCurrentProcessInfo) const;
    void CalculatePseudoBulkSurfaceDN_DXMatrix(MatrixType& rDN_DX, const IntegrationMethod& rIntegrationMethod, const IndexType PointNumber, const ProcessInfo& rCurrentProcessInfo) const;   
    void CalculateAvgSurfUnitNormal(VectorType & rNormal) const;
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
    //HelmholtzSurfShapeElement& operator=(const HelmholtzSurfShapeElement& rOther);

    /// Copy constructor.
    //HelmholtzSurfShapeElement(const HelmholtzSurfShapeElement& rOther);


    ///@}

}; // Class HelmholtzSurfShapeElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    HelmholtzSurfShapeElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const HelmholtzSurfShapeElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_HELMHOLTZ_SURF_SHAPE_ELEMENT_H_INCLUDED  defined


