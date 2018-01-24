//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//  Co-authors:      Jordi Cotela
//

#ifndef KRATOS_SYMBOLIC_NAVIER_STOKES_H
#define KRATOS_SYMBOLIC_NAVIER_STOKES_H

#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"

#include "includes/cfd_variables.h"
#include "custom_elements/fluid_element.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

///@addtogroup FluidDynamicsApplication
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

template< class TElementData >
class SymbolicNavierStokes : public FluidElement<TElementData>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SymbolicNavierStokes
    KRATOS_CLASS_POINTER_DEFINITION(SymbolicNavierStokes);

    /// Node type (default is: Node<3>)
    typedef Node<3> NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    /// Type for shape function values container
    typedef typename FluidElement<TElementData>::ShapeFunctionsType ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef typename FluidElement<TElementData>::ShapeFunctionDerivativesType ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef typename FluidElement<TElementData>::ShapeFunctionDerivativesArrayType ShapeFunctionDerivativesArrayType;

    constexpr static unsigned int Dim = FluidElement<TElementData>::Dim;
    constexpr static unsigned int NumNodes = FluidElement<TElementData>::NumNodes;
    constexpr static unsigned int BlockSize = FluidElement<TElementData>::BlockSize;
    constexpr static unsigned int LocalSize = FluidElement<TElementData>::LocalSize;

    constexpr static unsigned int StrainSize = (Dim*3)-3;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    SymbolicNavierStokes(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    SymbolicNavierStokes(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    SymbolicNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    SymbolicNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~SymbolicNavierStokes();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new SymbolicNavierStokes element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            Properties::Pointer pProperties) const override;

    /// Create a new element of this type using given geometry
    /**
     * Returns a pointer to a new FluidElement element, created using given input
     * @param NewId: the ID of the new element
     * @param pGeom: a pointer to the geomerty to be used to create the element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            Properties::Pointer pProperties) const override;

    void Initialize() override;

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

    ///@name Protected member Variables
    ///@{

    /// Constitutive law pointer
    ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void AddTimeIntegratedSystem(
        TElementData& rData,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    void AddTimeIntegratedLHS(
        TElementData& rData,
        MatrixType& rLHS) override;

    void AddTimeIntegratedRHS(
        TElementData& rData,
        VectorType& rRHS) override;

    virtual void AddBoundaryIntegral(TElementData& rData,
        const Vector& rUnitNormal, MatrixType& rLHS, VectorType& rRHS);

    void ComputeGaussPointLHSContribution(
        TElementData& rData,
        MatrixType& rLHS);

    void ComputeGaussPointRHSContribution(
        TElementData& rData,
        VectorType& rRHS);

    void ComputeConstitutiveResponse(TElementData& rData);

    void ComputeStrain(TElementData& rData);

    /**
    * This functions sets the auxiliar matrix to compute the normal projection in Voigt notation.
    * 2D version.
    * @param rUnitNormal: reference to Gauss pt. unit normal vector
    * @param rVoigtNormProjMatrix: reference to the computed normal projection auxiliar matrix
    */
    void SetVoigtNormalProjectionMatrix(
        const array_1d<double, 3>& rUnitNormal,
        bounded_matrix<double, 2, 3>& rVoigtNormProjMatrix) const;
        
     
    /**
    * This functions sets the auxiliar matrix to compute the normal projection in Voigt notation.
    * 3D version
    * @param rUnitNormal: reference to Gauss pt. unit normal vector
    * @param rVoigtNormProjMatrix: reference to the computed normal projection auxiliar matrix
    */ 
    void SetVoigtNormalProjectionMatrix(
        const array_1d<double, 3>& rUnitNormal,
        bounded_matrix<double, 3, 6>& rVoigtNormProjMatrix) const;  

    /**
    * This functions sets the B strain matrix (pressure columns are set to zero).
    * 2D version.
    * @param rDN_DX: reference to the current Gauss pt. shape function gradients
    * @param rB_matrix: reference to the computed B strain matrix
    */
    void SetInterfaceStrainMatrix(
        const bounded_matrix<double, NumNodes, 2>& rDN_DX,
        bounded_matrix<double, 3, NumNodes*3>& rB_matrix) const;
        
    /**
    * This functions sets the B strain matrix (pressure columns are set to zero).
    * 3D version.
    * @param rDN_DX: reference to the current Gauss pt. shape function gradients
    * @param rB_matrix: reference to the computed B strain matrix
    */
    void SetInterfaceStrainMatrix(
        const bounded_matrix<double, NumNodes, 3>& rDN_DX,
        bounded_matrix<double, 6, NumNodes*4>& rB_matrix) const;

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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SymbolicNavierStokes& operator=(SymbolicNavierStokes const& rOther);

    /// Copy constructor.
    SymbolicNavierStokes(SymbolicNavierStokes const& rOther);

    ///@}


}; // Class SymbolicNavierStokes

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TElementData >
inline std::istream& operator >>(std::istream& rIStream,
                                 SymbolicNavierStokes<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const SymbolicNavierStokes<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_SYMBOLIC_NAVIER_STOKES_H
