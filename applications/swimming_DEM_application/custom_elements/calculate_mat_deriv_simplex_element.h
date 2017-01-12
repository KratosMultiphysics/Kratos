//
//   Project Name:        Kratos
//   Last Modified by:    $Author: gcasas $
//   Date:                $Date: 2016-03-12
//

#if !defined(KRATOS_COMPUTE_MATERIAL_DERIVATIVE_SIMPLEX_ELEMENT_H_INCLUDED )
#define  KRATOS_COMPUTE_MATERIAL_DERIVATIVE_SIMPLEX_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "includes/variables.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

///@addtogroup SwimmingDEMApplication
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

/// A post-processing element to recover the Laplacian from the velocity solution.
/**
 */
template< unsigned int TDim,
          unsigned int TNumNodes = TDim + 1 >
class ComputeMaterialDerivativeSimplex : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeMaterialDerivativeSimplex
    KRATOS_CLASS_POINTER_DEFINITION(ComputeMaterialDerivativeSimplex);

    /// Node type (default is: Node<3>)
    typedef Node <3> NodeType;

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

    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    ComputeMaterialDerivativeSimplex(IndexType NewId = 0) :
        Element(NewId)
    {}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    ComputeMaterialDerivativeSimplex(IndexType NewId, const NodesArrayType& ThisNodes) :
        Element(NewId, ThisNodes)
    {}

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    ComputeMaterialDerivativeSimplex(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
    {}

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    ComputeMaterialDerivativeSimplex(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ComputeMaterialDerivativeSimplex()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new ComputeMaterialDerivativeSimplex element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new ComputeMaterialDerivativeSimplex(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Calculate the element's local contribution to the system for the current step.
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo)
    {
        //KRATOS_WATCH(this->Id())
        const unsigned int NumNodes(TDim+1), LocalSize(TDim * NumNodes);

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);

        for (unsigned int i=0; i<LocalSize; ++i){
            for (unsigned int j=0; j<LocalSize; ++j){
                rLeftHandSideMatrix(i, j) = 0.0;
            }
            rRightHandSideVector(i) = 0.0;
        }

        boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim > DN_DX;
        array_1d<double, TDim+1 > N;
        double Area;

        CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
        //getting data for the given geometry
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

        CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
    }


    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    virtual void EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo)
    {

        const unsigned int NumNodes(TDim+1), LocalSize(TDim * NumNodes);
        unsigned int LocalIndex = 0;
        unsigned int local_position = this->GetGeometry()[0].GetDofPosition(MATERIAL_ACCELERATION_X);

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(MATERIAL_ACCELERATION_X, local_position).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(MATERIAL_ACCELERATION_Y, local_position + 1).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(MATERIAL_ACCELERATION_Z, local_position + 2).EquationId();
        }
//Z
    }

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void GetDofList(DofsVectorType& rElementalDofList,
                            ProcessInfo& rCurrentProcessInfo)
    {

        const unsigned int NumNodes(TDim+1), LocalSize(TDim * NumNodes);

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(MATERIAL_ACCELERATION_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(MATERIAL_ACCELERATION_Y);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(MATERIAL_ACCELERATION_Z);
        }
//Z

    }



    /// Obtain an array_1d<double,3> elemental variable, evaluated on gauss points.
    /**
     * @param rVariable Kratos vector variable to get
     * @param Output Will be filled with the values of the variable on integrartion points
     * @param rCurrentProcessInfo Process info instance
     */
//    virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
//            std::vector<array_1d<double, 3 > >& rValues,
//            const ProcessInfo& rCurrentProcessInfo)
//    {
//
//    }


    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Elemental Data
    ///@{

    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    virtual int Check(const ProcessInfo& rCurrentProcessInfo)
{
        KRATOS_TRY

        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if(ErrorCode != 0) return ErrorCode;

        if(this->GetGeometry().size() != TDim+1)
            KRATOS_THROW_ERROR(std::invalid_argument,"wrong number of nodes for element",this->Id());

        // Check that all required variables have been registered
//G
//        if(DISTANCE.Key() == 0)

//            //KRATOS_THROW_ERROR(std::invalid_argument,"DISTANCE Key is 0. Check if the application was correctly registered.","");


//        // Checks on nodes

//        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
//        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
//        {
//            if(this->GetGeometry()[i].SolutionStepsDataHas(DISTANCE) == false)
//                KRATOS_THROW_ERROR(std::invalid_argument,"missing DISTANCE variable on solution step data for node ",this->GetGeometry()[i].Id());
//        }

        if (MATERIAL_ACCELERATION.Key() == 0){
            KRATOS_THROW_ERROR(std::invalid_argument,"MATERIAL_ACCELERATION Key is 0. Check if the application was correctly registered.","");
        }

        // Checks on nodes

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(MATERIAL_ACCELERATION) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing MATERIAL_ACCELERATION variable on solution step data for node ",this->GetGeometry()[i].Id());
        }
//Z
        return 0;

        KRATOS_CATCH("");
    }


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "ComputeMaterialDerivativeSimplex #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ComputeMaterialDerivativeSimplex" << TDim << "D";
    }

//        /// Print object's data.
//        virtual void PrintData(std::ostream& rOStream) const;

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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Private Operators
    ///@{

    /// Computes local contributions to the mass matrix
    /**
     * Provides the local contributions to the mass matrix, which is defined here
     * as the matrix associated to velocity derivatives. Note that the mass
     * matrix implemented here is lumped.
     * @param rMassMatrix Will be filled with the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int LocalSize = TDim * TNumNodes;

        // Resize and set to zero
        if (rMassMatrix.size1() != LocalSize)
            rMassMatrix.resize(LocalSize, LocalSize, false);

        rMassMatrix = ZeroMatrix(LocalSize, LocalSize);

        // Get the element's geometric parameters
        double Area;
        array_1d<double, TNumNodes> N;
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);


        // Add 'classical' mass matrix (lumped)
        if (rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] == 1){
            double Coeff = Area / TNumNodes; //Optimize!
            this->CalculateLumpedMassMatrix(rMassMatrix, Coeff);
        }

        else {
            // Add 'consistent' mass matrix
            MatrixType NContainer;
            ShapeFunctionDerivativesArrayType DN_DXContainer;
            VectorType GaussWeights;
            this->CalculateWeights(DN_DXContainer, NContainer, GaussWeights);
            const SizeType NumGauss = NContainer.size1();

            for (SizeType g = 0; g < NumGauss; g++){
                const double GaussWeight = GaussWeights[g];
                const ShapeFunctionsType& Ng = row(NContainer, g);
                this->AddConsistentMassMatrixContribution(rMassMatrix, Ng, GaussWeight);
            }
        }
    }

    void CalculateLumpedMassMatrix(MatrixType& rLHSMatrix,
                                   const double Mass)
    {
        unsigned int DofIndex = 0;
        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
        {
            for (unsigned int d = 0; d < TDim; ++d)
            {
                rLHSMatrix(DofIndex, DofIndex) += Mass;
                ++DofIndex;
            }
//G
//            ++DofIndex; // Skip pressure Dof
//Z
        }
    }

    virtual void CalculateRHS(VectorType& F, ProcessInfo& rCurrentProcessInfo)
    {
        // Get the element's geometric parameters
        double Area;
        array_1d<double, TNumNodes> N;
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        MatrixType NContainer;
        ShapeFunctionDerivativesArrayType DN_DXContainer;
        VectorType GaussWeights;
        this->CalculateWeights(DN_DXContainer, NContainer, GaussWeights);
        const SizeType NumGauss = NContainer.size1();

        for (SizeType g = 0; g < NumGauss; g++){
            const double GaussWeight = GaussWeights[g];
            const ShapeFunctionsType& Ng = row(NContainer, g);
            this->AddRHSMatAcceleration(F, Ng, DN_DX, GaussWeight);
        }

    }

    void AddConsistentMassMatrixContribution(MatrixType& rLHSMatrix,
            const array_1d<double,TNumNodes>& rShapeFunc,
            const double Weight)
    {
//G
        const unsigned int BlockSize = TDim;

        double Coef = Weight;
//Z
        unsigned int FirstRow(0), FirstCol(0);
        double K; // Temporary results

        // Note: Dof order is (vx,vy,[vz,]p) for each node
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            // Loop over columns
            for (unsigned int j = 0; j < TNumNodes; ++j)
            {
                K = Coef * rShapeFunc[i] * rShapeFunc[j];

                for (unsigned int d = 0; d < TDim; ++d) // iterate over dimensions for velocity Dofs in this node combination
                {
                    rLHSMatrix(FirstRow + d, FirstCol + d) += K;
                }
                // Update column index
                FirstCol += BlockSize;
            }
            // Update matrix indices
            FirstRow += BlockSize;
            FirstCol = 0;
        }
    }

    void CalculateWeights(ShapeFunctionDerivativesArrayType& rDN_DX, Matrix& rNContainer, Vector& rGaussWeights);
    void EvaluateInPoint(array_1d< double, 3 > & rResult,
                         const Variable< array_1d< double, 3 > >& rVariable,
                         const array_1d< double, TNumNodes >& rShapeFunc);

    void AddRHSMatAcceleration(VectorType& F,
                         const array_1d<double,TNumNodes>& rShapeFunc,
                         const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim>& rShapeDeriv,
                         const double Weight);
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
    ComputeMaterialDerivativeSimplex & operator=(ComputeMaterialDerivativeSimplex const& rOther);

    /// Copy constructor.
    ComputeMaterialDerivativeSimplex(ComputeMaterialDerivativeSimplex const& rOther);

    ///@}

}; // Class ComputeMaterialDerivativeSimplex

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim >
inline std::istream& operator >>(std::istream& rIStream,
                                 ComputeMaterialDerivativeSimplex<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const ComputeMaterialDerivativeSimplex<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_COMPUTE_MATERIAL_DERIVATIVE_SIMPLEX_ELEMENT_H_INCLUDED  defined
