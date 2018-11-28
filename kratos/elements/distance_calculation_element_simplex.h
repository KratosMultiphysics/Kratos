//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//

#if !defined(KRATOS_DISTANCE_CALCULATION_ELEMENT_H_INCLUDED )
#define  KRATOS_DISTANCE_CALCULATION_ELEMENT_H_INCLUDED

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

/// A stabilized element for the incompressible Navier-Stokes equations.
/**
 */
template< unsigned int TDim >
class DistanceCalculationElementSimplex : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceCalculationElementSimplex
    KRATOS_CLASS_POINTER_DEFINITION(DistanceCalculationElementSimplex);

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
    DistanceCalculationElementSimplex(IndexType NewId = 0) :
        Element(NewId)
    {}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    DistanceCalculationElementSimplex(IndexType NewId, const NodesArrayType& ThisNodes) :
        Element(NewId, ThisNodes)
    {}

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    DistanceCalculationElementSimplex(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
    {}

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    DistanceCalculationElementSimplex(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~DistanceCalculationElementSimplex() override
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new DistanceCalculationElementSimplex element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Element::Pointer(new DistanceCalculationElementSimplex(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Calculate the element's local contribution to the system for the current step.
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int number_of_points = TDim+1;

        BoundedMatrix<double, TDim+1, TDim > DN_DX;
        array_1d<double, TDim+1 > N;

        if (rLeftHandSideMatrix.size1() != number_of_points)
            rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);

        if (rRightHandSideVector.size() != number_of_points)
            rRightHandSideVector.resize(number_of_points, false);

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
        
        //get distances at the nodes
        array_1d<double, TDim+1 > distances;
        for(unsigned int i=0; i<number_of_points; i++)
        {
            distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
//			if (distances[i] >= 0.0 && distances[i] < 1e-30) distances[i] = 1e-30;
        }

       const double dgauss = inner_prod(N,distances);

        const unsigned int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    
        if(step == 1) //solve a poisson problem with a positive/negative heat source depending on the sign of the existing distance function
        {
            //compute distance on gauss point

            this->SetValue(DISTANCE,dgauss); //saving the distance, to see if it changed sign between iterations
            //compute LHS
            noalias(rLeftHandSideMatrix) = Area*prod(DN_DX,trans(DN_DX));
            
            //compute RHS
            double source = 1.0;
            if(dgauss < 0.0) source=-1.0;
            noalias(rRightHandSideVector) = source*Area*N;            
            
            noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,distances); 
            
            
            //impose that the normal gradient is 1 on outer faces
            unsigned int nboundary = 0;
            for(unsigned int i=0; i<TDim+1; i++)
                if(GetGeometry()[i].Is(BOUNDARY)) nboundary +=1;
            
            if(nboundary == TDim)
            {
	      array_1d<double,TDim> DN_out(TDim, 0.0);
                for(unsigned int i=0; i<TDim+1; i++)
                    if(GetGeometry()[i].IsNot(BOUNDARY))
                    {
                        noalias(DN_out) = row(DN_DX,i);
                        break;
                    }
                
                double normDn = norm_2(DN_out);
                for(unsigned int i=0; i<TDim+1; i++)
                    if(GetGeometry()[i].Is(BOUNDARY))
                    {
                        rRightHandSideVector[i] += 0.01*source*normDn*Area*(TDim-1); //TODO: check this! it should be TDim*(TDim-1)*N[i] with N[i] on the face and then equal to 1/TDim
                        // using the area as weighting factor is excessive. Reduced it to get a closer to constant gradient between the regions close and far away from the interface
                    }
                
            }   
            
                        
        }
        else //solve an optimization problem with the goal of achievieng a gradient of one for the distance function
        {
            //for debuggin purposes:
            if( dgauss * (this->GetValue(DISTANCE)) < 0.0 )
                std::cout << "Element " << this->Id() << " changed sign while redistancing!!" << std::endl;

            //compute the gradient of the distance
            const array_1d<double,TDim> grad = prod(trans(DN_DX),distances);
            double grad_norm = norm_2(grad);

            //compute RHS ad grad N_i \cdot ( 1/norm_grad * grad - grad) 
            //and multiply everything by grad_norm
            noalias(rRightHandSideVector) = Area*(1.0 - grad_norm)* prod(DN_DX,grad);
            
                        
            //compute the LHS as an approximation of the tangent. 
            //such approximation is taken as a laplacian, which comes from the assumption that the 
            //direction of n does not change when d changes
            //
            //
            //note that the exact tangent could be computed as "P1+P2" with
            //n = grad/grad_norm
            //P1 = (1.0 - 1.0 / (grad_norm + eps) )    * DN_DX * DN_DX.transpose()
            //P2 = 1.0/(grad_norm + eps) * dot(DN_DX * outer(n,n) * DN_DX.transpose() )
            //unfortunately the numerical experiments tell that this in too unstable to be used unless a very 
            //good initial approximation is used
//            noalias(rLeftHandSideMatrix) = (Area*(grad_norm - 1.0))*rod(DN_DX,trans(DN_DX) ); //RISKY!!
            noalias(rLeftHandSideMatrix) = Area*std::max(grad_norm,0.1)*prod( DN_DX,trans(DN_DX) );
        }
    }


    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo) override
    {

        unsigned int number_of_nodes = TDim+1;
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            rResult[i] = GetGeometry()[i].GetDof(DISTANCE).EquationId();
    }

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList,
                            ProcessInfo& rCurrentProcessInfo) override
    {
        unsigned int number_of_nodes = TDim+1;

        if (rElementalDofList.size() != number_of_nodes)
            rElementalDofList.resize(number_of_nodes);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            rElementalDofList[i] = GetGeometry()[i].pGetDof(DISTANCE);

    }


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
    int Check(const ProcessInfo& rCurrentProcessInfo) override
{
        KRATOS_TRY

        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if(ErrorCode != 0) return ErrorCode;
        
        if(this->GetGeometry().size() != TDim+1)
            KRATOS_THROW_ERROR(std::invalid_argument,"wrong number of nodes for element",this->Id());

        // Check that all required variables have been registered
        if(DISTANCE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISTANCE Key is 0. Check if the application was correctly registered.","");

        // Checks on nodes

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(DISTANCE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing DISTANCE variable on solution step data for node ",this->GetGeometry()[i].Id());
        }
 


        return 0;

        KRATOS_CATCH("");
    }

    void Initialize() override
    {
    }
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "should not enter here" << std::endl;
        if (rRightHandSideVector.size() != 0)
	  rRightHandSideVector.resize(0, false);
    }
    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DistanceCalculationElementSimplex #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DistanceCalculationElementSimplex" << TDim << "D";
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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
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
    DistanceCalculationElementSimplex & operator=(DistanceCalculationElementSimplex const& rOther);

    /// Copy constructor.
    DistanceCalculationElementSimplex(DistanceCalculationElementSimplex const& rOther);

    ///@}

}; // Class DistanceCalculationElementSimplex

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim >
inline std::istream& operator >>(std::istream& rIStream,
                                 DistanceCalculationElementSimplex<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const DistanceCalculationElementSimplex<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_DISTANCE_CALCULATION_ELEMENT_H_INCLUDED  defined
