/*
Kratos Multi-Physics

Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
		in the documentation and/or other materials provided with the distribution.
	-	All advertising materials mentioning features or use of this software must display the following acknowledgement: 
			This product includes Kratos Multi-Physics technology.
	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
	
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY 
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2015-09-12
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
    virtual ~DistanceCalculationElementSimplex()
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
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new DistanceCalculationElementSimplex(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Calculate the element's local contribution to the system for the current step.
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int number_of_points = TDim+1;

        boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim > DN_DX;
        array_1d<double, TDim+1 > N;

        if (rLeftHandSideMatrix.size1() != number_of_points)
            rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);

        if (rRightHandSideVector.size() != number_of_points)
            rRightHandSideVector.resize(number_of_points, false);

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
        
        double h = 0.0;
        if(TDim == 2)
            h = sqrt(2.0*Area);
        else
            h = pow(6.0*Area, 0.333333333);
        
        //get distances at the nodes
        array_1d<double, TDim+1 > distances;
        for(unsigned int i=0; i<number_of_points; i++)
        {
            distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }

        const unsigned int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    
        if(step == 1) //solve a poisson problem with a positive/negative heat source depending on the sign of the existing distance function
        {
            //compute distance on gauss point
            const double dgauss = inner_prod(N,distances);
            
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
                        rRightHandSideVector[i] += source*normDn*Area*(TDim-1); //TODO: check this! it should be TDim*(TDim-1)*N[i] with N[i] on the face and then equal to 1/TDim
                        
                    }
                
            }   
            
                        
        }
        else //solve an optimization problem with the goal of achievieng a gradient of one for the distance function
        {
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
            noalias(rLeftHandSideMatrix) = Area*std::max(grad_norm,1e-2/h)*prod( DN_DX,trans(DN_DX) );
        }
        
        
        //apply by penalizaton a constrain to force mantaining the distance to zero over the cut DistanceCalculationElementSimplex
        unsigned int positives = 0, negatives=0;
        for(unsigned int i=0; i<number_of_points; i++)
        {
            if(distances[i] >= 0) positives++;
            else negatives++;
        }
        
        if(positives> 0  && negatives>0) //the element is cut by the interface
        {   
            array_1d<double,3> x0; //point on the cut
            
            //****************************************************************
            //here impose the constraint that the zero is mantained
            //****************************************************************
            //compute a penalty factor by inspecting the diagonal of the LHS
            double penalty = 0.0;
            for(unsigned int i=0; i<TDim+1; i++)
            {
                penalty = std::max( penalty, fabs( rLeftHandSideMatrix(i,i) ) );
            }
            penalty *= 1e2; //1e6;
            
            //now loop over all the edges and 
            for(unsigned int i=0; i<TDim; i++)
            {
                for(unsigned int j=i+1; j<TDim+1; j++)
                {
                    if(distances[i]*distances[j] < 0) //edge is divided
                    {
                        const double Ni = fabs(distances[j])/(fabs(distances[i])+fabs(distances[j]) + 1e-30);
                        const double Nj = 1.0-Ni;
                        
                        rLeftHandSideMatrix(i,i) += penalty*Ni*Ni;
                        rLeftHandSideMatrix(i,j) += penalty*Ni*Nj;
                        rLeftHandSideMatrix(j,i) += penalty*Nj*Ni;
                        rLeftHandSideMatrix(j,j) += penalty*Nj*Nj;
                        
                        noalias(x0) = Ni*GetGeometry()[i].Coordinates() + Nj*GetGeometry()[j].Coordinates();
                    }
                }
            }
            /*
            if(step == 2)
            {*/
                //****************************************************************
                //here impose (more weakly) that exact distance are mantained
                //****************************************************************
                array_1d<double,TDim+1> dexact;
                
                //compute normal
                array_1d<double,TDim> n = prod(trans(DN_DX),distances);
                n /= (norm_2(n) + 1e-30);
                
                //find exact distances
                for(unsigned int i=0; i<TDim; i++)
                {
                    array_1d<double,3> dx = GetGeometry()[i].Coordinates() - x0;
                    dexact[i] = inner_prod(n, dx);
                }
                
                //impose constraint
                penalty *= 0.1;
                for(unsigned int i=0; i<TDim; i++)
                {
                    rLeftHandSideMatrix(i,i) += penalty;
                    rRightHandSideVector(i) += penalty*(dexact[i] - distances[i]);
                }
//             }
            
             ImposeBCs(rLeftHandSideMatrix, rRightHandSideVector, distances);
        }
        
        
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

        unsigned int number_of_nodes = TDim+1;
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            rResult[i] = GetGeometry()[i].GetDof(DISTANCE).EquationId();
        
        
//         for (unsigned int i = 0; i < number_of_nodes; i++)
//             std::cout << rResult[i] << " " ;
//         std::cout << std::endl;
    }

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void GetDofList(DofsVectorType& rElementalDofList,
                            ProcessInfo& rCurrentProcessInfo)
    {
        unsigned int number_of_nodes = TDim+1;

        if (rElementalDofList.size() != number_of_nodes)
            rElementalDofList.resize(number_of_nodes);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            rElementalDofList[i] = GetGeometry()[i].pGetDof(DISTANCE);

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

    void Initialize()
    {
    }
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo)
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DistanceCalculationElementSimplex #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
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
    void ImposeBCs(Matrix& lhs,
                   Vector& rhs,
                   const Vector& distances
         )
    {
        Matrix positive_side_lhs = ZeroMatrix(lhs.size1(), lhs.size2());
        Matrix negative_side_lhs = ZeroMatrix(lhs.size1(), lhs.size2());
        Vector positive_side_rhs = ZeroVector(rhs.size());
        Vector negative_side_rhs = ZeroVector(rhs.size());

        CondenseLocally(lhs,rhs, distances,  1.0, positive_side_lhs, positive_side_rhs);
        CondenseLocally(lhs,rhs, distances,  -1.0, negative_side_lhs, negative_side_rhs);
        noalias(lhs) = positive_side_lhs + negative_side_lhs;
        noalias(rhs) = positive_side_rhs + negative_side_rhs;
    }
    
    // | A11 A12 | x1  =   b1
    // | A21 A22 | x2  =   b2
    //
    //solving the second equation
    //x2 = A22inv * (b2 - A21*x1)
    //
    //substitutiting into the first
    //(A11 - A12*A22inv*A21)*x1 = b1 - A12*A22inv*b2
    //
    void CondenseLocally( const Matrix& A, 
                          const Vector& b,
                          const Vector& distances, 
                          const double sign,
                          Matrix& output,
                          Vector& output_rhs
        
    )
    {
        //count outer nodes to be condensed
        std::vector<unsigned int> inner_indices; 
        inner_indices.reserve(TDim+1);
        std::vector<unsigned int> outer_indices;
        outer_indices.reserve(TDim+1);
        for(unsigned int i=0; i<distances.size(); i++)
        {
            if(distances[i]*sign < 0.0)
                outer_indices.push_back(i);
            else
                inner_indices.push_back(i);
        }
        
        Matrix A11(inner_indices.size(),inner_indices.size(),0.0);
        Matrix A12(inner_indices.size(),outer_indices.size(),0.0);
        Matrix A21(outer_indices.size(),inner_indices.size(),0.0);
        Matrix A22(outer_indices.size(),outer_indices.size(),0.0);
        Vector b1(inner_indices.size());
        Vector b2(outer_indices.size());
        for(unsigned int i=0; i<inner_indices.size(); i++)
        {
            b1[i] = b[inner_indices[i]];
            for(unsigned int j=0; j<inner_indices.size(); j++)
            {
                A11(i,j) = A(inner_indices[i], inner_indices[j]);
                
            }
            for(unsigned int j=0; j<outer_indices.size(); j++)
            {
                A12(i,j) = A(inner_indices[i], outer_indices[j]);
            }
        }

        for(unsigned int i=0; i<outer_indices.size(); i++)
        {
            b2[i] = b[outer_indices[i]];
            for(unsigned int j=0; j<inner_indices.size(); j++)
            {
                A21(i,j) = A(outer_indices[i], inner_indices[j]);
            }
            for(unsigned int j=0; j<outer_indices.size(); j++)
            {
                A22(i,j) = A(outer_indices[i], outer_indices[j]);
            }
        }
        
        double detA22;
        Matrix A22inv(outer_indices.size(), outer_indices.size());
        MathUtils<double>::InvertMatrix(A22, A22inv, detA22); //here i should invert the matrix
        if(detA22 < 1e-30)
        {
            KRATOS_ERROR << "impossible to inverte the local matrix for element " << this->Id() << " matrix is " << A22 << " distances are: " << distances << std::endl;
        }
        
        //A11 -= A12*A22inv*A21;
        Matrix tmp = prod(A22inv,A21);
        noalias(A11) -= prod(A12,tmp);
        
        Vector btmp = prod(A22inv, b2);
        noalias(b1) -= prod( A12,btmp);
                
        //assemble to output matrix
        for(unsigned int i=0; i<inner_indices.size(); i++)
        {
            output_rhs(inner_indices[i]) = b1[i];
            for(unsigned int j=0; j<inner_indices.size(); j++)
            {
                output(inner_indices[i], inner_indices[j]) = A11(i,j);
            }
        }
        
    }

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
