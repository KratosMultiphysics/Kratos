/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2010-10-09 10:34:00 $
//   Revision:            $Revision: 0.1 $
//
//


#if !defined(KRATOS_FRACTIONAL_STEP_DISCONTINUOUS_H_INCLUDED )
#define  KRATOS_FRACTIONAL_STEP_DISCONTINUOUS_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "fractional_step.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "utilities/discont_utils.h"

// Application includes
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

    /// A stabilized element for the incompressible Navier-Stokes equations.
    /**
     */
    template< unsigned int TDim >
    class FractionalStepDiscontinuous : public FractionalStep<TDim>
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of FractionalStepDiscontinuous
        KRATOS_CLASS_POINTER_DEFINITION(FractionalStepDiscontinuous);

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
        FractionalStepDiscontinuous(IndexType NewId = 0) :
            FractionalStep<TDim>(NewId)
        {}

        /// Constructor using an array of nodes.
        /**
         * @param NewId Index of the new element
         * @param ThisNodes An array containing the nodes of the new element
         */
        FractionalStepDiscontinuous(IndexType NewId, const NodesArrayType& ThisNodes) :
            FractionalStep<TDim>(NewId, ThisNodes)
        {}

        /// Constructor using a geometry object.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         */
        FractionalStepDiscontinuous(IndexType NewId, GeometryType::Pointer pGeometry) :
            FractionalStep<TDim>(NewId, pGeometry)
        {}

        /// Constuctor using geometry and properties.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         * @param pProperties Pointer to the element's properties
         */
        FractionalStepDiscontinuous(IndexType NewId, GeometryType::Pointer pGeometry, Element::PropertiesType::Pointer pProperties) :
            FractionalStep<TDim>(NewId, pGeometry, pProperties)
        {}

        /// Destructor.
        virtual ~FractionalStepDiscontinuous()
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /// Create a new element of this type
        /**
         * Returns a pointer to a new FractionalStepDiscontinuous element, created using given input
         * @param NewId: the ID of the new element
         * @param ThisNodes: the nodes of the new element
         * @param pProperties: the properties assigned to the new element
         * @return a Pointer to the new element
         */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                Element::PropertiesType::Pointer pProperties) const
        {
            return Element::Pointer(new FractionalStepDiscontinuous<TDim>(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
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
/*        virtual int Check(const ProcessInfo& rCurrentProcessInfo);*/


 /**
         * @param rVariable Use DIVPROJ
         * @param rOutput (unused)
         * @param rCurrentProcessInfo Process info instance (unused)
         */
        virtual void Calculate(const Variable<double>& rVariable,
                               double& rOutput,
                               const ProcessInfo& rCurrentProcessInfo);
	
/**
         * @param rVariable Use ADVPROJ or VELOCITY
         * @param Output (unused)
         * @param rCurrentProcessInfo Process info instance (unused)
         */
        virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                               array_1d<double, 3 > & rOutput,
                               const ProcessInfo& rCurrentProcessInfo);

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
            buffer << "FractionalStepDiscontinuous #" << this->Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "FractionalStepDiscontinuous" << TDim << "D";
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
		array_1d<double,(TDim-1)*3 > medge_areas;


        ///@}
        ///@name Protected Operators
        ///@{
        virtual void CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo);

        virtual void AddMomentumSystemTerms(Matrix& rLHSMatrix,
                                    Vector& rRHSVector,
                                    const double Density,
                                    const Vector& rConvOperator,
                                    const array_1d<double,3>& rBodyForce,
                                    const double OldPressure,
                                    const double TauOne,
                                    const double TauTwo,
                                    const array_1d<double,3>& rMomentumProjection,
                                    const double MassProjection,
                                    const ShapeFunctionsType& rN,
                                    const ShapeFunctionDerivativesType& rDN_DX,
                                    const double Weight);
        ///@}
        ///@name Protected Operations
        ///@{

        /// Determine integration point weights and shape funcition derivatives from the element's geometry.
        virtual void CalculateGeometryData(ShapeFunctionDerivativesArrayType& rDN_DX,
                                           Matrix& rNContainer,
                                           Vector& rGaussWeights);





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
        FractionalStepDiscontinuous & operator=(FractionalStepDiscontinuous const& rOther);

        /// Copy constructor.
        FractionalStepDiscontinuous(FractionalStepDiscontinuous const& rOther);

        ///@}

    }; // Class FractionalStepDiscontinuous

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    template< unsigned int TDim >
    inline std::istream& operator >>(std::istream& rIStream,
                                     FractionalStepDiscontinuous<TDim>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const FractionalStepDiscontinuous<TDim>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_FRACTIONAL_STEP_DISCONTINUOUS_H_INCLUDED  defined
