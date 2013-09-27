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


#if !defined(KRATOS_FRACTIONAL_STEP_BINGHAM_H_INCLUDED )
#define  KRATOS_FRACTIONAL_STEP_BINGHAM_H_INCLUDED

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
#include "fractional_step.h"
#include "vms.h"
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
    class FractionalStepBingham : public FractionalStep<TDim>
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of FractionalStep
        KRATOS_CLASS_POINTER_DEFINITION(FractionalStepBingham);

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
        
        typedef FractionalStep<TDim> ElementBaseType;
        
        typedef Properties PropertiesType;

        ///@}
        ///@name Life Cycle
        ///@{

        //Constructors.

        /// Default constuctor.
        /**
         * @param NewId Index number of the new element (optional)
         */
        FractionalStepBingham(IndexType NewId = 0) :
            ElementBaseType(NewId)
        {}

        /// Constructor using an array of nodes.
        /**
         * @param NewId Index of the new element
         * @param ThisNodes An array containing the nodes of the new element
         */
        FractionalStepBingham(IndexType NewId, const NodesArrayType& ThisNodes) :
            ElementBaseType(NewId, ThisNodes)
        {}

        /// Constructor using a geometry object.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         */
        FractionalStepBingham(IndexType NewId, GeometryType::Pointer pGeometry) :
            ElementBaseType(NewId, pGeometry)
        {}

        /// Constuctor using geometry and properties.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         * @param pProperties Pointer to the element's properties
         */
        FractionalStepBingham(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
            ElementBaseType(NewId, pGeometry, pProperties)
        {}

        /// Destructor.
        virtual ~FractionalStepBingham()
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /// Create a new element of this type
        /**
         * Returns a pointer to a new FractionalStep element, created using given input
         * @param NewId: the ID of the new element
         * @param ThisNodes: the nodes of the new element
         * @param pProperties: the properties assigned to the new element
         * @return a Pointer to the new element
         */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const
        {
            return Element::Pointer(new FractionalStepBingham(NewId, (this->GetGeometry()).Create(ThisNodes), pProperties));
        }



        /// Obtain a double elemental variable, evaluated on gauss points.
        /**
         * Implemented for a single gauss point only, gets the value of chosen variable in the elemental database
         * @param rVariable Kratos vector variable to compute
         * @param Output Will be filled with the values of the variable on integrartion points
         * @param rCurrentProcessInfo Process info instance
         */
        virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                std::vector<double>& rValues,
                const ProcessInfo& rCurrentProcessInfo)
        {
            if (rVariable == TAUONE || rVariable == TAUTWO || rVariable == MU || rVariable==TAU || rVariable==EQ_STRAIN_RATE)
            {
//                double TauOne, TauTwo;
                double Area;
                
                array_1d<double, TDim+1> N;
                boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim> DN_DX;
                GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

                array_1d<double, 3 > AdvVel;
                this->EvaluateConvVelocity(AdvVel, N);

                double Density, KinViscosity;
                this->EvaluateInPoint(Density, DENSITY, N);
                this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

                double Viscosity;
                double ElemSize = this->ElementSize();
                Viscosity = this->EffectiveViscosity( N, DN_DX, ElemSize, rCurrentProcessInfo);

//                this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Density, Viscosity, rCurrentProcessInfo);

                rValues.resize(1, false);

                 if (rVariable == MU)
                {
                    rValues[0] = Viscosity;
                }
                else if (rVariable == TAU)
                {
                    double NormS = SymmetricGradientNorm(DN_DX);
                    rValues[0] = Density*Viscosity*NormS;
                }
                else if (rVariable == EQ_STRAIN_RATE)
                {
                    rValues[0] = SymmetricGradientNorm(DN_DX);
                }

            }
            else // Default behaviour (returns elemental data)
            {
                rValues.resize(1, false);
                /*
                 The cast is done to avoid modification of the element's data. Data modification
                 would happen if rVariable is not stored now (would initialize a pointer to &rVariable
                 with associated value of 0.0). This is catastrophic if the variable referenced
                 goes out of scope.
                 */
                const FractionalStepBingham<TDim>* const_this = static_cast<const FractionalStepBingham<TDim>*> (this);
                rValues[0] = const_this->GetValue(rVariable);
            }
        }





        ///@}
        ///@name Access
        ///@{

        ///@}
        ///@name Elemental Data
        ///@{



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
            buffer << "FractionalStep #" << this->Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "FractionalStepBingham" << TDim << "D";
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

 
        double EffectiveViscosity(const ShapeFunctionsType &rN,
                                  const ShapeFunctionDerivativesType &rDN_DX,
                                  double ElemSize,
                                  const ProcessInfo &rCurrentProcessInfo)
        {
        const double yield = rCurrentProcessInfo.GetValue(YIELD_STRESS);
        double TotalViscosity = 0.0;
        this->EvaluateInPoint(TotalViscosity,VISCOSITY,rN);
       // TotalViscosity = MolecularViscosity;
        double Density= 0.0;
        this->EvaluateInPoint(Density,DENSITY,rN);       
        if ( yield != 0)
        {
            const double NormS = SymmetricGradientNorm(rDN_DX);

            //non newtonian case
            const double mcoef = rCurrentProcessInfo.GetValue(M);
            double aux_1 = 1.0 - exp(-(mcoef * NormS));
            if(NormS < 1.0/(1000.0*mcoef))
                TotalViscosity += yield*mcoef;
            else
                TotalViscosity += (yield / NormS) * aux_1;
            TotalViscosity/=Density;
        }            
        return TotalViscosity;
            
        }
        
     double SymmetricGradientNorm(const boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim >& rShapeDeriv)
    {
        const unsigned int GradientSize = (TDim*(TDim+1))/2; // Number of different terms in the symmetric gradient matrix
        array_1d<double,GradientSize> GradientVector( GradientSize, 0.0 );
        unsigned int Index,TNumNodes;
        TNumNodes = TDim+1;     
        // Compute Symmetric Grad(u). Note that only the lower half of the matrix is calculated
        for (unsigned int k = 0; k < TNumNodes; ++k)
        {
            const array_1d< double, 3 > & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(VELOCITY);
            Index = 0;
            for (unsigned int i = 0; i < TDim; ++i)
            {
                for (unsigned int j = 0; j < i; ++j) // Off-diagonal
                    GradientVector[Index++] += 0.5 * (rShapeDeriv(k, j) * rNodeVel[i] + rShapeDeriv(k, i) * rNodeVel[j]);
                GradientVector[Index++] += rShapeDeriv(k, i) * rNodeVel[i]; // Diagonal
            }
        }

        // Norm[ Symmetric Grad(u) ] = ( 2 * Sij * Sij )^(1/2)
        Index = 0;
        double NormS(0.0);
        for (unsigned int i = 0; i < TDim; ++i)
        {
            for (unsigned int j = 0; j < i; ++j)
            {
                NormS += 2.0 * GradientVector[Index] * GradientVector[Index]; // Using symmetry, lower half terms of the matrix are added twice
                ++Index;
            }
            NormS += GradientVector[Index] * GradientVector[Index]; // Diagonal terms
            ++Index; // Diagonal terms
        }

        NormS = sqrt( 2.0 * NormS );
        return NormS;
    }       


        /// Helper function to print results on gauss points
        /** Reads a variable from the element's database and returns it in a format
          * that can be used by GetValueOnIntegrationPoints functions.
          * @see GetValueOnIntegrationPoints
          */

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
        FractionalStepBingham & operator=(FractionalStepBingham const& rOther);

        /// Copy constructor.
        FractionalStepBingham(FractionalStepBingham const& rOther);

        ///@}

    }; // Class FractionalStep

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
//    template< unsigned int TDim >
//    inline std::istream& operator >>(std::istream& rIStream,
//                                     FractionalStepBingham<TDim>& rThis)
//    {
//        return rIStream;
//    }
//
//    /// output stream function
//    template< unsigned int TDim >
//    inline std::ostream& operator <<(std::ostream& rOStream,
//                                     const FractionalStepBingham<TDim>& rThis)
//    {
//        rThis.PrintInfo(rOStream);
//        rOStream << std::endl;
//        rThis.PrintData(rOStream);
//
//        return rOStream;
//    }
    ///@}

    ///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_FRACTIONAL_STEP_H_INCLUDED  defined
