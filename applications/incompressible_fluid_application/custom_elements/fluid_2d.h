/*
==============================================================================
KratosIncompressibleFluidApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-06-28 09:32:36 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_FLUID_2D_H_INCLUDED )
#define  KRATOS_FLUID_2D_H_INCLUDED



// System includes  


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
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

    /// This element implements a Multi-stage element (2D case) to be used in conjuntion with @see FractionalStepStrategy

    /** The element is designed for the solution of the Navier-Stokes equations. Velocity components are considered to be uncoupled, and
     * laplacian formulation is used for the viscous term.
     * OSS (Orthogonal Sub-grid Scale) stabilization is used for both the incompressibility constraint and for the convective term.
     * smagorinsky turbulence model is optionally implemented and controlled by the value of the C_SMAGORINSKY constant, which is passed thorugh the
     * Processinfo.
     * The computation of the "tau" used in the stabilization allows the user to take in account a term depending on 1/Dt
     * this option is controlled by the variable ProcessInfo[DYNAMIC_TAU]. Setting it to 0.0 implies NOT considering a dependence
     * of tau on Dt.
     * The class is organized mainly in 3 stages
     * Stage1 - computes the velocity (designed for non-linear iteration)
     * Stage2 - computes the pressure
     * Stage3 - corrects the velocity taking in account the pressure variation computed in the second step
     */
    class Fluid2D
    : public Element
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer of Fluid2D
        KRATOS_CLASS_POINTER_DEFINITION(Fluid2D);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        Fluid2D(IndexType NewId, GeometryType::Pointer pGeometry);
        Fluid2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /// Destructor.
        virtual ~Fluid2D();


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

        void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

        void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

        void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

        void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo);

        void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

        void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);

        int Check(const ProcessInfo& rCurrentProcessInfo);

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

        virtual std::string Info() const
        {
            return "Fluid2D #";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << Info() << Id();
        }

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
        ///@name Private Operators
        ///@{
        void Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, unsigned int ComponentIndex);
        void Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
        inline double CalculateH(boost::numeric::ublas::bounded_matrix<double, 3,2 > & DN_DX, double Volume);
        inline double CalculateTau(boost::numeric::ublas::bounded_matrix<double, 3,2 > & DN_DX, array_1d<double, 2 > & vel_gauss, const double h, const double nu, const double norm_u, const ProcessInfo& CurrentProcessInfo);
        double ComputeSmagorinskyViscosity(const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX,
                const double& h,
                const double& C,
                const double nu
                );
				
		void AssignConsistentMassMatrixCoefficients(boost::numeric::ublas::bounded_matrix<double, 3,3 >& m);

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
        ///@name Serialization
        ///@{

        friend class Serializer;

        // A private default constructor necessary for serialization

        Fluid2D() : Element()
        {
        }

        virtual void save(Serializer& rSerializer) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        }

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        }

        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        //Fluid2D& operator=(const Fluid2D& rOther);

        /// Copy constructor.
        //Fluid2D(const Fluid2D& rOther);


        ///@}

    }; // Class Fluid2D 

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    /*  inline std::istream& operator >> (std::istream& rIStream,
                                        Fluid2D& rThis);
     */
    /// output stream function
    /*  inline std::ostream& operator << (std::ostream& rOStream,
                                        const Fluid2D& rThis)
        {
          rThis.PrintInfo(rOStream);
          rOStream << std::endl;
          rThis.PrintData(rOStream);

          return rOStream;
        }*/
    ///@}

} // namespace Kratos.

#endif // KRATOS_FLUID_2D_H_INCLUDED  defined 


