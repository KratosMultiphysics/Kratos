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
//   Date:                $Date: 2010-10-11 10:42:00 $
//   Revision:            $Revision: 0.1 $
//
//


#if !defined(KRATOS_VMS3D_H_INCLUDED )
#define  KRATOS_VMS3D_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "vms_base.h"
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/geometry_utilities.h"
#include "includes/variables.h"
#include "../fluid_dynamics_application_variables.h"

#define DIM 3 // Spatial dimension
#define NUMNODES 4 // NumNodes = Dim + 1 (trianlges)
#define BLOCKSIZE 4 // BlockSize = Dim + 1 (Dim Velocity Dofs + 1 Pressure Dof)
#define LOCALSIZE 16 // LocalSize = NumNodes * BlockSize (Size of elemental contributions)

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

    /// A stabilized element for the incompressible Navier-Stokes equations
    /**
     * This class implements a triangular 3D element with a stabilized formulation
     * in the Variational Multiscale framework. The the subscales can be modeled
     * using either Algebraic Subgird Scales (ASGS) or Orthogonal Subscales (OSS).
     * In the case of OSS, the projection terms are treated explicitly (computed
     * using the results of the previous iteration) and the subscales are not
     * tracked in time. The choice of subscale model is made based on the Process Info
     * variable OSS_SWITCH (OSS if 1.0, ASGS otherwise)
     */
    class VMS3D : public VMSBase<DIM>
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of VMS3D
        KRATOS_CLASS_POINTER_DEFINITION(VMS3D);

        ///base type: an IndexedObject that automatically has a unique number
        typedef IndexedObject BaseType;

        ///definition of node type (default is: Node<3>)
        typedef Node < 3 > NodeType;

        /**
         * Properties are used to store any parameters
         * related to the constitutive law
         */
        typedef Properties PropertiesType;

        ///definition of the geometry type with given NodeType
        typedef Geometry<NodeType> GeometryType;

        ///definition of nodes container type, redefined from GeometryType
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

        typedef Vector VectorType;

        typedef Matrix MatrixType;

        typedef std::size_t IndexType;

        typedef std::size_t SizeType;

        typedef std::vector<std::size_t> EquationIdVectorType;

        typedef std::vector< Dof<double>::Pointer > DofsVectorType;

        typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

        typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

        ///Type definition for integration methods
        typedef GeometryData::IntegrationMethod IntegrationMethod;

        ///@}
        ///@name Life Cycle
        ///@{

        ///Constructors.

        /**
         * Default constuctor.
         * @param NewId Index number of the new element (optional)
         */
        VMS3D(IndexType NewId = 0) :
            VMSBase<DIM>(NewId)
        {}

        /**
         * Constructor using an array of nodes.
         * @param NewId Index of the new element
         * @param ThisNodes An array containing the nodes of the new element
         */
        VMS3D(IndexType NewId, const NodesArrayType& ThisNodes) :
            VMSBase<DIM>(NewId, ThisNodes)
        {}

        /**
         * Constructor using a geometry object.
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         */
        VMS3D(IndexType NewId, GeometryType::Pointer pGeometry) :
            VMSBase<DIM>(NewId, pGeometry)
        {}

        /**
         * Constuctor using geometry and properties.
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         * @param pProperties Pointer to the element's properties
         */
        VMS3D(IndexType NewId, GeometryType::Pointer pGeometry,
                 PropertiesType::Pointer pProperties) :
            VMSBase<DIM>(NewId, pGeometry, pProperties)
        {}

        /// Destructor.
        virtual ~VMS3D()
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /// Create a new element of this type
        /**
         * Returns a pointer to a new VMS3D element, created using given input
         * @param NewId: the ID of the new element
         * @param ThisNodes: the nodes of the new element
         * @param pProperties: the properties assigned to the new element
         * @return a Pointer to the new element
         */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const
        {
            KRATOS_TRY
            return Element::Pointer(new VMS3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
            KRATOS_CATCH("");
        }

        /// Provides the global indices for each one of this element's local rows
        /**
         * this determines the elemental equation ID vector for all elemental
         * DOFs
         * @param rResult: A vector containing the global Id of each row
         * @param rCurrentProcessInfo: the current process info object
         */
        virtual void EquationIdVector(EquationIdVectorType& rResult, // Note: requires different implementations depending on space dimension
                                      ProcessInfo& rCurrentProcessInfo)
        {
            unsigned int LocalIndex = 0;

            if (rResult.size() != LOCALSIZE)
                rResult.resize(LOCALSIZE, false);

            for (unsigned int iNode = 0; iNode < NUMNODES; ++iNode)
            {
                rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
                rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
                rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
                rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
            }
        }

        /// Returns a list of the element's Dofs
        /**
         * @param ElementalDofList: the list of DOFs
         * @param rCurrentProcessInfo: the current process info instance
         */
        virtual void GetDofList(DofsVectorType& rElementalDofList, // Note: requires different implementations depending on space dimension
                                ProcessInfo& rCurrentProcessInfo)
        {
            if (rElementalDofList.size() != LOCALSIZE)
                rElementalDofList.resize(LOCALSIZE);

            unsigned int LocalIndex = 0;

            for (unsigned int iNode = 0; iNode < NUMNODES; ++iNode)
            {
                rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
                rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
                rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
                rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
            }
        }

        /// Returns vx, vy, p for each node
        virtual void GetFirstDerivativesVector(Vector& Values, int Step = 0) // Note: requires different implementations depending on space dimension
        {
            unsigned int LocalIndex = 0;

            if (Values.size() != LOCALSIZE)
                Values.resize(LOCALSIZE, false);

            for (unsigned int iNode = 0; iNode < NUMNODES; ++iNode)
            {
                Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY_X, Step);
                Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY_Y, Step);
                Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY_Z, Step);
                Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
            }
        }

        /// Returns ax, ay, 0 for each node
        virtual void GetSecondDerivativesVector(Vector& Values, int Step = 0) // Note: requires different implementations depending on space dimension
        {
            unsigned int LocalIndex = 0;

            if (Values.size() != LOCALSIZE)
                Values.resize(LOCALSIZE, false);

            for (unsigned int iNode = 0; iNode < NUMNODES; ++iNode)
            {
                Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION_X, Step);
                Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION_Y, Step);
                Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION_Z, Step);
                Values[LocalIndex++] = 0.0;
            }
        }

        /// Implementation of GetValueOnIntegrationPoints to obtain the vorticity
        /**
         * Computes the vorticity (rotational of the velocity) for the current velocity values
         * @param rVariable: Kratos vector variable to compute (only implemented for VORTICITY)
         * @param Output: Values of vorticity on integrartion points
         * @param rCurrentProcessInfo: Process info instance
         */
        virtual void GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable, // Note: requires different implementations depending on space dimension
                                                  std::vector<array_1d<double,3> >& rOutput,
                                                  const ProcessInfo& rCurrentProcessInfo)
        {
            if (rVariable == VORTICITY)
            {
                // Set output vector (for a single integration point)
                rOutput.resize(1);
                array_1d<double, 3 > & rVorticity = rOutput[0];
                rVorticity.resize(3);
                rVorticity[0] = 0.0;
                rVorticity[1] = 0.0;
                rVorticity[2] = 0.0;

                double Area;
                array_1d<double, NUMNODES> N;
                boost::numeric::ublas::bounded_matrix<double, NUMNODES, DIM> DN_DX;
                GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

                for (unsigned int iNode = 0; iNode < NUMNODES; ++iNode)
                {
                    const array_1d<double, 3 > & rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
                    rVorticity[0] += N[iNode] * ( DN_DX(iNode,1)*rVelocity[2] - DN_DX(iNode,2)*rVelocity[1] );
                    rVorticity[1] += N[iNode] * ( DN_DX(iNode,2)*rVelocity[0] - DN_DX(iNode,0)*rVelocity[2] );
                    rVorticity[2] += N[iNode] * ( DN_DX(iNode,0)*rVelocity[1] - DN_DX(iNode,1)*rVelocity[0] );
                }
                rVorticity *= 0.5; // vorticity = 1/2 (nabla x velocity)
            }
        }

        /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
        virtual void GetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                                  std::vector<double>& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo)
        {
            double TauOne,TauTwo;
            double Area;
            array_1d<double, NUMNODES> N;
            boost::numeric::ublas::bounded_matrix<double, NUMNODES, DIM> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            array_1d<double,3> AdvVel;
            GetAdvectiveVel(AdvVel,N);

            double KinViscosity;
            GetPointContribution(KinViscosity,VISCOSITY,N);

            CalculateTau(TauOne,TauTwo,AdvVel,Area,KinViscosity,rCurrentProcessInfo);

            rValues.resize(1);
            if (rVariable == TAUONE)
            {
                rValues[0] = TauOne;
            }
            else if (rVariable == TAUTWO)
            {
                rValues[0] = TauTwo;
            }
        }

        /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
        virtual void GetValueOnIntegrationPoints( const Variable<array_1d<double,6> >& rVariable,
                                                  std::vector<array_1d<double,6> >& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo)
        {}

        /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
        virtual void GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
                                                  std::vector<Vector>& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo)
        {}

        /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
        virtual void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
                                                  std::vector<Matrix>& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo)
        {}

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
            std::stringstream buffer;
            buffer << "VMS3D #" << Id();
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "VMS3D #" << Id();
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
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{

        double ElementSize(const double Volume)
        {
            return 0.60046878 * pow(Volume,0.333333333333333333333); // Diameter of sphere circumscribed to regular tetrahedron of given volume
        }


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
        VMS3D & operator=(VMS3D const& rOther);

        /// Copy constructor.
        VMS3D(VMS3D const& rOther);

        ///@}

    }; // Class VMS3D

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function

    inline std::istream & operator >>(std::istream& rIStream,
                                      VMS3D& rThis)
    {
        return rIStream;
    }

    /// output stream function

    inline std::ostream & operator <<(std::ostream& rOStream,
                                      const VMS3D& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}


} // namespace Kratos.

#undef DIM // Spatial dimension
#undef NUMNODES // NumNodes = Dim + 1 (trianlges)
#undef BLOCKSIZE // BlockSize = Dim + 1 (Dim Velocity Dofs + 1 Pressure Dof)
#undef LOCALSIZE // LocalSize = NumNodes * BlockSize (Size of elemental contributions)

#endif // KRATOS_VMS3D_H_INCLUDED  defined
