/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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

/* *********************************************************
*
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2008-01-29 13:07:22 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_FACEFORCE3D_H_INCLUDED )
#define  KRATOS_FACEFORCE3D_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


namespace Kratos
{

    class FaceForce3D
                : public Condition
    {
        public:

            // Counted pointer of FaceForce3D
            KRATOS_CLASS_POINTER_DEFINITION( FaceForce3D );


            // Constructor void
            FaceForce3D();

            // Constructor using an array of nodes
            FaceForce3D( IndexType NewId, GeometryType::Pointer pGeometry );

            // Constructor using an array of nodes with properties
            FaceForce3D( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

            // Destructor
            virtual ~FaceForce3D();


            // Name Operations

            Condition::Pointer Create(
                IndexType NewId,
                NodesArrayType const& ThisNodes,
                PropertiesType::Pointer pProperties ) const;

            void EquationIdVector(
                EquationIdVectorType& rResult,
                ProcessInfo& rCurrentProcessInfo );

            void GetDofList(
                DofsVectorType& ElementalDofList,
                ProcessInfo& rCurrentProcessInfo );

            void CalculateRightHandSide(
                VectorType& rRightHandSideVector,
                ProcessInfo& rCurrentProcessInfo );

            void CalculateLocalSystem(
                MatrixType& rLeftHandSideMatrix,
                VectorType& rRightHandSideVector,
                ProcessInfo& rCurrentProcessInfo );

            void MassMatrix(
                MatrixType& rMassMatrix,
                ProcessInfo& rCurrentProcessInfo );

            void DampMatrix(
                MatrixType& rDampMatrix,
                ProcessInfo& rCurrentProcessInfo );

            void GetValuesVector(
                Vector& values,
                int Step = 0 );

            void GetFirstDerivativesVector(
                Vector& values,
                int Step = 0 );

            void GetSecondDerivativesVector(
                Vector& values,
                int Step = 0 );

            /**
             * This function provides the place to perform checks on the completeness of the input.
             * It is designed to be called only once (or anyway, not often) typically at the beginning
             * of the calculations, so to verify that nothing is missing from the input
             * or that no common error is found.
             * @param rCurrentProcessInfo
             */
            virtual int Check( const ProcessInfo& rCurrentProcessInfo );


        protected:


        private:
            ///@name Static Member Variables

            /// privat variables


            // privat name Operations

            void CalculateAll(
                MatrixType& rLeftHandSideMatrix,
                VectorType& rRightHandSideVector,
                const ProcessInfo& rCurrentProcessInfo,
                bool CalculateStiffnessMatrixFlag,
                bool CalculateResidualVectorFlag );

            void CalculateAndSubKp(
                Matrix& K,
                array_1d<double, 3>& ge,
                array_1d<double, 3>& gn,
                const Matrix& DN_De,
                const Vector& N,
                double pressure,
                double weight );

            void MakeCrossMatrix(
                boost::numeric::ublas::bounded_matrix<double, 3, 3>& M,
                array_1d<double, 3>& U );

            void CrossProduct(
                array_1d<double, 3>& cross,
                array_1d<double, 3>& a,
                array_1d<double, 3>& b );

            void SubtractMatrix(
                MatrixType& Destination,
                boost::numeric::ublas::bounded_matrix<double, 3, 3>& InputMatrix,
                int InitialRow,
                int InitialCol );

            void ExpandReducedMatrix(
                Matrix& Destination,
                Matrix& ReducedMatrix );

            void CalculateAndAdd_PressureForce(
                VectorType& residualvector,
                const Vector& N,
                const array_1d<double, 3>& v3,
                double pressure,
                double weight,
                const ProcessInfo& rCurrentProcessInfo );

            ///@}
            ///@name Serialization
            ///@{

            friend class Serializer;

            // A private default constructor necessary for serialization

            virtual void save( Serializer& rSerializer ) const
            {
                KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
            }

            virtual void load( Serializer& rSerializer )
            {
                KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
            }

    }; // class FaceForce3D.

} // namespace Kratos.

#endif // KRATOS_FACEFORCE3D_H_INCLUDED  defined 
