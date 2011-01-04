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
//   Date:                $Date: 2010-12-10 11:43:00 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/openmp_utils.h"
#include "utilities/geometry_utilities.h"

#include "fluid_dynamics_application_variables.h"

#ifndef KRATOS_DYNAMIC_SMAGORINSKY_UTILITIES_H_INCLUDED
#define	KRATOS_DYNAMIC_SMAGORINSKY_UTILITIES_H_INCLUDED

namespace Kratos
{
    class DynamicSmagorinskyUtils
    {
    public:

        /// Default Constructor
        DynamicSmagorinskyUtils(ModelPart& rModelPart, ModelPart::ElementsContainerType& rCoarseMesh):
        mrModelPart(rModelPart),
        mrCoarseMesh(rCoarseMesh)
        {}

        /// Destructor
        ~DynamicSmagorinskyUtils()
        {}

        /// Operations

        /// Provide a value for the Smagorinsky coefficient using the Variational Germano Identity
        void CalculateC()
        {

            mrModelPart.GetProcessInfo()[C_SMAGORINSKY] = 0.0; // Disable Smagorinsky to calculate the residual terms

            // Loop over coarse mesh to evaluate all terms that do not involve the fine mesh
            const int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector ElementPartition;
            OpenMPUtils::DivideInPartitions(mrCoarseMesh.size(),NumThreads,ElementPartition);

            double ResC(0.0),ResF(0.0),SmaC(0.0),SmaF(0.0);

            #pragma omp parallel reduction(+:ResC,SmaC)
            {
                int k = OpenMPUtils::ThisThread();
                ModelPart::ElementsContainerType::iterator ElemBegin = mrCoarseMesh.begin() + ElementPartition[k];
                ModelPart::ElementsContainerType::iterator ElemEnd = mrCoarseMesh.begin() + ElementPartition[k+1];

                Vector LocalValues, LocalCoarseVel;
                Matrix LocalMassMatrix;
                ProcessInfo& rProcessInfo = mrModelPart.GetProcessInfo();

                double Residual,Model;

                const unsigned int Dim = ElemBegin->GetGeometry().WorkingSpaceDimension();

                if (Dim == 2)
                {
                    LocalValues.resize(9);
                    LocalCoarseVel.resize(9);
                    LocalMassMatrix.resize(9,9,false);
                    array_1d<double,3> N;
                    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
                    boost::numeric::ublas::bounded_matrix<double,2,2> dv_dx;

                    for( ModelPart::ElementsContainerType::iterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
                    {
                        this->GermanoTerms2D(*itElem,N,DN_DX,dv_dx,LocalValues,LocalCoarseVel,LocalMassMatrix,rProcessInfo,Residual,Model);

                        ResC += Residual;
                        SmaC += Model;
                    }
                }
                else // Dim == 3
                {
                    LocalValues.resize(16);
                    LocalCoarseVel.resize(16);
                    LocalMassMatrix.resize(16,16,false);
                    array_1d<double,4> N;
                    boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
                    boost::numeric::ublas::bounded_matrix<double,3,3> dv_dx;

                    for( ModelPart::ElementsContainerType::iterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
                    {
                        this->GermanoTerms3D(*itElem,N,DN_DX,dv_dx,LocalValues,LocalCoarseVel,LocalMassMatrix,rProcessInfo,Residual,Model);

                        ResC += Residual;
                        SmaC += Model;
                    }
                }
            }

            // Loop over fine mesh to evaluate remaining terms
            OpenMPUtils::DivideInPartitions(mrModelPart.Elements().size(),NumThreads,ElementPartition);

            #pragma omp parallel reduction(+:ResF,SmaF)
            {
                int k = OpenMPUtils::ThisThread();
                ModelPart::ElementsContainerType::iterator ElemBegin = mrModelPart.ElementsBegin() + ElementPartition[k];
                ModelPart::ElementsContainerType::iterator ElemEnd = mrModelPart.ElementsBegin() + ElementPartition[k+1];

                Vector LocalValues, LocalCoarseVel;
                Matrix LocalMassMatrix;
                ProcessInfo& rProcessInfo = mrModelPart.GetProcessInfo();

                double Residual,Model;

                const unsigned int Dim = ElemBegin->GetGeometry().WorkingSpaceDimension();

                if (Dim == 2)
                {
                    LocalValues.resize(9);
                    LocalCoarseVel.resize(9);
                    LocalMassMatrix.resize(9,9,false);
                    array_1d<double,3> N;
                    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
                    boost::numeric::ublas::bounded_matrix<double,2,2> dv_dx;

                    for( ModelPart::ElementsContainerType::iterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
                    {
                        this->GermanoTerms2D(*itElem,N,DN_DX,dv_dx,LocalValues,LocalCoarseVel,LocalMassMatrix,rProcessInfo,Residual,Model);

                        ResC += Residual;
                        SmaC += Model;
                    }
                }
                else // Dim == 3
                {
                    LocalValues.resize(16);
                    LocalCoarseVel.resize(16);
                    LocalMassMatrix.resize(16,16,false);
                    array_1d<double,4> N;
                    boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
                    boost::numeric::ublas::bounded_matrix<double,3,3> dv_dx;

                    for( ModelPart::ElementsContainerType::iterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
                    {
                        this->GermanoTerms3D(*itElem,N,DN_DX,dv_dx,LocalValues,LocalCoarseVel,LocalMassMatrix,rProcessInfo,Residual,Model);

                        ResC += Residual;
                        SmaC += Model;
                    }
                }
            }

            double Num = ResC - ResF;
            double Den = SmaC - SmaF;

            // Once we have all terms, compute the new Smagorinsky coefficient
            if (fabs(Den) < 1.0e-12) // If the result is infinity, disable Smagorinsky (set to zero)
                mrModelPart.GetProcessInfo()[C_SMAGORINSKY] = 0.0;
            else
                mrModelPart.GetProcessInfo()[C_SMAGORINSKY] = sqrt(0.5 * fabs( Num / Den ) );
            KRATOS_WATCH(mrModelPart.GetProcessInfo()[C_SMAGORINSKY])
        }

        /// Calculate the "Coarse Mesh" velocity
        /**
         The operations on the coarse mesh are evaluated on the fine mesh, but using an averaged velocity on the nodes that only exist on the fine mesh.
         Velocity gradients calculated on the fine mesh using this average velocity will be equal to those that would be obtained using the coarse mesh.
         This function assigns the "coarse" velocity value to all nodes
         */
        void SetCoarseVel()
        {
            /*
             Note: This loop can't be parallelized, as we are relying on the fact that refined nodes are at the end of the list and their parents will be updated before the refined nodes
             are reached. There is an alternative solution (always calculate the coarse mesh velocity from the historic database) which can be parallelized but won't work for multiple
             levels of refinement
             */
            for( ModelPart::NodeIterator itNode = mrModelPart.NodesBegin(); itNode != mrModelPart.NodesEnd(); ++itNode)
            {
                if( itNode->GetValue(FATHER_NODES).size() == 2 )
                {
//                    boost::shared_ptr< Node<3> > pParent1 = itNode->GetValue(FATHER_NODES)(0).lock();
                    Node<3>& rParent1 = itNode->GetValue(FATHER_NODES)[0];
//                    boost::shared_ptr< Node<3> > pParent2 = itNode->GetValue(FATHER_NODES)(1).lock();
                    Node<3>& rParent2 = itNode->GetValue(FATHER_NODES)[1];

//                   itNode->GetValue(COARSE_VELOCITY) = 0.5 * ( pParent1->FastGetSolutionStepValue(VELOCITY) + pParent2->FastGetSolutionStepValue(VELOCITY) );
                    itNode->GetValue(COARSE_VELOCITY) = 0.5 * ( rParent1.FastGetSolutionStepValue(VELOCITY) + rParent2.FastGetSolutionStepValue(VELOCITY) );
                }
                else
                {
                    itNode->GetValue(COARSE_VELOCITY) = itNode->FastGetSolutionStepValue(VELOCITY);
                }
            }
        }

        /// Operators

    private:

        void GermanoTerms2D(Element& rElem,
                            array_1d<double,3>& rShapeFunc,
                            boost::numeric::ublas::bounded_matrix<double,3,2>& rShapeDeriv,
                            boost::numeric::ublas::bounded_matrix<double,2,2>& rGradient,
                            Vector& rNodalResidualContainer,
                            Vector& rNodalVelocityContainer,
                            Matrix& rMassMatrix,
                            ProcessInfo& rProcessInfo,
                            double& rResidual,
                            double& rModel)
        {
            const double Dim = 2;
            const double NumNodes = 3;

            // Initialize
            double Area;
            double Density = 0.0;
            rGradient = ZeroMatrix(Dim,Dim);

            rResidual = 0.0;
            rModel = 0.0;

            double Gradij;

            // Calculate the residual
            this->CalculateResidual(rElem,rMassMatrix,rNodalVelocityContainer,rNodalResidualContainer,rProcessInfo); // We use rNodalVelocityContainer as an auxiliaty variable
            this->GetCoarseVelocity2D(rElem,rNodalVelocityContainer);

            for( Vector::iterator itRHS = rNodalResidualContainer.begin(), itVel = rNodalVelocityContainer.begin(); itRHS != rNodalResidualContainer.end(); ++itRHS, ++itVel)
                rResidual += (*itVel) * (*itRHS);

            // Calculate the model term
            GeometryUtils::CalculateGeometryData( rElem.GetGeometry(), rShapeDeriv, rShapeFunc, Area);

            // Compute Grad(u), Density and < Grad(w), Grad(u) >
            for (unsigned int j = 0; j < NumNodes; ++j) // Columns of <Grad(Ni),Grad(Nj)>
            {
                Density += rShapeFunc[j] * rElem.GetGeometry()[j].FastGetSolutionStepValue(DENSITY);
                const array_1d< double,3 >& rNodeVel = rElem.GetGeometry()[j].FastGetSolutionStepValue(VELOCITY); // Nodal velocity

                for (unsigned int i = 0; i < NumNodes; ++i) // Rows of <Grad(Ni),Grad(Nj)>
                {
                    const array_1d< double,3 >& rNodeTest = rElem.GetGeometry()[i].GetValue(COARSE_VELOCITY); // Test function (particularized to coarse velocity)
                    Gradij = 0.0;

                    for (unsigned int k = 0; k < Dim; ++k) // Space Dimensions
                    {
                        Gradij += rShapeDeriv(i,k) * rShapeDeriv(j,k);
                    }
                    rModel += rNodeTest[i] * Gradij * rNodeVel[j];
                }

                for (unsigned int m = 0; m < Dim; ++m)
                    for (unsigned int n = 0; n < Dim; ++n)
                        rGradient(m,n) += rShapeDeriv(j,n) * rNodeVel[m]; // Grad(u)
            }

            rModel *= Area;

            // Norm[ Grad(u) ]
            double SqNorm = 0.0;
            for (unsigned int i = 0; i < Dim; ++i)
                for (unsigned int j = 0; j < Dim; ++j)
                    SqNorm += rGradient(i,j) * rGradient(i,j);

            rModel *= Density * /*pow(Area, 2.0/3.0)*/Area * sqrt(SqNorm);
        }

        void GermanoTerms3D(Element& rElem,
                            array_1d<double,4>& rShapeFunc,
                            boost::numeric::ublas::bounded_matrix<double,4,3>& rShapeDeriv,
                            boost::numeric::ublas::bounded_matrix<double,3,3>& rGradient,
                            Vector& rNodalResidualContainer,
                            Vector& rNodalVelocityContainer,
                            Matrix& rMassMatrix,
                            ProcessInfo& rProcessInfo,
                            double& rResidual,
                            double& rModel)
        {
            const double Dim = 3;
            const double NumNodes = 4;

            // Initialize
            double Volume;
            double Density = 0.0;
            rGradient = ZeroMatrix(Dim,Dim);

            rResidual = 0.0;
            rModel = 0.0;

            double Gradij;

            // Calculate the residual
            this->CalculateResidual(rElem,rMassMatrix,rNodalVelocityContainer,rNodalResidualContainer,rProcessInfo); // We use rNodalVelocityContainer as an auxiliaty variable
            this->GetCoarseVelocity3D(rElem,rNodalVelocityContainer);

            for( Vector::iterator itRHS = rNodalResidualContainer.begin(), itVel = rNodalVelocityContainer.begin(); itRHS != rNodalResidualContainer.end(); ++itRHS, ++itVel)
                rResidual += (*itVel) * (*itRHS);

            // Calculate the model term
            GeometryUtils::CalculateGeometryData( rElem.GetGeometry(), rShapeDeriv, rShapeFunc, Volume);

            // Compute Grad(u), Density and < Grad(w), Grad(u) >
            for (unsigned int j = 0; j < NumNodes; ++j) // Columns of <Grad(Ni),Grad(Nj)>
            {
                Density += rShapeFunc[j] * rElem.GetGeometry()[j].FastGetSolutionStepValue(DENSITY);
                const array_1d< double,3 >& rNodeVel = rElem.GetGeometry()[j].FastGetSolutionStepValue(VELOCITY); // Nodal velocity

                for (unsigned int i = 0; i < NumNodes; ++i) // Rows of <Grad(Ni),Grad(Nj)>
                {
                    const array_1d< double,3 >& rNodeTest = rElem.GetGeometry()[i].GetValue(COARSE_VELOCITY); // Test function (particularized to coarse velocity)
                    Gradij = 0.0;

                    for (unsigned int k = 0; k < Dim; ++k) // Space Dimensions
                    {
                        Gradij += rShapeDeriv(i,k) * rShapeDeriv(j,k);
                    }
                    rModel += rNodeTest[i] * Gradij * rNodeVel[j];
                }

                for (unsigned int m = 0; m < Dim; ++m)
                    for (unsigned int n = 0; n < Dim; ++n)
                        rGradient(m,n) += rShapeDeriv(j,n) * rNodeVel[m]; // Grad(u)
            }

            rModel *= Volume;

            // Norm[ Grad(u) ]
            double SqNorm = 0.0;
            for (unsigned int i = 0; i < Dim; ++i)
                for (unsigned int j = 0; j < Dim; ++j)
                    SqNorm += rGradient(i,j) * rGradient(i,j);

            rModel *= Density * pow(Volume, 2.0/3.0) * sqrt(SqNorm);
        }

        /// 2D version
        void GetCoarseVelocity2D(Element& rElement,
                                 Vector& rVar)
        {
            unsigned int LocalIndex = 0;
            const Element::GeometryType& rGeom = rElement.GetGeometry();

            for (unsigned int itNode = 0; itNode < 3; ++itNode)
            {
                rVar[LocalIndex++] = rGeom[itNode].GetValue(COARSE_VELOCITY_X);
                rVar[LocalIndex++] = rGeom[itNode].GetValue(COARSE_VELOCITY_Y);
                rVar[LocalIndex++] = 0.0; // Pressure Dof
            }
        }

        /// 3D version
        void GetCoarseVelocity3D(Element& rElement,
                                 Vector& rVar)
        {
            unsigned int LocalIndex = 0;
            const Element::GeometryType& rGeom = rElement.GetGeometry();

            for (unsigned int itNode = 0; itNode < 4; ++itNode)
            {
                rVar[LocalIndex++] = rGeom[itNode].GetValue(COARSE_VELOCITY_X);
                rVar[LocalIndex++] = rGeom[itNode].GetValue(COARSE_VELOCITY_Y);
                rVar[LocalIndex++] = rGeom[itNode].GetValue(COARSE_VELOCITY_Z);
                rVar[LocalIndex++] = 0.0; // Pressure Dof
            }
        }

        void CalculateResidual(Element& rElement,
                               Matrix& rMassMatrix, ///@TODO This matrix and the next vector should be transformed to static members once we find a threadsafe way to do so
                               Vector& rAuxVector,
                               Vector& rResidual,
                               ProcessInfo& rCurrentProcessInfo)
        {
            rElement.InitializeNonLinearIteration(rCurrentProcessInfo);

            // Dynamic stabiilzation terms
            rElement.CalculateRightHandSide(rResidual,rCurrentProcessInfo);

            // Dynamic Terms
            rElement.MassMatrix(rMassMatrix,rCurrentProcessInfo);
            rElement.GetSecondDerivativesVector(rAuxVector,0);

            noalias(rResidual) -= prod(rMassMatrix,rAuxVector);

            // Velocity Terms
            rElement.CalculateLocalVelocityContribution(rMassMatrix,rResidual,rCurrentProcessInfo); // Note that once we are here, we no longer need the mass matrix
        }

        /// Data members
        ModelPart& mrModelPart;
        ModelPart::ElementsContainerType& mrCoarseMesh;

    };

}

#endif	/* KRATOS_DYNAMIC_SMAGORINSKY_UTILITIES_H_INCLUDED */

