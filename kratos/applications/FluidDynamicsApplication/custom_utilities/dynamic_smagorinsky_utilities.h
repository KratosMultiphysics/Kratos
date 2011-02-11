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
#include <vector>
#include <map>

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "utilities/openmp_utils.h"
#include "utilities/geometry_utilities.h"

#include "fluid_dynamics_application_variables.h"

#ifndef KRATOS_DYNAMIC_SMAGORINSKY_UTILITIES_H_INCLUDED
#define	KRATOS_DYNAMIC_SMAGORINSKY_UTILITIES_H_INCLUDED

namespace Kratos
{
    ///@addtogroup FluidDynamicsApplication
    ///@{

    ///@name Kratos Classes
    ///@{

    /// Helper class to dynamically determine a value for the Smagorinsly parameter.
    /**
     This class uses the Variational Germano Identity to determine a value for the
     Smagorinsky parameter. This value is stored in the elemental variable C_SMAGORINSKY,
     the element implementation is responsible for using it. The ability to assign
     different values to different patches of elements (identified by the PATCH_INDEX
     variable) is supported, although it tends to produce unreliable results due to
     a 0/0 indetermination in patches with smooth velocity fields.

     This class is based in Oberai, A.A. and Wanderer, J., Variational formulation
     of the Germano identity for the Navier Stokes equations, Journal of Turbulence,
     2005, vol 6. Note that the formulation described there requires a nested mesh.
     It takes the model part containing a coarse mesh as input and assumes that
     all elements will be subdivided before CalculateC() is called.

     Remember to call StoreCoarseMesh before refining the element, otherwise the coarse mesh will be lost.

     @see VMS for an element implementation that uses the Smagorinsky model.
     @see Local_Refine_Triangle_Mesh,Local_Refine_Tetrahedra_Mesh for the element refinement process.
     */
    class DynamicSmagorinskyUtils
    {
    public:

        ///@name Life Cycle
        ///@{

        /// Constructor
        /**
         @param rModelPart Reference to the model part containing the coarse mesh
         @param DomainSize Spatial dimension (2 or 3)
         */
        DynamicSmagorinskyUtils(ModelPart& rModelPart, unsigned int DomainSize):
            mrModelPart(rModelPart),
            mDomainSize(DomainSize),
            mCoarseMesh(),
            mPatchIndices()
        {}

        /// Destructor
        ~DynamicSmagorinskyUtils() {}

        ///@}
        ///@name Operations
        ///@{

        /// Store current mesh as coarse mesh. Call before refining.
        /**
         If you are refining more than once, this only has to be called before last refinement.
         */
        void StoreCoarseMesh()
        {
            // Clear existing mesh (if any)
            mCoarseMesh.clear();

            // Store current mesh
            for( ModelPart::ElementsContainerType::ptr_iterator itpElem = mrModelPart.Elements().ptr_begin();
                 itpElem != mrModelPart.Elements().ptr_end(); ++itpElem)
            {
//                (*itpElem)->GetValue(C_SMAGORINSKY) = 0.0; // Set the Smagorinsky parameter to zero for the coarse mesh (do this once to reset any input values)
                mCoarseMesh.push_back(*itpElem);
            }

            // Count the number of patches in the model (in parallel)
            const int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector ElementPartition;
            OpenMPUtils::DivideInPartitions(mCoarseMesh.size(),NumThreads,ElementPartition);

            std::vector< std::vector<int> > LocalIndices(NumThreads);

            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();
                ModelPart::ElementsContainerType::iterator ElemBegin = mCoarseMesh.begin() + ElementPartition[k];
                ModelPart::ElementsContainerType::iterator ElemEnd = mCoarseMesh.begin() + ElementPartition[k+1];

                for( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
                {
                    this->AddNewIndex(LocalIndices[k],itElem->GetValue(PATCH_INDEX));
                }
            }

            // Combine the partial lists and create a map for PATCH_INDEX -> Vector position
            unsigned int Counter = 0;
            std::pair<int, unsigned int> NewVal;
            std::pair< std::map<int, unsigned int>::iterator, bool > Result;
            for( std::vector< std::vector<int> >::iterator itList = LocalIndices.begin(); itList != LocalIndices.end(); ++itList )
            {
                for( std::vector<int>::iterator itIndex = itList->begin(); itIndex != itList->end(); ++itIndex)
                {
                    // Note that instering in map already sorts and checks for uniqueness
                    NewVal.first = *itIndex;
                    NewVal.second = Counter;
                    Result = mPatchIndices.insert(NewVal);
                    if (Result.second)
                        ++Counter;
                }
            }
        }

        /// Provide a value for the Smagorinsky coefficient using the Variational Germano Identity
        void CalculateC()
        {
            // Update the velocity values for the terms that belong to the coarse mesh
            this->SetCoarseVel();

            // Partitioning
            const int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector CoarseElementPartition,FineElementPartition;
            OpenMPUtils::DivideInPartitions(mCoarseMesh.size(),NumThreads,CoarseElementPartition);
            OpenMPUtils::DivideInPartitions(mrModelPart.Elements().size(),NumThreads,FineElementPartition);

            // Initialize temporary containers
            unsigned int PatchNumber = mPatchIndices.size();

            std::vector< std::vector<double> > GlobalPatchNum(NumThreads); // Numerator on each patch
            std::vector< std::vector<double> > GlobalPatchDen(NumThreads); // Denominator on each patch

            const double EnergyTol = 0.05;
            double TotalDissipation = 0;

            #pragma omp parallel reduction(+:TotalDissipation)
            {
                int k = OpenMPUtils::ThisThread();

                // Initialize the iterator boundaries for this thread
                ModelPart::ElementsContainerType::iterator CoarseElemBegin = mCoarseMesh.begin() + CoarseElementPartition[k];
                ModelPart::ElementsContainerType::iterator CoarseElemEnd = mCoarseMesh.begin() + CoarseElementPartition[k+1];

                ModelPart::ElementsContainerType::iterator FineElemBegin = mrModelPart.ElementsBegin() + FineElementPartition[k];
                ModelPart::ElementsContainerType::iterator FineElemEnd = mrModelPart.ElementsBegin() + FineElementPartition[k+1];

                // Initialize some thread-local variables
                Vector LocalValues, LocalCoarseVel;
                Matrix LocalMassMatrix;
                ProcessInfo& rProcessInfo = mrModelPart.GetProcessInfo();

                double Residual,Model;
                unsigned int PatchPosition;

                // Thread-local containers for the values in each patch
                std::vector<double>& rPatchNum = GlobalPatchNum[k];
                std::vector<double>& rPatchDen = GlobalPatchDen[k];
                rPatchNum.resize(PatchNumber,0.0);// Fill with zeros
                rPatchDen.resize(PatchNumber,0.0);

                if (mDomainSize == 2)
                {
                    LocalValues.resize(9);
                    LocalCoarseVel.resize(9);
                    LocalMassMatrix.resize(9,9,false);
                    array_1d<double,3> N;
                    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
                    boost::numeric::ublas::bounded_matrix<double,2,2> dv_dx;

                    // Evaluate the N-S and model terms in each coarse element
                    for( ModelPart::ElementsContainerType::iterator itElem = CoarseElemBegin; itElem != CoarseElemEnd; ++itElem)
                    {
                        PatchPosition = mPatchIndices[ itElem->GetValue(PATCH_INDEX) ];
                        this->GermanoTerms2D(*itElem,N,DN_DX,dv_dx,LocalValues,LocalCoarseVel,LocalMassMatrix,rProcessInfo,Residual,Model);

                        rPatchNum[PatchPosition] += Residual;
                        rPatchDen[PatchPosition] += Model;
                        TotalDissipation += Residual;
                    }

                    // Now evaluate the corresponding terms in the fine mesh
                    for( ModelPart::ElementsContainerType::iterator itElem = FineElemBegin; itElem != FineElemEnd; ++itElem)
                    {
                        // Deactivate Smagorinsky to compute the residual of galerkin+stabilization terms only
                        itElem->GetValue(C_SMAGORINSKY) = 0.0;

                        PatchPosition = mPatchIndices[ itElem->GetValue(PATCH_INDEX) ];
                        this->GermanoTerms2D(*itElem,N,DN_DX,dv_dx,LocalValues,LocalCoarseVel,LocalMassMatrix,rProcessInfo,Residual,Model);

                        rPatchNum[PatchPosition] -= Residual;
                        rPatchDen[PatchPosition] -= Model;
                    }
                }
                else // mDomainSize == 3
                {
                    LocalValues.resize(16);
                    LocalCoarseVel.resize(16);
                    LocalMassMatrix.resize(16,16,false);
                    array_1d<double,4> N;
                    boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
                    boost::numeric::ublas::bounded_matrix<double,3,3> dv_dx;

                    // Evaluate the N-S and model terms in each coarse element
                    for( ModelPart::ElementsContainerType::iterator itElem = CoarseElemBegin; itElem != CoarseElemEnd; ++itElem)
                    {
                        PatchPosition = mPatchIndices[ itElem->GetValue(PATCH_INDEX) ];
                        this->GermanoTerms3D(*itElem,N,DN_DX,dv_dx,LocalValues,LocalCoarseVel,LocalMassMatrix,rProcessInfo,Residual,Model);

                        rPatchNum[PatchPosition] += Residual;
                        rPatchDen[PatchPosition] += Model;
                        TotalDissipation += Residual;
                    }

                    // Now evaluate the corresponding terms in the fine mesh
                    for( ModelPart::ElementsContainerType::iterator itElem = FineElemBegin; itElem != FineElemEnd; ++itElem)
                    {
                        // Deactivate Smagorinsky to compute the residual of galerkin+stabilization terms only
                        itElem->GetValue(C_SMAGORINSKY) = 0.0;

                        PatchPosition = mPatchIndices[ itElem->GetValue(PATCH_INDEX) ];
                        this->GermanoTerms3D(*itElem,N,DN_DX,dv_dx,LocalValues,LocalCoarseVel,LocalMassMatrix,rProcessInfo,Residual,Model);

                        rPatchNum[PatchPosition] -= Residual;
                        rPatchDen[PatchPosition] -= Model;
                    }
                }
            }

            // Combine the results of each thread in position 0
            for( std::vector< std::vector<double> >::iterator itNum = GlobalPatchNum.begin()+1, itDen = GlobalPatchDen.begin()+1;
                 itNum != GlobalPatchNum.end(); ++itNum, ++itDen)
            {
                for( std::vector<double>::iterator TotalNum = GlobalPatchNum[0].begin(), LocalNum = itNum->begin(),
                        TotalDen = GlobalPatchDen[0].begin(), LocalDen = itDen->begin();
                     TotalNum != GlobalPatchNum[0].end(); ++TotalNum,++LocalNum,++TotalDen,++LocalDen)
                {
                    *TotalNum += *LocalNum;
                    *TotalDen += *LocalDen;
                }
            }

            // Compute the smagorinsky coefficient for each patch by combining the values from each thread
            std::vector<double> PatchC(PatchNumber);
            double NumTol = EnergyTol * fabs(TotalDissipation);
            for( std::vector<double>::iterator itNum = GlobalPatchNum[0].begin(), itDen = GlobalPatchDen[0].begin(), itC = PatchC.begin();
                 itC != PatchC.end(); ++itNum, ++itDen, ++itC)
            {
                // If the dissipation we are "missing" by not considering Smagorinsky is small, do not use Smagorinsky (this avoids a division by ~0, as the denominator should go to zero too)
                if ( (fabs(*itNum) < NumTol) )//|| (fabs(*itDen) < 1.0e-12) )
                    *itC = 0.0;
                else
                    *itC = sqrt( 0.5 * fabs( *itNum / *itDen ) );
            }

            // Finally, assign each element its new smagorinsky value
            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();
                ModelPart::ElementsContainerType::iterator ElemBegin = mrModelPart.ElementsBegin() + FineElementPartition[k];
                ModelPart::ElementsContainerType::iterator ElemEnd = mrModelPart.ElementsBegin() + FineElementPartition[k+1];

                unsigned int PatchPosition;

                for( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
                {
                    PatchPosition = mPatchIndices[ itElem->GetValue(PATCH_INDEX) ];
                    itElem->GetValue(C_SMAGORINSKY) = PatchC[PatchPosition];
                }
            }
        }

        /// For the bridge analysis problem, correct the boundary flag after the refinement.
        /**
         Remember to run this AFTER EACH REFINEMENT STEP
         Possible values for the variable: 1.0 inlet, 2.0 bridge surface, 3.0 outlet, 0.0 otherwise
         @param rThisVariable The Kratos variable used to identify the boundary
         */
        void CorrectFlagValues(Variable<double>& rThisVariable = FLAG_VARIABLE)
        {
            // Loop over coarse mesh to evaluate all terms that do not involve the fine mesh
            const int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector NodePartition;
            OpenMPUtils::DivideInPartitions(mrModelPart.NumberOfNodes(),NumThreads,NodePartition);

            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();
                ModelPart::NodeIterator NodesBegin = mrModelPart.NodesBegin() + NodePartition[k];
                ModelPart::NodeIterator NodesEnd = mrModelPart.NodesBegin() + NodePartition[k+1];

                double Value0, Value1;

                for( ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
                {
                    if( itNode->GetValue(FATHER_NODES).size() == 2 ) // If the node is refined
                    {
                        Value0 = itNode->GetValue(FATHER_NODES)[0].FastGetSolutionStepValue(rThisVariable);
                        Value1 = itNode->GetValue(FATHER_NODES)[1].FastGetSolutionStepValue(rThisVariable);

                        if( Value0 != Value1 ) // If this node is problematic
                        {
                            if ( Value0 == 0.0 || Value1 == 0.0 )
                            {
                                // if either of the parents is not on the boundary, this node is not on the boundary
                                itNode->FastGetSolutionStepValue(rThisVariable) = 0.0;
                            }
                            /* All remaining cases are unlikely in well-posed problems,
                             I'm arbitrarily giving priority to the outlet,
                             so that the node is only inlet or bridge surface if both parents are
                             */
                            else if( Value0 == 3.0 )
                            {
                                itNode->FastGetSolutionStepValue(rThisVariable) = Value0;
                            }
                            else if( Value1 == 3.0 )
                            {
                                // The node is only bridge surface if both parents are
                                itNode->FastGetSolutionStepValue(rThisVariable) = Value1;
                            }
                            else // Default behaviour: Parent 0 takes precedence
                            {
                                itNode->FastGetSolutionStepValue(rThisVariable) = Value0;
                            }
                        }
                    }
                }
            }
        }

        ///@}

    private:

        ///@name Member Variables
        ///@{

        /// ModelPart of the fluid problem
        ModelPart& mrModelPart;
        /// Spatial dimenstion
        unsigned int mDomainSize;
        /// Container for the coarse mesh (the fine mesh is stored by the model part)
        ModelPart::ElementsContainerType mCoarseMesh;
        /// A map relating patch indices to positions in the internal storage arrays
        std::map<int, unsigned int> mPatchIndices;

        ///@name Private Operations
        ///@{

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
                    Node<3>& rParent1 = itNode->GetValue(FATHER_NODES)[0];
                    Node<3>& rParent2 = itNode->GetValue(FATHER_NODES)[1];

                    itNode->GetValue(COARSE_VELOCITY) = 0.5 * ( rParent1.FastGetSolutionStepValue(VELOCITY) + rParent2.FastGetSolutionStepValue(VELOCITY) );
                }
                else
                {
                    itNode->GetValue(COARSE_VELOCITY) = itNode->FastGetSolutionStepValue(VELOCITY);
                }
            }
        }

        /// Return the Galerkin (+stabilization) and Model terms for this element (2D version)
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

                    for (unsigned int k = 0; k < Dim; ++k) // Space Dimensions
                        rModel += rNodeTest[k] * rShapeDeriv(i,k) * rShapeDeriv(j,k) * rNodeVel[k];
                }

                for (unsigned int m = 0; m < Dim; ++m) // Calculate symmetric gradient
                {
                    for (unsigned int n = 0; n < m; ++n) // Off-diagonal
                        rGradient(m,n) += 0.5 * (rShapeDeriv(j,n) * rNodeVel[m] + rShapeDeriv(j,m) * rNodeVel[n]); // Symmetric gradient, only lower half is written
                    rGradient(m,m) += rShapeDeriv(j,m) * rNodeVel[m]; // Diagonal
                }
            }

            rModel *= Area; // To this point, rModel contains the integral over the element of Grad(U_coarse):Grad(U)

            // Norm[ Grad(u) ]
            double SqNorm = 0.0;
            for (unsigned int i = 0; i < Dim; ++i)
            {
                for (unsigned int j = 0; j < i; ++j)
                    SqNorm += 2.0 * rGradient(i,j) * rGradient(i,j); // Adding off-diagonal terms (twice, as matrix is symmetric)
                SqNorm += rGradient(i,i) * rGradient(i,i); // Diagonal terms
            }

            // "Fixed" part of Smagorinsky viscosity: Density * FilterWidth^2 * Norm(SymmetricGrad(U)). 2*C^2 is accounted for in the caller function
            const double sqH = 2*Area;
            rModel *= Density * sqH * sqrt(SqNorm);
        }

        /// Return the Galerkin (+stabilization) and Model terms for this element (3D version)
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

                    for (unsigned int k = 0; k < Dim; ++k) // Space Dimensions
                        rModel += rNodeTest[k] * rShapeDeriv(i,k) * rShapeDeriv(j,k) * rNodeVel[k];
                }

                for (unsigned int m = 0; m < Dim; ++m) // Calculate symmetric gradient
                {
                    for (unsigned int n = 0; n < m; ++n) // Off-diagonal
                        rGradient(m,n) += 0.5 * (rShapeDeriv(j,n) * rNodeVel[m] + rShapeDeriv(j,m) * rNodeVel[n]); // Symmetric gradient, only lower half is written
                    rGradient(m,m) += rShapeDeriv(j,m) * rNodeVel[m]; // Diagonal
                }
            }

            rModel *= Volume; // To this point, rModel contains the integral over the element of Grad(U_coarse):Grad(U)

            // Norm[ Grad(u) ]
            double SqNorm = 0.0;
            for (unsigned int i = 0; i < Dim; ++i)
            {
                for (unsigned int j = 0; j < i; ++j)
                    SqNorm += 2.0 * rGradient(i,j) * rGradient(i,j); // Adding off-diagonal terms (twice, as matrix is symmetric)
                SqNorm += rGradient(i,i) * rGradient(i,i); // Diagonal terms
            }

            const double cubeH = 6*Volume;
            rModel *= Density * pow(cubeH, 2.0/3.0) * sqrt(SqNorm);
        }

        /// Equivalent to VMS2DSmagorinsky::GetFirstDerivativesVector(), using the velocity evaluated on the coarse mesh
        void GetCoarseVelocity2D(Element& rElement,
                                 Vector& rVar)
        {
            unsigned int LocalIndex = 0;
            const Element::GeometryType& rGeom = rElement.GetGeometry();

            for (unsigned int itNode = 0; itNode < 3; ++itNode)
            {
                const array_1d< double,3>& rCoarseVel = rGeom[itNode].GetValue(COARSE_VELOCITY);
                rVar[LocalIndex++] = rCoarseVel[0];
                rVar[LocalIndex++] = rCoarseVel[1];
                rVar[LocalIndex++] = 0.0; // Pressure Dof
            }
        }

        /// Equivalent to VMS3DSmagorinsky::GetFirstDerivativesVector(), using the velocity evaluated on the coarse mesh
        void GetCoarseVelocity3D(Element& rElement,
                                 Vector& rVar)
        {
            unsigned int LocalIndex = 0;
            const Element::GeometryType& rGeom = rElement.GetGeometry();

            for (unsigned int itNode = 0; itNode < 4; ++itNode)
            {
                const array_1d< double,3>& rCoarseVel = rGeom[itNode].GetValue(COARSE_VELOCITY);
                rVar[LocalIndex++] = rCoarseVel[0];
                rVar[LocalIndex++] = rCoarseVel[1];
                rVar[LocalIndex++] = rCoarseVel[2];
                rVar[LocalIndex++] = 0.0; // Pressure Dof
            }
        }

        /// Call the element's member functions to obtain its residual
        void CalculateResidual(Element& rElement,
                               Matrix& rMassMatrix, ///@todo This matrix and the next vector should be transformed to static members once we find a threadsafe way to do so
                               Vector& rAuxVector,
                               Vector& rResidual,
                               ProcessInfo& rCurrentProcessInfo)
        {
            rElement.InitializeNonLinearIteration(rCurrentProcessInfo);

            // Dynamic stabilization terms
            rElement.CalculateRightHandSide(rResidual,rCurrentProcessInfo);

            // Dynamic Terms
            rElement.MassMatrix(rMassMatrix,rCurrentProcessInfo);
            rElement.GetSecondDerivativesVector(rAuxVector,0);

            noalias(rResidual) -= prod(rMassMatrix,rAuxVector);

            // Velocity Terms
            rElement.CalculateLocalVelocityContribution(rMassMatrix,rResidual,rCurrentProcessInfo); // Note that once we are here, we no longer need the mass matrix
        }

        /// Check if a patch index is known
        void AddNewIndex( std::vector<int>& rIndices,
                          int ThisIndex )
        {
            bool IsNew = true;
            for( std::vector<int>::iterator itIndex = rIndices.begin(); itIndex != rIndices.end(); ++itIndex)
            {
                if( ThisIndex == *itIndex)
                {
                    IsNew = false;
                    break;
                }
            }

            if (IsNew)
                rIndices.push_back(ThisIndex);
        }

        ///@} // Private operations

    };

    ///@} Kratos classes

    ///@} Application group

} // namespace Kratos

#endif	/* KRATOS_DYNAMIC_SMAGORINSKY_UTILITIES_H_INCLUDED */
