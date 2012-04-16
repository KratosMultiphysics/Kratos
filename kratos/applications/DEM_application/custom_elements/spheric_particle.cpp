//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelspon $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "spheric_particle.h"


namespace Kratos
{

      /// Default constructor.
      SphericParticle::SphericParticle(){}
      
      SphericParticle::SphericParticle(IndexType NewId, const Node<3>::Pointer rNode, PropertiesType::Pointer pProperties){}
      
      SphericParticle::SphericParticle(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties){}
           
      Element::Pointer SphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const 
      {
         return Element::Pointer(new SphericParticle(NewId, ThisNodes, pProperties) );
      } 
      

      /// Destructor.
      SphericParticle::~SphericParticle(){}
      

      void SphericParticle::Initialize(){}
      void SphericParticle::CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::SphericParticle::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo){}
      void SphericParticle::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo){}
      void SphericParticle::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo){}
      void SphericParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo){}
      
      void SphericParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
  
  
}  // namespace Kratos.



