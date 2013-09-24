#include "wall_condition_werner_wengle.h"

namespace Kratos
{
  
  ///@name Specialized implementation of VMS for functions that depend on TDim
  ///@{
  
  /**
   * @see WallConditionWernerWengle::EquationIdVector
   */
  template <>
  void WallConditionWernerWengle<2,2>::EquationIdVector(EquationIdVectorType& rResult,
					     ProcessInfo& rCurrentProcessInfo)
  {
    if ( rCurrentProcessInfo[FRACTIONAL_STEP] == 1 )
      {
        const unsigned int NumNodes = 2;
        const unsigned int LocalSize = 4;
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
	  rResult.resize(LocalSize, false);
	
        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
	  {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
	  }
      }
    else
      {
        rResult.resize(0,false);
      }
  }

  /**
   * @see WallConditionWernerWengle::EquationIdVector
   */
  template <>
  void WallConditionWernerWengle<3,3>::EquationIdVector(EquationIdVectorType& rResult,
					     ProcessInfo& rCurrentProcessInfo)
  {
    if ( rCurrentProcessInfo[FRACTIONAL_STEP] == 1 )
      {
        const SizeType NumNodes = 3;
        const SizeType LocalSize = 9;
        unsigned int LocalIndex = 0;
	
        if (rResult.size() != LocalSize)
	  rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
	  rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
	  rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
	  rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        }
      }
    else
      {
        rResult.resize(0,false);
      }
  }

  /**
   * @see WallConditionWernerWengle::GetDofList
   */
  template <>
  void WallConditionWernerWengle<2,2>::GetDofList(DofsVectorType& rElementalDofList,
				       ProcessInfo& rCurrentProcessInfo)
  {
    if ( rCurrentProcessInfo[FRACTIONAL_STEP] == 1 )
      {
        const SizeType NumNodes = 2;
        const SizeType LocalSize = 4;
	
        if (rElementalDofList.size() != LocalSize)
	  rElementalDofList.resize(LocalSize);
	
        unsigned int LocalIndex = 0;
	
        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
	  {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
	  }
      }
    else
      {
        rElementalDofList.resize(0);
      }
  }

  /**
   * @see WallConditionWernerWengle::GetDofList
   */
  template <>
  void WallConditionWernerWengle<3,3>::GetDofList(DofsVectorType& rElementalDofList,
				       ProcessInfo& rCurrentProcessInfo)
  {
    if ( rCurrentProcessInfo[FRACTIONAL_STEP] == 1 )
      {
        const SizeType NumNodes = 3;
        const SizeType LocalSize = 9;
	
        if (rElementalDofList.size() != LocalSize)
	  rElementalDofList.resize(LocalSize);
	
        unsigned int LocalIndex = 0;
	
        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
	  {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
	  }
      }
    else
      {
        rElementalDofList.resize(0);
      }
  }
  
} // namespace Kratos
