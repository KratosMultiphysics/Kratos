#include "autoslip_inlet_3d.h"
#include "pfem_2_application_variables.h"


namespace Kratos
{


void MonolithicAutoSlipInlet3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo)
    {
        const SizeType BlockSize = 1;
        unsigned int TNumNodes = GetGeometry().size();
        const SizeType LocalSize = BlockSize * TNumNodes;

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize,LocalSize);

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);
        
        if (this->GetValue(IS_INLET) == 1.0) //for inlet faces:
		{
			 const double factor = 1.0/3.0; //or 2 in 2d (2 nodes in the condition)
			
			 double velocity = GetValue(INLET_VELOCITY);
			 array_1d<double,3> normal =  ZeroVector(3);
			 CalculateNormal(normal);
			 
			 double area = 0.0;
			 for (unsigned int i=0; i!=3 ; i++)  //or 2 in 2d
			 	 area+=pow(normal(i),2);
			 
			 area = sqrt(area);
			 
			 const double added_RHS_boundary_term = - velocity * area; //inlet!
			 	 
			 for (unsigned int i=0; i!=3 ; i++)  //3 nodes of the triangle
			     rRightHandSideVector(i) = factor * added_RHS_boundary_term;
				 
			 
		}
    }

void MonolithicAutoSlipInlet3D::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        const ProcessInfo& rCurrentProcessInfo)
{
        const SizeType BlockSize = 1;
        unsigned int TNumNodes = GetGeometry().size();
        const SizeType LocalSize = BlockSize * TNumNodes;

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize);

        noalias(rRightHandSideVector) = ZeroVector(LocalSize);
}
    
	
void MonolithicAutoSlipInlet3D::EquationIdVector(EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 3;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        //rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        //rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        //rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

void MonolithicAutoSlipInlet3D::GetDofList(DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 3;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        //rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        //rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        //rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}


// protected funcions

void MonolithicAutoSlipInlet3D::CalculateNormal(array_1d<double,3>& An )
{
    Geometry<Node >& pGeometry = this->GetGeometry();
    

    array_1d<double,3> v1,v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;
    
    //now checking orientation using the normal:
    const SizeType NumNodes = 3;
    array_1d<double,3> nodal_normal =  ZeroVector(3);
    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
		nodal_normal += pGeometry[iNode].FastGetSolutionStepValue(NORMAL);
    
    double dot_prod = nodal_normal[0]*An[0] + nodal_normal[1]*An[1] + nodal_normal[2]*An[2];
    if (dot_prod<0.0)
    {
		//std::cout << "inverting the normal of a MonolithicAutoSlipInlet3D!!" << std::endl;
		An *= -1.0; // inverting the direction of the normal!!!
	}
    
}

} // namespace Kratos
