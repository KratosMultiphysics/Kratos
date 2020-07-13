//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:               June 2018 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_element.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

template <unsigned int TDim>
Element::Pointer TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{
  KRATOS_TRY;

  TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
  return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement(NewElement));

  KRATOS_CATCH("");
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
  KRATOS_TRY;

  // std::cout<<"InitializeNonLinearIteration "<<std::endl;
  GeometryType &rGeom = this->GetGeometry();
  const unsigned int NumNodes = rGeom.PointsNumber();
  // Shape functions and integration points
  ShapeFunctionDerivativesArrayType DN_DX;
  Matrix NContainer;
  VectorType GaussWeights;
  this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
  // ElementalVariables rElementalVariables;
  // this->InitializeElementalVariables(rElementalVariables);
  // Loop on integration points
  const ShapeFunctionDerivativesType &rDN_DX = DN_DX[0];

  // const unsigned int NumGauss = GaussWeights.size();
  // if(NumGauss==1){
  // 	this->CalcElementalStrains(rElementalVariables,rCurrentProcessInfo,rDN_DX);
  // }else{
  // 	std::cout<<"a different structure is required for more gauss points"<<std::endl;
  // }

  double meanElementEdgesLength = this->ElementSize();
  double elementVolume = 0;
  if (TDim == 3)
  {
    elementVolume = rGeom.Volume() * 0.25;
  }
  else if (TDim == 2)
  {
    elementVolume = rGeom.Area() / 3.0;
  }

  if (this->Is(SOLID))
  {
    for (unsigned int i = 0; i < NumNodes; i++)
    {
      // if(rGeom[i].FastGetSolutionStepValue(INTERFACE_NODE)==false){
      //   VectorType solidNodalSFDneighbours=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
      //   unsigned int solidNodalSFDneighboursSize=solidNodalSFDneighbours.size();
      //   // std::cout<<"SOLID_NODAL_SFD_NEIGHBOURS_ORDER) "<<rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER)<<std::endl;

      //   if(solidNodalSFDneighboursSize>1)
      //   {
      //     double& solidMeshSize = rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
      //     ElementWeakPtrVectorType& neighb_elems = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
      //     double numberOfNeighElems=double(neighb_elems.size());
      //     solidMeshSize+=meanElementEdgesLength/numberOfNeighElems;

      //     double solidNodalVolume=rGeom[i].FastGetSolutionStepValue(NODAL_VOLUME);

      //     for (unsigned int j = 0; j< NumNodes; j++)
      //     {
      //         unsigned int idNodeOfConsideredElement=rGeom[j].Id();
      //         unsigned int solidSFDposition=0;
      //         for (unsigned int k = 0; k< solidNodalSFDneighboursSize; k++)
      //         {
      //           if(idNodeOfConsideredElement==solidNodalSFDneighbours[k])
      //           {
      //             rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[solidSFDposition]   += rDN_DX(j,0)*elementVolume/solidNodalVolume;
      //             rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[solidSFDposition+1] += rDN_DX(j,1)*elementVolume/solidNodalVolume;
      //             if(TDim==3){
      //               rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[solidSFDposition+2] += rDN_DX(j,2)*elementVolume/solidNodalVolume;
      //             }
      //             break;
      //           }
      //           solidSFDposition+=TDim;
      //         }
      //     }
      //   }
      //   else
      //   {
      //       std::cout<<rGeom[i].Id()<<"  this solid node is isolated!!! "<<std::endl;
      //       for (unsigned int k = 0; k< TDim; k++)
      //       {
      //         rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD)(k,k)=1.0;
      //       }
      //   }
      // }
      // else
      // {
      VectorType solidNodalSFDneighbours = rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
      unsigned int solidNodalSFDneighboursSize = solidNodalSFDneighbours.size();
      // std::cout<<"SOLID_NODAL_SFD_NEIGHBOURS_ORDER) "<<rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER)<<std::endl;

      if (solidNodalSFDneighboursSize > 1)
      {
        double &solidMeshSize = rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_MEAN_MESH_SIZE);
        ElementWeakPtrVectorType &neighb_elems = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
        double numberOfNeighElems = double(neighb_elems.size());
        solidMeshSize += meanElementEdgesLength / numberOfNeighElems;

        double solidNodalVolume = rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_VOLUME);

        for (unsigned int j = 0; j < NumNodes; j++)
        {
          unsigned int idNodeOfConsideredElement = rGeom[j].Id();
          unsigned int solidSFDposition = 0;
          for (unsigned int k = 0; k < solidNodalSFDneighboursSize; k++)
          {
            if (idNodeOfConsideredElement == solidNodalSFDneighbours[k])
            {
              rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[solidSFDposition] += rDN_DX(j, 0) * elementVolume / solidNodalVolume;
              rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[solidSFDposition + 1] += rDN_DX(j, 1) * elementVolume / solidNodalVolume;
              if (TDim == 3)
              {
                rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[solidSFDposition + 2] += rDN_DX(j, 2) * elementVolume / solidNodalVolume;
              }
              break;
            }
            solidSFDposition += TDim;
          }
        }
      }
      else
      {
        std::cout << rGeom[i].Id() << "  this solid node is isolated!!! " << std::endl;
        for (unsigned int k = 0; k < TDim; k++)
        {
          rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_DEFORMATION_GRAD)(k, k) = 1.0;
        }
      }
      //}
    }
  }
  else
  //  if(this->Is(FLUID) && this->IsNot(SOLID))
  {
    for (unsigned int i = 0; i < NumNodes; i++)
    {
      VectorType nodalSFDneighbours = rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
      // std::cout<<"rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER) "<<rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER)<<std::endl;
      unsigned int nodalSFDneighboursSize = nodalSFDneighbours.size();

      if (nodalSFDneighboursSize > 1)
      {
        double &meshSize = rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
        ElementWeakPtrVectorType &neighb_elems = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
        double numberOfNeighElems = double(neighb_elems.size());
        meshSize += meanElementEdgesLength / numberOfNeighElems;

        if (rGeom[i].Is(FREE_SURFACE))
        {
          this->NodalFreeSurfaceLength(i);
        }

        double nodalVolume = rGeom[i].FastGetSolutionStepValue(NODAL_VOLUME);

        for (unsigned int j = 0; j < NumNodes; j++)
        {
          unsigned int idNodeOfConsideredElement = rGeom[j].Id();
          unsigned int SFDposition = 0;
          for (unsigned int k = 0; k < nodalSFDneighboursSize; k++)
          {
            if (idNodeOfConsideredElement == nodalSFDneighbours[k])
            {
              rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition] += rDN_DX(j, 0) * elementVolume / nodalVolume;
              rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition + 1] += rDN_DX(j, 1) * elementVolume / nodalVolume;
              if (TDim == 3)
              {
                rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition + 2] += rDN_DX(j, 2) * elementVolume / nodalVolume;
              }
              break;
            }
            SFDposition += TDim;
          }
        }
      }
      else
      {
        std::cout << rGeom[i].Id() << "  this fluid node is isolated!!! " << std::endl;
        for (unsigned int k = 0; k < TDim; k++)
        {
          rGeom[i].FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD)(k, k) = 1.0;
        }
      }
    }
  }

  // for (unsigned int i = 0; i < NumNodes; i++)
  // 	{

  //     if(rGeom[i].FastGetSolutionStepValue(INTERFACE_NODE)==true && this->Is(SOLID))
  //     // if(this->Is(SOLID))
  //     {
  //       VectorType solidNodalSFDneighbours=rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
  //       unsigned int solidNodalSFDneighboursSize=solidNodalSFDneighbours.size();
  //       // std::cout<<"SOLID_NODAL_SFD_NEIGHBOURS_ORDER) "<<rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER)<<std::endl;

  //       if(solidNodalSFDneighboursSize>1)
  //       {
  //         double& solidMeshSize = rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_MEAN_MESH_SIZE);
  //         ElementWeakPtrVectorType& neighb_elems = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
  //         double numberOfNeighElems=double(neighb_elems.size());
  //         solidMeshSize+=meanElementEdgesLength/numberOfNeighElems;

  //         double solidNodalVolume=rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_VOLUME);

  //         for (unsigned int j = 0; j< NumNodes; j++)
  //         {
  //             unsigned int idNodeOfConsideredElement=rGeom[j].Id();
  //             unsigned int solidSFDposition=0;
  //             for (unsigned int k = 0; k< solidNodalSFDneighboursSize; k++)
  //             {
  //               if(idNodeOfConsideredElement==solidNodalSFDneighbours[k])
  //               {
  //                 rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[solidSFDposition]   += rDN_DX(j,0)*elementVolume/solidNodalVolume;
  //                 rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[solidSFDposition+1] += rDN_DX(j,1)*elementVolume/solidNodalVolume;
  //                 if(TDim==3){
  //                   rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[solidSFDposition+2] += rDN_DX(j,2)*elementVolume/solidNodalVolume;
  //                 }
  //                 break;
  //               }
  //               solidSFDposition+=TDim;
  //           }
  //         }
  //       }
  //       else
  //       {
  //           std::cout<<rGeom[i].Id()<<"  this solid node is isolated!!! "<<std::endl;
  //           for (unsigned int k = 0; k< TDim; k++)
  //           {
  //             rGeom[i].FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD)(k,k)=1.0;
  //           }
  //       }
  //     }
  //     else if(this->Is(FLUID)){

  //       VectorType nodalSFDneighbours=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
  //       // std::cout<<"rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER) "<<rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER)<<std::endl;
  //       unsigned int nodalSFDneighboursSize=nodalSFDneighbours.size();

  //       if(nodalSFDneighboursSize>1)
  //       {
  //         double& meshSize = rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
  //         ElementWeakPtrVectorType& neighb_elems = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
  //         double numberOfNeighElems=double(neighb_elems.size());
  //         meshSize+=meanElementEdgesLength/numberOfNeighElems;

  //         if(rGeom[i].Is(FREE_SURFACE)){
  //           this->NodalFreeSurfaceLength(i);
  //         }

  //       double nodalVolume=rGeom[i].FastGetSolutionStepValue(NODAL_VOLUME);

  //       for (unsigned int j = 0; j< NumNodes; j++)
  //       {
  //           unsigned int idNodeOfConsideredElement=rGeom[j].Id();
  //           unsigned int SFDposition=0;
  //           for (unsigned int k = 0; k< nodalSFDneighboursSize; k++)
  //           {
  //             if(idNodeOfConsideredElement==nodalSFDneighbours[k])
  //             {
  //               rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition]   += rDN_DX(j,0)*elementVolume/nodalVolume;
  //               rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+1] += rDN_DX(j,1)*elementVolume/nodalVolume;
  //               if(TDim==3){
  //                 rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+2] += rDN_DX(j,2)*elementVolume/nodalVolume;
  //               }
  //               break;
  //             }
  //             SFDposition+=TDim;
  //         }
  //       }
  //     }
  //     else
  //     {
  //         std::cout<<rGeom[i].Id()<<"  this fluid node is isolated!!! "<<std::endl;
  //         for (unsigned int k = 0; k< TDim; k++)
  //         {
  //           rGeom[i].FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD)(k,k)=1.0;
  //         }
  //     }
  //     // if(rGeom[i].FastGetSolutionStepValue(INTERFACE_NODE)==true)
  //     // {

  //     //   if(this->Is(SOLID)){
  //     //     std::cout<<"        SOLID_NODAL_SFD_NEIGHBOURS_ORDER "<<rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER)<<std::endl;
  //     //     std::cout<<"        SOLID_NODAL_SFD_NEIGHBOURS "<<rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)<<std::endl;
  //     //   }else{
  //     //     std::cout<<"NODAL_SFD_NEIGHBOURS_ORDER "<<rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER)<<std::endl;
  //     //     std::cout<<"NODAL_SFD_NEIGHBOURS "<<rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)<<std::endl;
  //     //   }

  //     // }
  //   }

  // }

  // for (unsigned int i = 0; i < NumNodes; i++)
  // 	{

  //     VectorType nodalSFDneighbours=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
  //     unsigned int nodalSFDneighboursSize=nodalSFDneighbours.size();

  //     if(nodalSFDneighboursSize>1)
  //     {
  //       double& meshSize = rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
  //       ElementWeakPtrVectorType& neighb_elems = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
  //       double numberOfNeighElems=double(neighb_elems.size());
  //       meshSize+=meanElementEdgesLength/numberOfNeighElems;

  //       if(rGeom[i].Is(FREE_SURFACE)){
  //         this->NodalFreeSurfaceLength(i);
  //       }

  //     if(rGeom[i].FastGetSolutionStepValue(INTERFACE_NODE)==true && this->Is(SOLID))
  //     {
  //       VectorType solidNodalSFDneighbours=rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
  //       unsigned int solidNodalSFDneighboursSize=solidNodalSFDneighbours.size();
  //       double& solidMeshSize = rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_MEAN_MESH_SIZE);
  //       solidMeshSize=meshSize;

  //       double solidNodalVolume=rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_VOLUME);

  //       for (unsigned int j = 0; j< NumNodes; j++)
  //       {
  //           unsigned int idNodeOfConsideredElement=rGeom[j].Id();
  //           unsigned int solidSFDposition=0;
  //           for (unsigned int k = 0; k< solidNodalSFDneighboursSize; k++)
  //           {
  //             if(idNodeOfConsideredElement==solidNodalSFDneighbours[k])
  //             {
  //               rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[solidSFDposition]   += rDN_DX(j,0)*elementVolume/solidNodalVolume;
  //               rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[solidSFDposition+1] += rDN_DX(j,1)*elementVolume/solidNodalVolume;
  //               if(TDim==3){
  //                 rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[solidSFDposition+2] += rDN_DX(j,2)*elementVolume/solidNodalVolume;
  //               }
  //               break;
  //             }
  //             solidSFDposition+=TDim;
  //          }
  //       }
  //     }else{

  //       double nodalVolume=rGeom[i].FastGetSolutionStepValue(NODAL_VOLUME);

  //       for (unsigned int j = 0; j< NumNodes; j++)
  //       {
  //           unsigned int idNodeOfConsideredElement=rGeom[j].Id();
  //           unsigned int SFDposition=0;
  //           for (unsigned int k = 0; k< nodalSFDneighboursSize; k++)
  //           {
  //             if(idNodeOfConsideredElement==nodalSFDneighbours[k])
  //             {
  //               rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition]   += rDN_DX(j,0)*elementVolume/nodalVolume;
  //               rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+1] += rDN_DX(j,1)*elementVolume/nodalVolume;
  //               if(TDim==3){
  //                 rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+2] += rDN_DX(j,2)*elementVolume/nodalVolume;
  //               }
  //               break;
  //             }
  //             SFDposition+=TDim;
  //         }
  //       }
  //     }
  //     // if(rGeom[i].FastGetSolutionStepValue(INTERFACE_NODE)==true)
  //     // {

  //     //   if(this->Is(SOLID)){
  //     //     std::cout<<"        SOLID_NODAL_SFD_NEIGHBOURS_ORDER "<<rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER)<<std::endl;
  //     //     std::cout<<"        SOLID_NODAL_SFD_NEIGHBOURS "<<rGeom[i].FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)<<std::endl;
  //     //   }else{
  //     //     std::cout<<"NODAL_SFD_NEIGHBOURS_ORDER "<<rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER)<<std::endl;
  //     //     std::cout<<"NODAL_SFD_NEIGHBOURS "<<rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)<<std::endl;
  //     //   }

  //     // }
  //   }
  //   else
  //   {
  //       std::cout<<rGeom[i].Id()<<"  this node is isolated!!! "<<std::endl;
  //       for (unsigned int k = 0; k< TDim; k++)
  //       {
  //       	rGeom[i].FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD)(k,k)=1.0;
  //       }
  //   }
  // }

  KRATOS_CATCH("");
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::NodalFreeSurfaceLength(unsigned int nodeIndex)
{
  GeometryType &rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  array_1d<double, 3> Edge(3, 0.0);
  // unsigned int countFreeSurface=1;
  for (SizeType i = 0; i < NumNodes; i++)
  {

    if ((rGeom[i].Is(FREE_SURFACE) || (rGeom[i].Is(SOLID) && rGeom[i].Is(BOUNDARY))) && i != nodeIndex)
    {
      noalias(Edge) = rGeom[nodeIndex].Coordinates() - rGeom[i].Coordinates();
      rGeom[nodeIndex].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += sqrt(Edge[0] * Edge[0] + Edge[1] * Edge[1]) / 2.0;
      // countFreeSurface+=1;
    }
  }
  // if(countFreeSurface==NumNodes){
  //   ElementWeakPtrVectorType& neighb_elemsA = rGeom[0].GetValue(NEIGHBOUR_ELEMENTS);
  //   ElementWeakPtrVectorType& neighb_elemsB = rGeom[1].GetValue(NEIGHBOUR_ELEMENTS);
  //   ElementWeakPtrVectorType& neighb_elemsC = rGeom[2].GetValue(NEIGHBOUR_ELEMENTS);
  //   if(neighb_elemsA.size()==1 && neighb_elemsB.size()==1 && neighb_elemsC.size()==1){
  // 	rGeom[nodeIndex].Set(ISOLATED);
  //   }
  // }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::NodalFreeSurfaceLength(unsigned int nodeIndex)
{

  GeometryType &rGeom = this->GetGeometry();
  // const SizeType NumNodes = rGeom.PointsNumber();
  array_1d<double, 3> Edge(3, 0.0);

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {

    const double a = MathUtils<double>::Norm3(rGeom.GetPoint(0) - rGeom.GetPoint(1));
    const double b = MathUtils<double>::Norm3(rGeom.GetPoint(1) - rGeom.GetPoint(2));
    const double c = MathUtils<double>::Norm3(rGeom.GetPoint(2) - rGeom.GetPoint(0));

    const double s = (a + b + c) / 2.0;

    rGeom[nodeIndex].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += std::sqrt(s * (s - a) * (s - b) * (s - c)) / 3.0;
  }
  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    const double a = MathUtils<double>::Norm3(rGeom.GetPoint(0) - rGeom.GetPoint(1));
    const double b = MathUtils<double>::Norm3(rGeom.GetPoint(1) - rGeom.GetPoint(3));
    const double c = MathUtils<double>::Norm3(rGeom.GetPoint(3) - rGeom.GetPoint(0));

    const double s = (a + b + c) / 2.0;

    rGeom[nodeIndex].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += std::sqrt(s * (s - a) * (s - b) * (s - c)) / 3.0;
  }
  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    const double a = MathUtils<double>::Norm3(rGeom.GetPoint(0) - rGeom.GetPoint(2));
    const double b = MathUtils<double>::Norm3(rGeom.GetPoint(2) - rGeom.GetPoint(3));
    const double c = MathUtils<double>::Norm3(rGeom.GetPoint(3) - rGeom.GetPoint(0));

    const double s = (a + b + c) / 2.0;

    rGeom[nodeIndex].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += std::sqrt(s * (s - a) * (s - b) * (s - c)) / 3.0;
  }
  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    const double a = MathUtils<double>::Norm3(rGeom.GetPoint(1) - rGeom.GetPoint(2));
    const double b = MathUtils<double>::Norm3(rGeom.GetPoint(2) - rGeom.GetPoint(3));
    const double c = MathUtils<double>::Norm3(rGeom.GetPoint(3) - rGeom.GetPoint(1));

    const double s = (a + b + c) / 2.0;

    rGeom[nodeIndex].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += std::sqrt(s * (s - a) * (s - b) * (s - c)) / 3.0;
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::GetNodesPosition(Vector &rValues, const ProcessInfo &rCurrentProcessInfo, double theta)
{
  GeometryType &rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = 2 * NumNodes;

  if (rValues.size() != LocalSize)
    rValues.resize(LocalSize);

  SizeType Index = 0;

  for (SizeType i = 0; i < NumNodes; ++i)
  {
    rValues[Index++] = rGeom[i].X();
    rValues[Index++] = rGeom[i].Y();
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::GetNodesPosition(Vector &rValues, const ProcessInfo &rCurrentProcessInfo, double theta)
{
  GeometryType &rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = 3 * NumNodes;

  if (rValues.size() != LocalSize)
    rValues.resize(LocalSize);

  SizeType Index = 0;

  for (SizeType i = 0; i < NumNodes; ++i)
  {
    rValues[Index++] = rGeom[i].X();
    rValues[Index++] = rGeom[i].Y();
    rValues[Index++] = rGeom[i].Z();
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::CalcElementalStrains(ElementalVariables &rElementalVariables,
                                                                                         const ProcessInfo &rCurrentProcessInfo,
                                                                                         const ShapeFunctionDerivativesType &rDN_DX)
{
  unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
  GeometryType &rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = dimension * NumNodes;
  VectorType NodesPosition = ZeroVector(LocalSize);
  VectorType VelocityValues = ZeroVector(LocalSize);
  VectorType RHSVelocities = ZeroVector(LocalSize);
  double theta = 0.5;
  // if(rGeom[0].Is(SOLID) && rGeom[1].Is(SOLID) && rGeom[2].Is(SOLID)){
  //   theta=1.0;
  // }
  this->GetNodesPosition(NodesPosition, rCurrentProcessInfo, theta);
  this->GetVelocityValues(RHSVelocities, 0);
  RHSVelocities *= theta;
  this->GetVelocityValues(VelocityValues, 1);
  RHSVelocities += VelocityValues * (1.0 - theta);
  rElementalVariables.Fgrad = ZeroMatrix(dimension, dimension);
  rElementalVariables.FgradVel = ZeroMatrix(dimension, dimension);
  for (SizeType i = 0; i < dimension; i++)
  {
    for (SizeType j = 0; j < dimension; j++)
    {
      for (SizeType k = 0; k < NumNodes; k++)
      {
        rElementalVariables.Fgrad(i, j) += NodesPosition[dimension * k + i] * rDN_DX(k, j);
        rElementalVariables.FgradVel(i, j) += RHSVelocities[dimension * k + i] * rDN_DX(k, j);
      }
    }
  }

  //Inverse
  rElementalVariables.InvFgrad = ZeroMatrix(dimension, dimension);
  rElementalVariables.DetFgrad = 1;
  MathUtils<double>::InvertMatrix2(rElementalVariables.Fgrad,
                                   rElementalVariables.InvFgrad,
                                   rElementalVariables.DetFgrad);

  //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
  rElementalVariables.SpatialVelocityGrad.resize(dimension, dimension);
  rElementalVariables.SpatialVelocityGrad = prod(rElementalVariables.FgradVel, rElementalVariables.InvFgrad);

  rElementalVariables.VolumetricDefRate = 0;
  for (SizeType i = 0; i < dimension; i++)
  {
    rElementalVariables.VolumetricDefRate += rElementalVariables.SpatialVelocityGrad(i, i);
  }

  rElementalVariables.SpatialDefRate[0] = rElementalVariables.SpatialVelocityGrad(0, 0);
  rElementalVariables.SpatialDefRate[1] = rElementalVariables.SpatialVelocityGrad(1, 1);
  rElementalVariables.SpatialDefRate[2] = 0.5 * (rElementalVariables.SpatialVelocityGrad(1, 0) + rElementalVariables.SpatialVelocityGrad(0, 1));

  double aThird = 1.0 / 3.0;
  double dev_X = rElementalVariables.SpatialDefRate[0] -
                 (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1]) * aThird;
  double dev_Y = rElementalVariables.SpatialDefRate[1] -
                 (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1]) * aThird;
  rElementalVariables.DeviatoricInvariant = sqrt(2 * (dev_X * dev_X + dev_Y * dev_Y +
                                                      rElementalVariables.SpatialDefRate[2] * rElementalVariables.SpatialDefRate[2]));

  rElementalVariables.EquivalentStrainRate = sqrt((2.0 * rElementalVariables.SpatialDefRate[0] * rElementalVariables.SpatialDefRate[0] +
                                                   2.0 * rElementalVariables.SpatialDefRate[1] * rElementalVariables.SpatialDefRate[1] +
                                                   4.0 * rElementalVariables.SpatialDefRate[2] * rElementalVariables.SpatialDefRate[2]));
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::CalcElementalStrains(ElementalVariables &rElementalVariables,
                                                                                         const ProcessInfo &rCurrentProcessInfo,
                                                                                         const ShapeFunctionDerivativesType &rDN_DX)
{
  unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
  GeometryType &rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = dimension * NumNodes;
  VectorType NodesPosition = ZeroVector(LocalSize);
  VectorType VelocityValues = ZeroVector(LocalSize);
  VectorType RHSVelocities = ZeroVector(LocalSize);
  double theta = 0.5;
  if (rGeom[0].Is(SOLID) && rGeom[1].Is(SOLID) && rGeom[2].Is(SOLID) && rGeom[3].Is(SOLID))
  {
    theta = 1.0;
  }
  this->GetNodesPosition(NodesPosition, rCurrentProcessInfo, theta);
  this->GetVelocityValues(RHSVelocities, 0);
  RHSVelocities *= theta;
  this->GetVelocityValues(VelocityValues, 1);
  RHSVelocities += VelocityValues * (1.0 - theta);

  rElementalVariables.Fgrad = ZeroMatrix(dimension, dimension);
  rElementalVariables.FgradVel = ZeroMatrix(dimension, dimension);
  for (SizeType i = 0; i < dimension; i++)
  {
    for (SizeType j = 0; j < dimension; j++)
    {
      for (SizeType k = 0; k < NumNodes; k++)
      {
        rElementalVariables.Fgrad(i, j) += NodesPosition[dimension * k + i] * rDN_DX(k, j);
        rElementalVariables.FgradVel(i, j) += RHSVelocities[dimension * k + i] * rDN_DX(k, j);
      }
    }
  }

  //Inverse
  rElementalVariables.InvFgrad = ZeroMatrix(dimension, dimension);
  rElementalVariables.DetFgrad = 1;
  MathUtils<double>::InvertMatrix3(rElementalVariables.Fgrad,
                                   rElementalVariables.InvFgrad,
                                   rElementalVariables.DetFgrad);

  //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
  rElementalVariables.SpatialVelocityGrad.resize(dimension, dimension);
  rElementalVariables.SpatialVelocityGrad = prod(rElementalVariables.FgradVel, rElementalVariables.InvFgrad);

  rElementalVariables.VolumetricDefRate = 0;
  for (SizeType i = 0; i < dimension; i++)
  {
    rElementalVariables.VolumetricDefRate += rElementalVariables.SpatialVelocityGrad(i, i);
  }

  rElementalVariables.SpatialDefRate[0] = rElementalVariables.SpatialVelocityGrad(0, 0);
  rElementalVariables.SpatialDefRate[1] = rElementalVariables.SpatialVelocityGrad(1, 1);
  rElementalVariables.SpatialDefRate[2] = rElementalVariables.SpatialVelocityGrad(2, 2);
  rElementalVariables.SpatialDefRate[3] = 0.5 * (rElementalVariables.SpatialVelocityGrad(1, 0) + rElementalVariables.SpatialVelocityGrad(0, 1));
  rElementalVariables.SpatialDefRate[4] = 0.5 * (rElementalVariables.SpatialVelocityGrad(2, 0) + rElementalVariables.SpatialVelocityGrad(0, 2));
  rElementalVariables.SpatialDefRate[5] = 0.5 * (rElementalVariables.SpatialVelocityGrad(2, 1) + rElementalVariables.SpatialVelocityGrad(1, 2));
  // computeElement=CheckStrain3(rElementalVariables.SpatialDefRate,rElementalVariables.SpatialVelocityGrad);

  double aThird = 1.0 / 3.0;
  double dev_X = rElementalVariables.SpatialDefRate[0] -
                 (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1] + rElementalVariables.SpatialDefRate[2]) * aThird;
  double dev_Y = rElementalVariables.SpatialDefRate[1] -
                 (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1] + rElementalVariables.SpatialDefRate[2]) * aThird;
  double dev_Z = rElementalVariables.SpatialDefRate[2] -
                 (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1] + rElementalVariables.SpatialDefRate[2]) * aThird;
  rElementalVariables.DeviatoricInvariant = sqrt(2 * (dev_X * dev_X + dev_Y * dev_Y + dev_Z * dev_Z +
                                                      rElementalVariables.SpatialDefRate[3] * rElementalVariables.SpatialDefRate[3] +
                                                      rElementalVariables.SpatialDefRate[4] * rElementalVariables.SpatialDefRate[4] +
                                                      rElementalVariables.SpatialDefRate[5] * rElementalVariables.SpatialDefRate[5]));

  rElementalVariables.EquivalentStrainRate = sqrt(2.0 * (rElementalVariables.SpatialDefRate[0] * rElementalVariables.SpatialDefRate[0] +
                                                         rElementalVariables.SpatialDefRate[1] * rElementalVariables.SpatialDefRate[1] +
                                                         rElementalVariables.SpatialDefRate[2] * rElementalVariables.SpatialDefRate[2] +
                                                         2.0 * rElementalVariables.SpatialDefRate[3] * rElementalVariables.SpatialDefRate[3] +
                                                         2.0 * rElementalVariables.SpatialDefRate[4] * rElementalVariables.SpatialDefRate[4] +
                                                         2.0 * rElementalVariables.SpatialDefRate[5] * rElementalVariables.SpatialDefRate[5]));
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateGeometryData(ShapeFunctionDerivativesArrayType &rDN_DX,
                                                                                             Matrix &NContainer,
                                                                                             Vector &rGaussWeights)
{
  const GeometryType &rGeom = this->GetGeometry();
  Vector DetJ;
  rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, GeometryData::GI_GAUSS_1);
  NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
  const GeometryType::IntegrationPointsArrayType &IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1);

  rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1), false);

  for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1); ++g)
  {
    // rGaussWeights[g] = fabs(DetJ[g] * IntegrationPoints[g].Weight());
    rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
    // if(rGaussWeights[g]<0)
    // 	std::cout<<"NEGATIVE GAUSS WEIGHT "<<rGaussWeights[g]<<std::endl;
  }
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::GetValueOnIntegrationPoints(const Variable<double> &rVariable,
                                                                                                   std::vector<double> &rValues,
                                                                                                   const ProcessInfo &rCurrentProcessInfo)
{
  if (rVariable == YIELDED)
  {
    rValues[0] = this->GetValue(YIELDED);
  }
  if (rVariable == FLOW_INDEX)
  {
    rValues[0] = this->GetValue(FLOW_INDEX);
  }
}

template <unsigned int TDim>
int TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
{
  KRATOS_TRY;

  // Base class checks for positive Jacobian and Id > 0
  int ierr = Element::Check(rCurrentProcessInfo);
  if (ierr != 0)
    return ierr;

  // Check that all required variables have been registered
  if (VELOCITY.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "VELOCITY Key is 0. Check that the application was correctly registered.", "");
  if (ACCELERATION.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "ACCELERATION Key is 0. Check that the application was correctly registered.", "");
  if (PRESSURE.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "PRESSURE Key is 0. Check that the application was correctly registered.", "");
  if (BODY_FORCE.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "BODY_FORCE Key is 0. Check that the application was correctly registered.", "");
  if (DENSITY.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "DENSITY Key is 0. Check that the application was correctly registered.", "");
  if (DYNAMIC_VISCOSITY.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "DYNAMIC_VISCOSITY Key is 0. Check that the application was correctly registered.", "");
  if (DELTA_TIME.Key() == 0)
    KRATOS_THROW_ERROR(std::invalid_argument, "DELTA_TIME Key is 0. Check that the application was correctly registered.", "");

  // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
  for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
  {
    if (this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY variable on solution step data for node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE variable on solution step data for node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing BODY_FORCE variable on solution step data for node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing DENSITY variable on solution step data for node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].SolutionStepsDataHas(DYNAMIC_VISCOSITY) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing DYNAMIC_VISCOSITY variable on solution step data for node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
        this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
        this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY component degree of freedom on node ", this->GetGeometry()[i].Id());
    if (this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
      KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE component degree of freedom on node ", this->GetGeometry()[i].Id());
  }

  // If this is a 2D problem, check that nodes are in XY plane
  if (this->GetGeometry().WorkingSpaceDimension() == 2)
  {
    for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
    {
      if (this->GetGeometry()[i].Z() != 0.0)
        KRATOS_THROW_ERROR(std::invalid_argument, "Node with non-zero Z coordinate found. Id: ", this->GetGeometry()[i].Id());
    }
  }

  return ierr;

  KRATOS_CATCH("");
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateElementalLaplacian(MatrixType &rLeftHandSideMatrix,
                                                                                                   VectorType &rRightHandSideVector,
                                                                                                   ProcessInfo &rCurrentProcessInfo)
{

  GeometryType &rGeom = this->GetGeometry();
  const unsigned int NumNodes = rGeom.PointsNumber();

  // Check sizes and initialize
  if (rLeftHandSideMatrix.size1() != NumNodes)
    rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

  rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

  if (rRightHandSideVector.size() != NumNodes)
    rRightHandSideVector.resize(NumNodes);

  rRightHandSideVector = ZeroVector(NumNodes);

  ShapeFunctionDerivativesArrayType DN_DX;
  Matrix NContainer;
  VectorType GaussWeights;
  this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
  const unsigned int NumGauss = GaussWeights.size();
  double Tau = rGeom[0].FastGetSolutionStepValue(NODAL_TAU);

  for (unsigned int i = 1; i < NumNodes; i++)
  {
    Tau += rGeom[i].FastGetSolutionStepValue(NODAL_TAU);
  }

  Tau *= (1.0 / double(NumNodes));

  for (unsigned int g = 0; g < NumGauss; ++g)
  {
    const double GaussWeight = GaussWeights[g];
    const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
    double StabLaplacianWeight = Tau * GaussWeight;
    this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix, rDN_DX, StabLaplacianWeight);
  }

  VectorType PressureValuesForRHS = ZeroVector(NumNodes);
  this->GetPressureValues(PressureValuesForRHS, 0);
  noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateElementalLaplacianAndTau(MatrixType &rLeftHandSideMatrix,
                                                                                                         VectorType &rRightHandSideVector,
                                                                                                         ProcessInfo &rCurrentProcessInfo)
{

  GeometryType &rGeom = this->GetGeometry();
  const unsigned int NumNodes = rGeom.PointsNumber();

  // Check sizes and initialize
  if (rLeftHandSideMatrix.size1() != NumNodes)
    rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

  rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

  if (rRightHandSideVector.size() != NumNodes)
    rRightHandSideVector.resize(NumNodes);

  rRightHandSideVector = ZeroVector(NumNodes);

  ShapeFunctionDerivativesArrayType DN_DX;
  Matrix NContainer;
  VectorType GaussWeights;
  this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
  const unsigned int NumGauss = GaussWeights.size();
  // double Tau=rGeom[0].FastGetSolutionStepValue(NODAL_TAU);
  double factor = 1.0 / double(NumNodes);
  double ElemSize = rGeom[0].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
  double Density = rGeom[0].FastGetSolutionStepValue(DENSITY);
  double DeviatoricCoeff = rGeom[0].FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);

  for (unsigned int i = 1; i < NumNodes; i++)
  {
    ElemSize += rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
    Density += rGeom[i].FastGetSolutionStepValue(DENSITY);
    DeviatoricCoeff += rGeom[i].FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);
    // Tau += rGeom[i].FastGetSolutionStepValue(NODAL_TAU);
  }

  ElemSize *= factor;
  Density *= factor;
  DeviatoricCoeff *= factor;
  // Tau*=(1.0/double(NumNodes));

  double maxViscousValueForStabilization = 0.1;
  if (DeviatoricCoeff > maxViscousValueForStabilization)
  {
    DeviatoricCoeff = maxViscousValueForStabilization;
  }

  double Tau = 0;
  this->CalculateTauFIC(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);
  for (unsigned int g = 0; g < NumGauss; ++g)
  {
    const double GaussWeight = GaussWeights[g];
    const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
    double StabLaplacianWeight = Tau * GaussWeight;
    this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix, rDN_DX, StabLaplacianWeight);
  }

  VectorType PressureValuesForRHS = ZeroVector(NumNodes);
  this->GetPressureValues(PressureValuesForRHS, 0);
  noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateVolumetricStabilizedTerms(MatrixType &rLeftHandSideMatrix,
                                                                                                          VectorType &rRightHandSideVector,
                                                                                                          ProcessInfo &rCurrentProcessInfo)
{

  GeometryType &rGeom = this->GetGeometry();
  const unsigned int NumNodes = rGeom.PointsNumber();

  // Check sizes and initialize
  if (rLeftHandSideMatrix.size1() != NumNodes)
    rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

  rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

  if (rRightHandSideVector.size() != NumNodes)
    rRightHandSideVector.resize(NumNodes);

  rRightHandSideVector = ZeroVector(NumNodes);

  ShapeFunctionDerivativesArrayType DN_DX;
  Matrix NContainer;
  VectorType GaussWeights;
  this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
  const unsigned int NumGauss = GaussWeights.size();
  // double Tau=rGeom[0].FastGetSolutionStepValue(NODAL_TAU);
  double factor = 1.0 / double(NumNodes);
  double ElemSize = rGeom[0].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
  double Density = rGeom[0].FastGetSolutionStepValue(DENSITY);
  double DeviatoricCoeff = rGeom[0].FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);

  for (unsigned int i = 1; i < NumNodes; i++)
  {
    ElemSize += rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
    Density += rGeom[i].FastGetSolutionStepValue(DENSITY);
    DeviatoricCoeff += rGeom[i].FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);
    // Tau += rGeom[i].FastGetSolutionStepValue(NODAL_TAU);
  }

  ElemSize *= factor;
  Density *= factor;
  DeviatoricCoeff *= factor;
  //  Tau*=(1.0/double(NumNodes));

  double maxViscousValueForStabilization = 0.1;
  if (DeviatoricCoeff > maxViscousValueForStabilization)
  {
    DeviatoricCoeff = maxViscousValueForStabilization;
  }

  double Tau = 0;
  this->CalculateTauFIC(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);
  for (unsigned int g = 0; g < NumGauss; ++g)
  {
    const double GaussWeight = GaussWeights[g];
    const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
    double StabLaplacianWeight = Tau * GaussWeight;
    this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix, rDN_DX, StabLaplacianWeight);
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      this->AddStabilizationNodalTermsRHS(rRightHandSideVector, Tau, Density, GaussWeight, rDN_DX, i);
    }
  }

  VectorType PressureValuesForRHS = ZeroVector(NumNodes);
  this->GetPressureValues(PressureValuesForRHS, 0);
  noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateElementalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
                                                                                                                 VectorType &rRightHandSideVector,
                                                                                                                 ProcessInfo &rCurrentProcessInfo)
{
  GeometryType &rGeom = this->GetGeometry();
  const unsigned int NumNodes = rGeom.PointsNumber();

  // Check sizes and initialize
  if (rLeftHandSideMatrix.size1() != NumNodes)
    rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

  rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

  if (rRightHandSideVector.size() != NumNodes)
    rRightHandSideVector.resize(NumNodes);

  rRightHandSideVector = ZeroVector(NumNodes);

  // Shape functions and integration points
  ShapeFunctionDerivativesArrayType DN_DX;
  Matrix NContainer;
  VectorType GaussWeights;
  this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
  const unsigned int NumGauss = GaussWeights.size();

  double TimeStep = rCurrentProcessInfo[DELTA_TIME];

  double theta = 1.0;
  ElementalVariables rElementalVariables;
  this->InitializeElementalVariables(rElementalVariables);

  double factor = 1.0 / double(NumNodes);
  array_1d<double, 3> Normal = rGeom[0].FastGetSolutionStepValue(NORMAL);
  Vector SpatialDefRate = rGeom[0].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);

  double ElemSize = rGeom[0].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
  double Density = rGeom[0].FastGetSolutionStepValue(DENSITY);
  double DeviatoricCoeff = rGeom[0].FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);
  double VolumetricCoeff = TimeStep * rGeom[0].FastGetSolutionStepValue(BULK_MODULUS);
  double volumetricDefRate = rGeom[0].GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE);
  double elementalNormalProjDefRate = 0;

  if (TDim == 2)
  {
    elementalNormalProjDefRate += Normal[0] * SpatialDefRate[0] * Normal[0] + Normal[1] * SpatialDefRate[1] * Normal[1] + 2 * Normal[0] * SpatialDefRate[2] * Normal[1];
  }
  else if (TDim == 3)
  {
    elementalNormalProjDefRate += Normal[0] * SpatialDefRate[0] * Normal[0] + Normal[1] * SpatialDefRate[1] * Normal[1] + Normal[2] * SpatialDefRate[2] * Normal[2] +
                                  2 * Normal[0] * SpatialDefRate[3] * Normal[1] + 2 * Normal[0] * SpatialDefRate[4] * Normal[2] + 2 * Normal[1] * SpatialDefRate[5] * Normal[2];
  }
  for (unsigned int i = 1; i < NumNodes; i++)
  {
    Normal = rGeom[i].FastGetSolutionStepValue(NORMAL);
    SpatialDefRate = rGeom[i].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);

    ElemSize += rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
    Density += rGeom[i].FastGetSolutionStepValue(DENSITY);
    DeviatoricCoeff += rGeom[i].FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);
    VolumetricCoeff += TimeStep * rGeom[i].FastGetSolutionStepValue(BULK_MODULUS);
    volumetricDefRate += rGeom[i].GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE);

    if (TDim == 2)
    {
      elementalNormalProjDefRate += Normal[0] * SpatialDefRate[0] * Normal[0] + Normal[1] * SpatialDefRate[1] * Normal[1] + 2 * Normal[0] * SpatialDefRate[2] * Normal[1];
    }
    else if (TDim == 3)
    {
      elementalNormalProjDefRate += Normal[0] * SpatialDefRate[0] * Normal[0] + Normal[1] * SpatialDefRate[1] * Normal[1] + Normal[2] * SpatialDefRate[2] * Normal[2] +
                                    2 * Normal[0] * SpatialDefRate[3] * Normal[1] + 2 * Normal[0] * SpatialDefRate[4] * Normal[2] + 2 * Normal[1] * SpatialDefRate[5] * Normal[2];
    }
  }

  ElemSize *= factor;
  Density *= factor;
  DeviatoricCoeff *= factor;
  VolumetricCoeff *= factor;
  volumetricDefRate *= factor;
  elementalNormalProjDefRate *= factor;

  double maxViscousValueForStabilization = 0.1;
  if (DeviatoricCoeff > maxViscousValueForStabilization)
  {
    DeviatoricCoeff = maxViscousValueForStabilization;
  }

  double Tau = 0;
  this->CalculateTauFIC(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);

  double totalVolume = 0;
  bool computeElement = false;
  for (unsigned int g = 0; g < NumGauss; ++g)
  {
    const double GaussWeight = GaussWeights[g];
    totalVolume += GaussWeight;
    const ShapeFunctionsType &N = row(NContainer, g);
    const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
    computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);

    double BoundLHSCoeff = Tau * 4.0 * GaussWeight / (ElemSize * ElemSize);
    this->ComputeBoundLHSMatrix(rLeftHandSideMatrix, N, BoundLHSCoeff);

    double BoundRHSCoeffAcc = Tau * Density * 2 * GaussWeight / ElemSize;
    // double BoundRHSCoeffDev=elementalNormalProjDefRate*Tau*8.0*DeviatoricCoeff*GaussWeight/(ElemSize*ElemSize);
    double BoundRHSCoeffDev = Tau * 8.0 * DeviatoricCoeff * GaussWeight / (ElemSize * ElemSize);

    // this->ComputeElementalBoundRHSVector(rRightHandSideVector,TimeStep,BoundRHSCoeffAcc,BoundRHSCoeffDev);
    this->ComputeBoundRHSVectorComplete(rRightHandSideVector, TimeStep, BoundRHSCoeffAcc, BoundRHSCoeffDev, rElementalVariables.SpatialDefRate);

    double StabLaplacianWeight = Tau * GaussWeight;
    this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix, rDN_DX, StabLaplacianWeight);
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      // rRightHandSideVector[i] += GaussWeight * N[i] * rGeom[i].GetSolutionStepValue(NODAL_VOLUMETRIC_DEF_RATE);
      rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;
      // rRightHandSideVector[i] += GaussWeight * N[i] * volumetricDefRate;
      this->AddStabilizationNodalTermsRHS(rRightHandSideVector, Tau, Density, GaussWeight, rDN_DX, i);
    }
  }

  computeElement = true;

  if (computeElement == true)
  {

    VectorType PressureValues = ZeroVector(NumNodes);
    VectorType PressureValuesForRHS = ZeroVector(NumNodes);
    this->GetPressureValues(PressureValuesForRHS, 0);
    //the LHS matrix up to now just contains the laplacian term and the bound term
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);

    this->GetPressureValues(PressureValues, 1);
    noalias(PressureValuesForRHS) += -PressureValues;
    MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
    MatrixType BulkMatrixConsistent = ZeroMatrix(NumNodes, NumNodes);
    double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);

    this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
    noalias(rLeftHandSideMatrix) += BulkMatrix;
    noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);

    double lumpedBulkStabCoeff = lumpedBulkCoeff * Tau * Density / TimeStep;
    this->GetPressureVelocityValues(PressureValues, 0);
    noalias(PressureValuesForRHS) += -PressureValues * TimeStep;
    noalias(BulkMatrix) = ZeroMatrix(NumNodes, NumNodes);
    this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkStabCoeff);
    noalias(rLeftHandSideMatrix) += BulkMatrix;
    noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
  }
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
                                                                                                            VectorType &rRightHandSideVector,
                                                                                                            ProcessInfo &rCurrentProcessInfo)
{

  GeometryType &rGeom = this->GetGeometry();
  const unsigned int NumNodes = rGeom.PointsNumber();

  // Check sizes and initialize
  if (rLeftHandSideMatrix.size1() != NumNodes)
    rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

  rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

  if (rRightHandSideVector.size() != NumNodes)
    rRightHandSideVector.resize(NumNodes);

  rRightHandSideVector = ZeroVector(NumNodes);

  ShapeFunctionDerivativesArrayType DN_DX;
  Matrix NContainer;
  VectorType GaussWeights;
  this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
  double TimeStep = rCurrentProcessInfo[DELTA_TIME];
  const unsigned int NumGauss = GaussWeights.size();
  // double Tau=rGeom[0].FastGetSolutionStepValue(NODAL_TAU);
  double factor = 1.0 / double(NumNodes);
  double ElemSize = rGeom[0].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
  double Density = rGeom[0].FastGetSolutionStepValue(DENSITY);
  double DeviatoricCoeff = rGeom[0].FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);
  for (unsigned int i = 1; i < NumNodes; i++)
  {
    ElemSize += rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
    Density += rGeom[i].FastGetSolutionStepValue(DENSITY);
    DeviatoricCoeff += rGeom[i].FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);
  }

  ElemSize *= factor * 0.5;
  Density *= factor;
  DeviatoricCoeff *= factor;
  //  Tau*=(1.0/double(NumNodes));

  double maxViscousValueForStabilization = 0.1;
  if (DeviatoricCoeff > maxViscousValueForStabilization)
  {
    DeviatoricCoeff = maxViscousValueForStabilization;
  }

  double Tau = 0;
  this->CalculateTauFIC(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);

  for (unsigned int g = 0; g < NumGauss; ++g)
  {
    const double GaussWeight = GaussWeights[g];
    const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

    // ElemSize*= 0.5;

    const ShapeFunctionsType &N = row(NContainer, g);
    double BoundLHSCoeff = Tau * 4.0 * GaussWeight / (ElemSize * ElemSize);

    this->ComputeBoundLHSMatrix(rLeftHandSideMatrix, N, BoundLHSCoeff);

    //double BoundRHSCoeffDev=elementalNormalProjDefRate*Tau*8.0*DeviatoricCoeff*GaussWeight/(ElemSize*ElemSize);
    double BoundRHSCoeffDev = Tau * 8.0 * DeviatoricCoeff * GaussWeight / (ElemSize * ElemSize);

    double BoundRHSCoeffAcc = Tau * Density * 2 * GaussWeight / ElemSize;
    this->ComputeElementalBoundRHSVector(rRightHandSideVector, TimeStep, BoundRHSCoeffAcc, BoundRHSCoeffDev);

    double StabLaplacianWeight = Tau * GaussWeight;
    this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix, rDN_DX, StabLaplacianWeight);
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      this->AddStabilizationNodalTermsRHS(rRightHandSideVector, Tau, Density, GaussWeight, rDN_DX, i);
    }
  }

  VectorType PressureValuesForRHS = ZeroVector(NumNodes);
  this->GetPressureValues(PressureValuesForRHS, 0);
  noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);
}



template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateStabilizingTermsContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
                                                                                                                        VectorType &rRightHandSideVector,
                                                                                                                        ProcessInfo &rCurrentProcessInfo)
{

  GeometryType &rGeom = this->GetGeometry();
  const unsigned int NumNodes = rGeom.PointsNumber();

  // Check sizes and initialize
  if (rLeftHandSideMatrix.size1() != NumNodes)
    rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

  rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

  if (rRightHandSideVector.size() != NumNodes)
    rRightHandSideVector.resize(NumNodes);

  rRightHandSideVector = ZeroVector(NumNodes);

  ShapeFunctionDerivativesArrayType DN_DX;
  Matrix NContainer;
  VectorType GaussWeights;
  this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
  double TimeStep = rCurrentProcessInfo[DELTA_TIME];
  const unsigned int NumGauss = GaussWeights.size();
  // double Tau=rGeom[0].FastGetSolutionStepValue(NODAL_TAU);
  double factor = 1.0 / double(NumNodes);
  double ElemSize = rGeom[0].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
  double Density = rGeom[0].FastGetSolutionStepValue(DENSITY);
  double DeviatoricCoeff = rGeom[0].FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);
  for (unsigned int i = 1; i < NumNodes; i++)
  {
    ElemSize += rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
    Density += rGeom[i].FastGetSolutionStepValue(DENSITY);
    DeviatoricCoeff += rGeom[i].FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);
  }

  ElemSize *= factor * 0.5;
  Density *= factor;
  DeviatoricCoeff *= factor;
  //  Tau*=(1.0/double(NumNodes));

  double maxViscousValueForStabilization = 0.1;
  if (DeviatoricCoeff > maxViscousValueForStabilization)
  {
    DeviatoricCoeff = maxViscousValueForStabilization;
  }

  double Tau = 0;
  this->CalculateTauFIC(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);

  for (unsigned int g = 0; g < NumGauss; ++g)
  {
    const double GaussWeight = GaussWeights[g];
    const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

    // ElemSize*= 0.5;

    const ShapeFunctionsType &N = row(NContainer, g);
    double BoundLHSCoeff = Tau * 4.0 * GaussWeight / (ElemSize * ElemSize);

    this->ComputeBoundLHSMatrix(rLeftHandSideMatrix, N, BoundLHSCoeff);

    //double BoundRHSCoeffDev=elementalNormalProjDefRate*Tau*8.0*DeviatoricCoeff*GaussWeight/(ElemSize*ElemSize);
    double BoundRHSCoeffDev = Tau * 8.0 * DeviatoricCoeff * GaussWeight / (ElemSize * ElemSize);

    double BoundRHSCoeffAcc = Tau * Density * 2 * GaussWeight / ElemSize;
    this->ComputeElementalBoundRHSVector(rRightHandSideVector, TimeStep, BoundRHSCoeffAcc, BoundRHSCoeffDev);

    double StabLaplacianWeight = Tau * GaussWeight;
    this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix, rDN_DX, StabLaplacianWeight);
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      this->AddStabilizationNodalTermsRHS(rRightHandSideVector, Tau, Density, GaussWeight, rDN_DX, i);
    }
  }

  VectorType PressureValuesForRHS = ZeroVector(NumNodes);
  this->GetPressureValues(PressureValuesForRHS, 0);
  noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateTauFIC(double &Tau,
                                                                                       double ElemSize,
                                                                                       const double Density,
                                                                                       const double Viscosity,
                                                                                       const ProcessInfo &rCurrentProcessInfo)
{
  double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
  if (rCurrentProcessInfo.GetValue(DELTA_TIME) < rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME))
  {
    DeltaTime = 0.5 * rCurrentProcessInfo.GetValue(DELTA_TIME) + 0.5 * rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME);
  }

  double MeanVelocity = 0;
  this->CalcMeanVelocity(MeanVelocity, 0);

  Tau = (ElemSize * ElemSize * DeltaTime) / (Density * MeanVelocity * DeltaTime * ElemSize + Density * ElemSize * ElemSize + 8.0 * Viscosity * DeltaTime);
  if (MeanVelocity == 0)
  {
    Tau = 0;
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::ComputeBoundLHSMatrix(Matrix &BoundLHSMatrix,
                                                                                          const ShapeFunctionsType &rN,
                                                                                          const double Weight)
{
  GeometryType &rGeom = this->GetGeometry();
  double coeff = 1.0 / 3.0;

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
  {
    if (rGeom[0].IsNot(INLET))
      BoundLHSMatrix(0, 0) += Weight * coeff;
    if (rGeom[1].IsNot(INLET))
      BoundLHSMatrix(1, 1) += Weight * coeff;
  }
  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {
    if (rGeom[0].IsNot(INLET))
      BoundLHSMatrix(0, 0) += Weight * coeff;
    if (rGeom[2].IsNot(INLET))
      BoundLHSMatrix(2, 2) += Weight * coeff;
  }
  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {
    if (rGeom[1].IsNot(INLET))
      BoundLHSMatrix(1, 1) += Weight * coeff;
    if (rGeom[2].IsNot(INLET))
      BoundLHSMatrix(2, 2) += Weight * coeff;
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::ComputeBoundLHSMatrix(Matrix &BoundLHSMatrix,
                                                                                          const ShapeFunctionsType &rN,
                                                                                          const double Weight)
{
  GeometryType &rGeom = this->GetGeometry();

  double coeff = 1.0 / 4.0;

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {
    if (rGeom[0].IsNot(INLET))
      BoundLHSMatrix(0, 0) += Weight * coeff;
    if (rGeom[1].IsNot(INLET))
      BoundLHSMatrix(1, 1) += Weight * coeff;
    if (rGeom[2].IsNot(INLET))
      BoundLHSMatrix(2, 2) += Weight * coeff;
  }
  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {
    if (rGeom[0].IsNot(INLET))
      BoundLHSMatrix(0, 0) += Weight * coeff;
    if (rGeom[1].IsNot(INLET))
      BoundLHSMatrix(1, 1) += Weight * coeff;
    if (rGeom[3].IsNot(INLET))
      BoundLHSMatrix(3, 3) += Weight * coeff;
  }
  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {
    if (rGeom[0].IsNot(INLET))
      BoundLHSMatrix(0, 0) += Weight * coeff;
    if (rGeom[2].IsNot(INLET))
      BoundLHSMatrix(2, 2) += Weight * coeff;
    if (rGeom[3].IsNot(INLET))
      BoundLHSMatrix(3, 3) += Weight * coeff;
  }
  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {
    if (rGeom[1].IsNot(INLET))
      BoundLHSMatrix(1, 1) += Weight * coeff;
    if (rGeom[2].IsNot(INLET))
      BoundLHSMatrix(2, 2) += Weight * coeff;
    if (rGeom[3].IsNot(INLET))
      BoundLHSMatrix(3, 3) += Weight * coeff;
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::ComputeElementalBoundRHSVector(VectorType &BoundRHSVector,
                                                                                                   const double TimeStep,
                                                                                                   const double BoundRHSCoeffAcc,
                                                                                                   const double BoundRHSCoeffDev)
{
  GeometryType &rGeom = this->GetGeometry();

  // Vector SpatialDefRate =  1.0/3.0*(rGeom[0].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE) + rGeom[1].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE) +  rGeom[2].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE));

  //   array_1d<double, 3> AccAa(3,0.0);
  //   array_1d<double, 3> AccBb(3,0.0);
  //   array_1d<double, 3> AccCc(3,0.0);
  //   array_1d<double, 3> MeanAcc(3,0.0);

  //   noalias(AccAa)= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);
  //   noalias(AccBb)= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);
  //   noalias(AccCc)= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);

  //   noalias(MeanAcc)= (AccAa + AccBb + AccCc)/3.0;

  //   const double coeff = 1.0/3.0;
  const double coeff = 1.0 / 3.0;
  const double timeFactor = 0.5 / TimeStep;

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
  {
    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);

    this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 1, 2);

    double elementalNormalProjDefRate = 0;
    Vector SpatialDefRate = rGeom[0].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += 0.5 * (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + 2 * NormalVector[0] * SpatialDefRate[2] * NormalVector[1]);
    SpatialDefRate = rGeom[1].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += 0.5 * (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + 2 * NormalVector[0] * SpatialDefRate[2] * NormalVector[1]);

    noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(MeanAcc) = 0.5 * AccA + 0.5 * AccB;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
      BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);
  }

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);

    this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 2, 1);

    double elementalNormalProjDefRate = 0;
    Vector SpatialDefRate = rGeom[0].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += 0.5 * (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + 2 * NormalVector[0] * SpatialDefRate[2] * NormalVector[1]);
    SpatialDefRate = rGeom[2].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += 0.5 * (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + 2 * NormalVector[0] * SpatialDefRate[2] * NormalVector[1]);

    noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(MeanAcc) = 0.5 * AccA + 0.5 * AccB;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
      BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);
  }

  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);

    this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 1, 2, 0);

    double elementalNormalProjDefRate = 0;
    Vector SpatialDefRate = rGeom[1].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += 0.5 * (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + 2 * NormalVector[0] * SpatialDefRate[2] * NormalVector[1]);
    SpatialDefRate = rGeom[2].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += 0.5 * (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + 2 * NormalVector[0] * SpatialDefRate[2] * NormalVector[1]);

    noalias(AccA) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(MeanAcc) = 0.5 * AccA + 0.5 * AccB;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);
  }
}

///TODO AS IN 2D MAKE THE 3D!!!
template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::ComputeElementalBoundRHSVector(VectorType &BoundRHSVector,
                                                                                                   const double TimeStep,
                                                                                                   const double BoundRHSCoeffAcc,
                                                                                                   const double BoundRHSCoeffDev)
{
  GeometryType &rGeom = this->GetGeometry();

  Vector SpatialDefRate = rGeom[0].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE) * 0.25 + rGeom[1].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE) * 0.25;
  SpatialDefRate += rGeom[2].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE) * 0.25 + rGeom[3].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE) * 0.25;

  array_1d<double, 3> AccAa(3, 0.0);
  array_1d<double, 3> AccBb(3, 0.0);
  array_1d<double, 3> AccCc(3, 0.0);
  array_1d<double, 3> AccDd(3, 0.0);
  array_1d<double, 3> MeanAcc(3, 0.0);

  noalias(AccAa) = 0.5 / TimeStep * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
  noalias(AccBb) = 0.5 / TimeStep * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
  noalias(AccCc) = 0.5 / TimeStep * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
  noalias(AccDd) = 0.5 / TimeStep * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

  noalias(MeanAcc) = 0.25 * (AccAa + AccBb + AccCc + AccDd);

  const double coeff = 0.25;
  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {

    // array_1d<double, 3> AccA(3,0.0);
    // array_1d<double, 3> AccB(3,0.0);
    // array_1d<double, 3> AccC(3,0.0);
    // array_1d<double, 3> MeanAcc(3,0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    // const double factor = 0.5/TimeStep;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 2, 3);

    // noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);
    // noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);
    // noalias(AccC)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);

    // noalias(MeanAcc)= one_third*AccA + one_third*AccB + one_third*AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    double elementalNormalProjDefRate = 0;
    // Vector SpatialDefRate =  rGeom[0].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    // SpatialDefRate =  rGeom[1].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    // SpatialDefRate =  rGeom[2].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    //     if(rGeom[0].IsNot(INLET))
    // BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*elementalNormalProjDefRate);

    //     if(rGeom[1].IsNot(INLET))
    // BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*elementalNormalProjDefRate);

    //     if(rGeom[2].IsNot(INLET))
    // BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*elementalNormalProjDefRate);

    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);
  }

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    // array_1d<double, 3> AccA(3,0.0);
    // array_1d<double, 3> AccB(3,0.0);
    // array_1d<double, 3> AccC(3,0.0);
    // array_1d<double, 3> MeanAcc(3,0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    // const double factor = 0.5/TimeStep;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 3, 2);

    // noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);
    // noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);
    // noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1);

    // noalias(MeanAcc)= one_third*AccA + one_third*AccB + one_third*AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    double elementalNormalProjDefRate = 0;
    // Vector SpatialDefRate =  rGeom[0].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    // SpatialDefRate =  rGeom[1].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    // SpatialDefRate =  rGeom[3].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);
  }

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    // array_1d<double, 3> AccA(3,0.0);
    // array_1d<double, 3> AccB(3,0.0);
    // array_1d<double, 3> AccC(3,0.0);
    // array_1d<double, 3> MeanAcc(3,0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    // const double factor = 0.5/TimeStep;
    // const double one_third = 1.0/3.0;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 2, 3, 1);

    // noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);
    // noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);
    // noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1);

    // noalias(MeanAcc)= one_third*AccA + one_third*AccB + one_third*AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    double elementalNormalProjDefRate = 0;
    // Vector SpatialDefRate =  rGeom[0].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    // SpatialDefRate =  rGeom[2].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    // SpatialDefRate =  rGeom[3].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);
  }

  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    // array_1d<double, 3> AccA(3,0.0);
    // array_1d<double, 3> AccB(3,0.0);
    // array_1d<double, 3> AccC(3,0.0);
    // array_1d<double, 3> MeanAcc(3,0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    // const double factor = 0.5/TimeStep;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 1, 2, 3, 0);

    // noalias(AccA)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);
    // noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);
    // noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1);

    // noalias(MeanAcc)= one_third*AccA + one_third*AccB + one_third*AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    double elementalNormalProjDefRate = 0;
    // Vector SpatialDefRate =  rGeom[1].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    // SpatialDefRate =  rGeom[2].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    // SpatialDefRate =  rGeom[3].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE);
    elementalNormalProjDefRate += (NormalVector[0] * SpatialDefRate[0] * NormalVector[0] + NormalVector[1] * SpatialDefRate[1] * NormalVector[1] + NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                   2.0 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] + 2.0 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] + 2.0 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2]) /
                                  3.0;

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);

    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * elementalNormalProjDefRate);
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
                                                                                                  const double TimeStep,
                                                                                                  const double BoundRHSCoeffAcc,
                                                                                                  const double BoundRHSCoeffDev,
                                                                                                  const VectorType SpatialDefRate)
{
  GeometryType &rGeom = this->GetGeometry();

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
  {
    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_half = 1.0 / 2.0;

    this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 1, 2);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);
    // double SpatialDefRateNormalProjection=0;
    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(MeanAcc) = 0.5 * AccA + 0.5 * AccB;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
      BoundRHSVector[0] += one_half * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += one_half * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_half = 1.0 / 2.0;

    this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 2, 1);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);
    // double SpatialDefRateNormalProjection=0;

    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(MeanAcc) = 0.5 * AccA + 0.5 * AccB;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
      BoundRHSVector[0] += one_half * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += one_half * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }

  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_half = 1.0 / 2.0;

    this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 1, 2, 0);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);
    // double SpatialDefRateNormalProjection=0;

    noalias(AccA) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(MeanAcc) = 0.5 * AccA + 0.5 * AccB;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += one_half * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += one_half * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }
}

template <>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
                                                                                                  const double TimeStep,
                                                                                                  const double BoundRHSCoeffAcc,
                                                                                                  const double BoundRHSCoeffDev,
                                                                                                  const VectorType SpatialDefRate)
{
  GeometryType &rGeom = this->GetGeometry();

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> AccC(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 2, 3);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);

    noalias(MeanAcc) = one_third * AccA + one_third * AccB + one_third * AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> AccC(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 3, 2);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

    noalias(MeanAcc) = one_third * AccA + one_third * AccB + one_third * AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }

  if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> AccC(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 2, 3, 1);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

    noalias(MeanAcc) = one_third * AccA + one_third * AccB + one_third * AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    if (rGeom[0].IsNot(INLET))
      BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }

  if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
  {

    array_1d<double, 3> AccA(3, 0.0);
    array_1d<double, 3> AccB(3, 0.0);
    array_1d<double, 3> AccC(3, 0.0);
    array_1d<double, 3> MeanAcc(3, 0.0);
    array_1d<double, 3> NormalVector(3, 0.0);
    const double factor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    this->GetOutwardsUnitNormalForThreePoints(NormalVector, 1, 2, 3, 0);

    double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    noalias(AccA) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

    noalias(MeanAcc) = one_third * AccA + one_third * AccB + one_third * AccC;

    const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    if (rGeom[1].IsNot(INLET))
      BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[2].IsNot(INLET))
      BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    if (rGeom[3].IsNot(INLET))
      BoundRHSVector[3] += one_third * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
  }
}

// template< unsigned int TDim >
// void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::ComputeStabLaplacianMatrix(MatrixType& StabLaplacianMatrix,
// 										const ShapeFunctionDerivativesType& rDN_DX,
// 										const double Weight)

// {
//   // LHS contribution
//   GeometryType& rGeom = this->GetGeometry();
//   const SizeType NumNodes =rGeom.PointsNumber();
//   unsigned int SFDposition=0;

//   double elementVolume=0;
//   if(TDim==3){
//    	elementVolume=rGeom.Volume()*0.25;
//   }else if(TDim==2){
//     elementVolume=rGeom.Area()/3.0;
//   }

//   for (SizeType i = 0; i < NumNodes; ++i)
//     {
//       double nodalVolumeI=rGeom[i].FastGetSolutionStepValue(NODAL_VOLUME);
//       VectorType nodalSFDneighboursI=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
//       unsigned int nodalSFDneighboursISize=nodalSFDneighboursI.size();
//       if(nodalSFDneighboursISize>1)
//     {
// 	    double nodal_dNdXi=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[0];
//       double nodal_dNdYi=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[1];
//       double nodal_dNdZi=0;
//       if(TDim==3){
//         nodal_dNdZi=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[2];
//       }
//       double elemental_dNdXi=0;
//       double elemental_dNdYi=0;
//       double elemental_dNdZi=0;

//       for (unsigned int ii = 0; ii< NumNodes; ii++)
//       {
//         // elemental_dNdXi=+rDN_DX(ii,0)*elementVolume/nodalVolumeI;
//         // elemental_dNdYi=+rDN_DX(ii,1)*elementVolume/nodalVolumeI;
//         unsigned int idNodeOfConsideredElement=rGeom[ii].Id();
//         SFDposition=0;
//         for (unsigned int k = 0; k< nodalSFDneighboursISize; k++)
// 	      {
//           if(idNodeOfConsideredElement==nodalSFDneighboursI[k])
//           {
//             elemental_dNdXi+=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition]/(double(NumNodes));
//             elemental_dNdYi+=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+1]/(double(NumNodes));
//               // std::cout<<"       elemental_dNdXi:"<<elemental_dNdXi<<"      elemental_dNdYi:"<<elemental_dNdYi<<std::endl;

//             if(TDim==3){
//               elemental_dNdZi+=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+2]/(double(NumNodes));
//             }
//             break;
//           }
//           SFDposition+=TDim;
//         }
//       }

//       for (SizeType j = 0; j < NumNodes; ++j)
//         {
//           double nodalVolumeJ=rGeom[j].FastGetSolutionStepValue(NODAL_VOLUME);
//           VectorType nodalSFDneighboursJ=rGeom[j].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
//           unsigned int nodalSFDneighboursJSize=nodalSFDneighboursJ.size();
//           if(nodalSFDneighboursJSize>1)
//           {

//           double nodal_dNdXj=rGeom[j].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[0];
//           double nodal_dNdYj=rGeom[j].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[1];
//           double nodal_dNdZj=0;
//           if(TDim==3){
//             nodal_dNdZj=rGeom[j].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[2];
//           }
//           double elemental_dNdXj=0;
//           double elemental_dNdYj=0;
//           double elemental_dNdZj=0;
//           for (unsigned int jj = 0; jj< NumNodes; jj++)
//            {
//             //  elemental_dNdXj=+rDN_DX(jj,0)*elementVolume/nodalVolumeJ;
//             //  elemental_dNdYj=+rDN_DX(jj,1)*elementVolume/nodalVolumeJ;
//              unsigned int idNodeOfConsideredElement=rGeom[jj].Id();
//              SFDposition=0;
//              for (unsigned int k = 0; k< nodalSFDneighboursJSize; k++)
//              {
//               if(idNodeOfConsideredElement==nodalSFDneighboursJ[k])
//               {
//                 elemental_dNdXj+=rGeom[j].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition]/(double(NumNodes));
//                 elemental_dNdYj+=rGeom[j].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+1]/(double(NumNodes));
//                 if(TDim==3){
//                   elemental_dNdZj+=rGeom[j].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+2]/(double(NumNodes));
//                 }
//                 break;
//               }
//               SFDposition+=TDim;
//             }
//           }

//           // std::cout<<"                 nodal_dNdXi:"<<nodal_dNdXi<<" nodal_dNdyi:"<<nodal_dNdYi<<std::endl;
//           // std::cout<<" elemental_dNdXi:"<<elemental_dNdXi<<" elemental_dNdyi:"<<elemental_dNdYi<<std::endl;
//           // std::cout<<" nodal_dNdXi:"<<nodal_dNdXj<<" nodal_dNdyJ:"<<nodal_dNdYj<<std::endl;
//           // std::cout<<" elemental_dNdXi:"<<elemental_dNdXj<<" elemental_dNdyJ:"<<elemental_dNdYj<<std::endl;

//           // StabLaplacianMatrix(i,j) += ( (nodal_dNdXi-rDN_DX(i,0))*(nodal_dNdXj-rDN_DX(j,0)) +
//           //                               (nodal_dNdYi-rDN_DX(i,1))*(nodal_dNdYj-rDN_DX(j,1)) ) * Weight;

//           StabLaplacianMatrix(i,j) += ( (nodal_dNdXi-elemental_dNdXi)*(nodal_dNdXj-elemental_dNdXj) +
//                                          (nodal_dNdYi-elemental_dNdYi)*(nodal_dNdYj-elemental_dNdYj) ) * Weight;

//           if(TDim==3){
//             StabLaplacianMatrix(i,j) += (nodal_dNdZi-elemental_dNdZi)*(nodal_dNdZj-elemental_dNdZj) * Weight ;
//           }

//           // double Lij = 0.0;
//           // for (SizeType d = 0; d < TDim; ++d){
//           //   Lij += rDN_DX(i,d) * rDN_DX(j,d);
//           // }
//           // StabLaplacianMatrix(i,j) += Weight * Lij ;
//         }
//         }
//     }
//     }
// }

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::ComputeStabLaplacianMatrix(MatrixType &StabLaplacianMatrix,
                                                                                                  const ShapeFunctionDerivativesType &rDN_DX,
                                                                                                  const double Weight)

{
  // LHS contribution
  const SizeType NumNodes = this->GetGeometry().PointsNumber();
  for (SizeType i = 0; i < NumNodes; ++i)
  {
    for (SizeType j = 0; j < NumNodes; ++j)
    {
      double Lij = 0.0;
      for (SizeType d = 0; d < TDim; ++d)
      {
        Lij += rDN_DX(i, d) * rDN_DX(j, d);
      }
      StabLaplacianMatrix(i, j) += Weight * Lij;
    }
  }
}

template <unsigned int TDim>
void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::AddStabilizationNodalTermsRHS(VectorType &rRightHandSideVector,
                                                                                                     const double Tau,
                                                                                                     const double Density,
                                                                                                     const double Weight,
                                                                                                     const ShapeFunctionDerivativesType &rDN_DX,
                                                                                                     const SizeType i)
{

  double RHSi = 0;
  if (this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION))
  { // it must be checked once at the begining only
    array_1d<double, 3> &VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

    // double posX=(this->GetGeometry()[0].X() + this->GetGeometry()[1].X() + this->GetGeometry()[2].X())/3.0;

    // double posY=(this->GetGeometry()[0].Y() + this->GetGeometry()[1].Y() + this->GetGeometry()[2].Y())/3.0;

    // double coeffX =(12.0-24.0*posY)*pow(posX,4);

    // coeffX += (-24.0+48.0*posY)*pow(posX,3);

    // coeffX += (-48.0*posY+72.0*pow(posY,2)-48.0*pow(posY,3)+12.0)*pow(posX,2);

    // coeffX += (-2.0+24.0*posY-72.0*pow(posY,2)+48.0*pow(posY,3))*posX;

    // coeffX += 1.0-4.0*posY+12.0*pow(posY,2)-8.0*pow(posY,3);

    // double coeffY =(8.0-48.0*posY+48.0*pow(posY,2))*pow(posX,3);

    // coeffY += (-12.0+72.0*posY-72.0*pow(posY,2))*pow(posX,2);

    // coeffY += (4.0-24.0*posY+48.0*pow(posY,2)-48.0*pow(posY,3)+24.0*pow(posY,4))*posX;

    // coeffY += -12.0*pow(posY,2)+24.0*pow(posY,3)-12.0*pow(posY,4);

    // RHSi += - Density  * Tau * ( rDN_DX(i,0) * VolumeAcceleration[0]*coeffX +  rDN_DX(i,1) * VolumeAcceleration[1] * coeffY );

    for (SizeType d = 0; d < TDim; ++d)
    {
      RHSi += -rDN_DX(i, d) * Tau * (Density * VolumeAcceleration[d]);
    }
  }
  rRightHandSideVector[i] += Weight * RHSi;
}

template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>;
template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>;

} // namespace Kratos
