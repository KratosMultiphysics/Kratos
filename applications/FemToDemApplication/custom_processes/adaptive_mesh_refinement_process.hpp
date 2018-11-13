// KRATOS FEM2DEM Application
//
//  License:         BSD License
//                   license: fem_to_dem_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Ignasi Pouplana
//

#if !defined(KRATOS_ADAPTIVE_MESH_REFINEMENT_PROCESS )
#define  KRATOS_ADAPTIVE_MESH_REFINEMENT_PROCESS

#include <fstream>
#include <cmath>

#include "includes/model_part.h"
//#include "boost/smart_ptr.hpp"
#include "processes/process.h"

#include "fem_to_dem_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/define.h"

namespace Kratos
{

//Only for Triangles2D3N or Quadrilateral2D4N

class AdaptiveMeshRefinementProcess : public Process
{
    
protected:

    struct NodeStresses
    {
        Vector EffectiveStressVector;
        int NElems;
      
        NodeStresses()
        {
            EffectiveStressVector = ZeroVector(3); //2D: sigma=[sx,sy,sxy]
            NElems = 0;
        }
    };

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

public:
  
    typedef ModelPart::ElementsContainerType ElementsArrayType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Constructor
    AdaptiveMeshRefinementProcess(ModelPart& r_model_part,std::string plane_state,std::string problem_name,
        std::string problem_path,std::string mesh_optimality_criteria,double permissible_error,
         int number_of_refinements) : mr_model_part(r_model_part)
    {
        mNNodes = mr_model_part.NumberOfNodes();
        mNElements = mr_model_part.NumberOfElements();
        mplane_state = plane_state;
        mproblem_name = problem_name;
        mproblem_path = problem_path;
        mmesh_optimality_criteria = mesh_optimality_criteria;
        mpermissible_error = permissible_error;
        mnumber_of_refinements = number_of_refinements;
    }
    
    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~AdaptiveMeshRefinementProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Generate a new mesh after a used mesh
    
    void Execute()
    {
        NodeStresses* NodeStressesVector = new NodeStresses[mNNodes];
        double* ElementError = new double[mNElements];
        double* ElementRefinementParameter = new double[mNElements];
        double* NewElementDimension = new double[mNElements];
        double GlobalError = 0.0, GlobalStrainEnergy = 0.0;
    
        this->StressExtrapolationAndSmoothing(NodeStressesVector);
        
        this->ErrorEstimationAndStrainEnergy(NodeStressesVector,ElementError,GlobalError,GlobalStrainEnergy);
        
        this->RefinementParameters(ElementError,GlobalError,GlobalStrainEnergy,ElementRefinementParameter,NewElementDimension);
        
        this->WriteNewMeshFiles(NewElementDimension);
        
        this->WriteLastMeshInfo(NodeStressesVector,ElementRefinementParameter,GlobalError,GlobalStrainEnergy);
        
        //Dellocate vectors:
        delete[] NodeStressesVector;
        delete[] ElementError;
        delete[] ElementRefinementParameter;
        delete[] NewElementDimension;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Generate a new mesh after an unused mesh
    
    void ExecuteAfterOutputStep()
    {
        NodeStresses* NodeStressesVector = new NodeStresses[mNNodes];
        double* ElementError = new double[mNElements];
        double* ElementRefinementParameter = new double[mNElements];
        double* NewElementDimension = new double[mNElements];
        double GlobalError = 0.0, GlobalStrainEnergy = 0.0;
        
        this->StressExtrapolationAndSmoothing(NodeStressesVector);
        
        this->ErrorEstimationAndStrainEnergy(NodeStressesVector,ElementError,GlobalError,GlobalStrainEnergy);
        
        this->RefinementParameters(ElementError,GlobalError,GlobalStrainEnergy,ElementRefinementParameter,NewElementDimension);
        
        this->WriteNewMeshFiles(NewElementDimension);
        
        //Dellocate vectors:
        delete[] NodeStressesVector;
        delete[] ElementError;
        delete[] ElementRefinementParameter;
        delete[] NewElementDimension;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Show information about the last mesh of the problem
    
    void ExecuteFinalize()
    {
        NodeStresses* NodeStressesVector = new NodeStresses[mNNodes];
        double* ElementError = new double[mNElements];
        double* ElementRefinementParameter = new double[mNElements];
        double* NewElementDimension = new double[mNElements];
        double GlobalError = 0.0, GlobalStrainEnergy = 0.0;
        
        this->StressExtrapolationAndSmoothing(NodeStressesVector);
        
        this->ErrorEstimationAndStrainEnergy(NodeStressesVector,ElementError,GlobalError,GlobalStrainEnergy);
        
        this->RefinementParameters(ElementError,GlobalError,GlobalStrainEnergy,ElementRefinementParameter,NewElementDimension);
        
        this->WriteLastMeshInfo(NodeStressesVector,ElementRefinementParameter,GlobalError,GlobalStrainEnergy);
        
        //Dellocate vectors:
        delete[] NodeStressesVector;
        delete[] ElementError;
        delete[] ElementRefinementParameter;
        delete[] NewElementDimension;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mr_model_part;

    std::string mplane_state, mproblem_name, mproblem_path, mmesh_optimality_criteria;
    double mpermissible_error;
    unsigned int mNNodes,mNElements;
    int mnumber_of_refinements;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    std::vector<double> ShapeFunctions(const double local_x,const double local_y)
    {
        std::vector<double> SF(4);

        SF[0] = (1-local_x-local_y+local_x*local_y)/4;
        SF[1] = (1+local_x-local_y-local_x*local_y)/4;
        SF[2] = (1+local_x+local_y+local_x*local_y)/4;
        SF[3] = (1-local_x+local_y-local_x*local_y)/4;

        return SF;
    }
  
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Matrix ShapeFunctionDerivativesMatrix(const double local_x,const double local_y)
    {
        Matrix SFD_Matrix = ZeroMatrix(2,4);

        SFD_Matrix(0,0) = (-1+local_y)/4;
        SFD_Matrix(0,1) = (1-local_y)/4;
        SFD_Matrix(0,2) = (1+local_y)/4;
        SFD_Matrix(0,3) = (-1-local_y)/4;

        SFD_Matrix(1,0) = (-1+local_x)/4;
        SFD_Matrix(1,1) = (-1-local_x)/4;
        SFD_Matrix(1,2) = (1+local_x)/4;
        SFD_Matrix(1,3) = (1-local_x)/4;

        return SFD_Matrix;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void StressExtrapolationAndSmoothing(NodeStresses* pNodeStressesVector)
    {
        //Basic declarations
        Vector GaussPointsStresses = ZeroVector(3);

        //Quadrilateral declarations
        Vector NodeCoord = ZeroVector(2);
        std::vector<Vector> NodeLocalCoordinates(4);
        std::vector<double> ExtrapolationComponents;
        Matrix Extrapolation = ZeroMatrix(4,4);
        Matrix QuadGaussStresses = ZeroMatrix(4,3);
        Matrix QuadNodeStresses = ZeroMatrix(4,3);

        NodeCoord[0] = -sqrt(3);
        NodeCoord[1] = -sqrt(3);
        NodeLocalCoordinates[0] = NodeCoord;
        NodeCoord[0] = sqrt(3);
        NodeCoord[1] = -sqrt(3);
        NodeLocalCoordinates[1] = NodeCoord;
        NodeCoord[0] = sqrt(3);
        NodeCoord[1] = sqrt(3);
        NodeLocalCoordinates[2] = NodeCoord;
        NodeCoord[0] = -sqrt(3);
        NodeCoord[1] = sqrt(3);
        NodeLocalCoordinates[3] = NodeCoord;

        for(int i = 0; i < 4; i++)
        {
            ExtrapolationComponents = this->ShapeFunctions( (NodeLocalCoordinates[i])[0] , (NodeLocalCoordinates[i])[1] );
            Extrapolation(i,0) = ExtrapolationComponents[0];
            Extrapolation(i,1) = ExtrapolationComponents[1];
            Extrapolation(i,2) = ExtrapolationComponents[2];
            Extrapolation(i,3) = ExtrapolationComponents[3];
        }
          
        
        for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
        {

			bool condition_is_active = true;
            if ((*it)->IsDefined(ACTIVE))
            {
                condition_is_active = (*it)->Is(ACTIVE);
            }

            if (condition_is_active)
            {
               GaussPointsStresses = (*it)->GetValue(STRESS_VECTOR);

                //Triangles2D3N
                if((*it)->GetGeometry().PointsNumber() == 3)
                {
                    for(int i = 0; i < 3; i++)
                    {   //KRATOS_WATCH(pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[0])
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[0] += GaussPointsStresses[0];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[1] += GaussPointsStresses[1];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[2] += GaussPointsStresses[2];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].NElems += 1;
                    }
                }
            }     
			//else // Inactive
			// {
			//	 if ((*it)->GetGeometry().PointsNumber() == 3)
			//	 {
			//		 for (int i = 0; i < 3; i++)
			//		 {   
			//			 pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].NElems += 1;
			//		 }
			//	 }
			// }
                // else //Quadrilateral2D4N
            // {
            //     for(int i=0;i<4;i++)
            //     {
            //         QuadGaussStresses(i,0) = (GaussPointsStresses[i])[0];
            //         QuadGaussStresses(i,1) = (GaussPointsStresses[i])[1];
            //         QuadGaussStresses(i,2) = (GaussPointsStresses[i])[2];
            //     }
       
            //     QuadNodeStresses = prod(Extrapolation,QuadGaussStresses);
       
            //     for(int i=0; i<4; i++)
            //     {
            //         pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[0] += QuadNodeStresses(i,0);
            //         pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[1] += QuadNodeStresses(i,1);
            //         pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[2] += QuadNodeStresses(i,2);
            //         pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].NElems += 1;
            //     }
            // }
        }

        for(unsigned int i = 0; i < mNNodes; i++)
        {
			if (pNodeStressesVector[i].NElems == 0) pNodeStressesVector[i].NElems = 1; // node surrounded by inactive elems

            pNodeStressesVector[i].EffectiveStressVector = pNodeStressesVector[i].EffectiveStressVector / pNodeStressesVector[i].NElems;
			//KRATOS_WATCH(i)
			//KRATOS_WATCH(pNodeStressesVector[i].EffectiveStressVector)
        }


        
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ErrorEstimationAndStrainEnergy(const NodeStresses* NodeStressesVector,double* pElementError,double& rGlobalError,double& rGlobalStrainEnergy)
    {
        //Basic delarations
        Vector GaussPointsStresses = ZeroVector(3);
        //std::vector<double>  GaussPointsStresses;
        Matrix InvElasticConstitutiveMatrix = ZeroMatrix(3,3);//Only 2D
        double young,poisson,aux1,aux2,aux3;
        double EError,EStrainEnergy;
        int Elem_it = 0;
        
        //Triangle declarations
        Vector SmoothedStress = ZeroVector(3);
        Vector StressDifference = ZeroVector(3);
        
        //Quadrilateral declarations
        Vector GPCoord = ZeroVector(2);
        std::vector<Vector> GPLocalCoordinates(4);
        double GPWeight = 1;
        std::vector<double> GaussPointWeights(2);
        std::vector<double> InterpolationComponents;
        Matrix Interpolation = ZeroMatrix(4,4);
        std::vector<Matrix> SFD_Matrices(4);
        Matrix QuadNodeStresses = ZeroMatrix(4,3);
        Matrix NodeCoordinates = ZeroMatrix(4,2);
        Matrix GPSmoothedStresses = ZeroMatrix(4,3);
        Matrix Jacobian = ZeroMatrix(2,2);
        double Det_Jacobian;
        std::vector<double> D_Area(4);
        std::vector<Vector> QuadStressDifference(4);
        std::vector<Vector> QuadSmoothedStress(4);

        GPCoord[0] = -sqrt(3)/3;
        GPCoord[1] = -sqrt(3)/3;
        GPLocalCoordinates[0] = GPCoord;
        GPCoord[0] = sqrt(3)/3;
        GPCoord[1] = -sqrt(3)/3;
        GPLocalCoordinates[1] = GPCoord;
        GPCoord[0] = sqrt(3)/3;
        GPCoord[1] = sqrt(3)/3;
        GPLocalCoordinates[2] = GPCoord;
        GPCoord[0] = -sqrt(3)/3;
        GPCoord[1] = sqrt(3)/3;
        GPLocalCoordinates[3] = GPCoord;

        GaussPointWeights[0] = GPWeight;
        GaussPointWeights[1] = GPWeight;

        for(int i = 0; i < 4; i++)
        {
            InterpolationComponents = this->ShapeFunctions( (GPLocalCoordinates[i])[0] , (GPLocalCoordinates[i])[1] );
            Interpolation(i,0) = InterpolationComponents[0];
            Interpolation(i,1) = InterpolationComponents[1];
            Interpolation(i,2) = InterpolationComponents[2];
            Interpolation(i,3) = InterpolationComponents[3];

            SFD_Matrices[i] = this->ShapeFunctionDerivativesMatrix( (GPLocalCoordinates[i])[0],(GPLocalCoordinates[i])[1] );
        }
          
        for (ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
        {
			bool condition_is_active = true;
			if ((*it)->IsDefined(ACTIVE))
			{
				condition_is_active = (*it)->Is(ACTIVE);
			}

			if (condition_is_active)
			{
				GaussPointsStresses = (*it)->GetValue(STRESS_VECTOR); // Cornejo

				young = (*it)->GetProperties()[YOUNG_MODULUS];
				poisson = (*it)->GetProperties()[POISSON_RATIO];

				if (mplane_state == "Plane_Stress")
				{
					aux1 = young / (1 - poisson*poisson);
					aux2 = poisson*aux1;
					aux3 = young / (2 * (1 + poisson));
				}
				else
				{
					aux1 = young*(1 - poisson) / ((1 + poisson)*(1 - 2 * poisson));
					aux2 = aux1*poisson / (1 - poisson);
					aux3 = young / (2 * (1 + poisson));
				}

				InvElasticConstitutiveMatrix(0, 0) = aux1 / (aux1*aux1 - aux2*aux2);
				InvElasticConstitutiveMatrix(0, 1) = -aux2 / (aux1*aux1 - aux2*aux2);
				InvElasticConstitutiveMatrix(1, 0) = -aux2 / (aux1*aux1 - aux2*aux2);
				InvElasticConstitutiveMatrix(1, 1) = aux1 / (aux1*aux1 - aux2*aux2);
				InvElasticConstitutiveMatrix(2, 2) = 1 / aux3;

				//Triangles2D3N
				if ((*it)->GetGeometry().PointsNumber() == 3)
				{
					SmoothedStress = ZeroVector(3);

					for (int i = 0; i < 3; i++)
					{
						SmoothedStress += NodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].EffectiveStressVector;
					}

					SmoothedStress /= 3;

					StressDifference = SmoothedStress - GaussPointsStresses;
					//KRATOS_WATCH(StressDifference)
					EError = sqrt(inner_prod((prod(StressDifference, InvElasticConstitutiveMatrix)), StressDifference)*(*it)->GetGeometry().Area());
					//KRATOS_WATCH(EError)

						EStrainEnergy = sqrt(inner_prod((prod(SmoothedStress, InvElasticConstitutiveMatrix)), SmoothedStress)*(*it)->GetGeometry().Area());
				}
				// else //Quadrilaterals2D4N
				// {
				//     for(int i=0; i<4; i++)
				//     {
				//         QuadNodeStresses(i,0) = NodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[0];
				//         QuadNodeStresses(i,1) = NodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[1];
				//         QuadNodeStresses(i,2) = NodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id()-1].EffectiveStressVector[2];

				//         NodeCoordinates(i,0) = (*it)->GetGeometry().GetPoint(i).X();
				//         NodeCoordinates(i,1) = (*it)->GetGeometry().GetPoint(i).Y();

				//         QuadStressDifference[i] = ZeroVector(3);
				//         QuadSmoothedStress[i] = ZeroVector(3);
				//     }

				//     GPSmoothedStresses = prod(Interpolation,QuadNodeStresses);

				//     for(int i=0; i<4; i++)
				//     {
				//         Jacobian = prod(SFD_Matrices[i],NodeCoordinates);
				//         Det_Jacobian = Jacobian(0,0)*Jacobian(1,1) - Jacobian(0,1)*Jacobian(1,0);
				//         D_Area[i] = Det_Jacobian*GaussPointWeights[0]*GaussPointWeights[1];

				//         (QuadStressDifference[i])[0] = GPSmoothedStresses(i,0) - (GaussPointsStresses[i])[0];
				//         (QuadStressDifference[i])[1] = GPSmoothedStresses(i,1) - (GaussPointsStresses[i])[1];
				//         (QuadStressDifference[i])[2] = GPSmoothedStresses(i,2) - (GaussPointsStresses[i])[2];          

				//         (QuadSmoothedStress[i])[0] = GPSmoothedStresses(i,0);
				//         (QuadSmoothedStress[i])[1] = GPSmoothedStresses(i,1);
				//         (QuadSmoothedStress[i])[2] = GPSmoothedStresses(i,2);
				//     }

				//     EError = sqrt(inner_prod((prod(QuadStressDifference[0],InvElasticConstitutiveMatrix)),QuadStressDifference[0])*D_Area[0] +
				//           inner_prod((prod(QuadStressDifference[1],InvElasticConstitutiveMatrix)),QuadStressDifference[1])*D_Area[1] +
				//           inner_prod((prod(QuadStressDifference[2],InvElasticConstitutiveMatrix)),QuadStressDifference[2])*D_Area[2] +
				//           inner_prod((prod(QuadStressDifference[3],InvElasticConstitutiveMatrix)),QuadStressDifference[3])*D_Area[3]);

				//     EStrainEnergy = sqrt(inner_prod((prod(QuadSmoothedStress[0],InvElasticConstitutiveMatrix)),QuadSmoothedStress[0])*D_Area[0] +
				//                  inner_prod((prod(QuadSmoothedStress[1],InvElasticConstitutiveMatrix)),QuadSmoothedStress[1])*D_Area[1] +
				//                  inner_prod((prod(QuadSmoothedStress[2],InvElasticConstitutiveMatrix)),QuadSmoothedStress[2])*D_Area[2] +
				//                  inner_prod((prod(QuadSmoothedStress[3],InvElasticConstitutiveMatrix)),QuadSmoothedStress[3])*D_Area[3]);
				// }

				pElementError[Elem_it] = EError;
				rGlobalError += EError*EError;
				rGlobalStrainEnergy += EStrainEnergy*EStrainEnergy;

				Elem_it += 1;
			}
			else Elem_it += 1;
        }

        rGlobalError = sqrt(rGlobalError);
        rGlobalStrainEnergy = sqrt(rGlobalStrainEnergy);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void RefinementParameters(const double* ElementError,const double& GlobalError,const double& GlobalStrainEnergy,double* pElementRefinementParameter,double* pNewElementDimension)
    {
        //Basic declarations
        double AuxEDim, EDim;
        double* ElementDimension = new double[mNElements];
        double m_coef = 1;
        double d_coef = 2;
        double q_coef = (2*m_coef + d_coef) / 2;
        int Elem_it = 0;
        double TotalArea = 0.0;

        for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
        {
            EDim = sqrt(((*it)->GetGeometry().GetPoint((*it)->GetGeometry().PointsNumber()-1).X()-(*it)->GetGeometry().GetPoint(0).X())*
                  ((*it)->GetGeometry().GetPoint((*it)->GetGeometry().PointsNumber()-1).X()-(*it)->GetGeometry().GetPoint(0).X()) +
                  ((*it)->GetGeometry().GetPoint((*it)->GetGeometry().PointsNumber()-1).Y()-(*it)->GetGeometry().GetPoint(0).Y())*
                  ((*it)->GetGeometry().GetPoint((*it)->GetGeometry().PointsNumber()-1).Y()-(*it)->GetGeometry().GetPoint(0).Y()));
      
            for(unsigned int i = 0; i < ((*it)->GetGeometry().PointsNumber()-1); i++)
            {
                AuxEDim = sqrt(((*it)->GetGeometry().GetPoint(i).X()-(*it)->GetGeometry().GetPoint(i+1).X())*
                       ((*it)->GetGeometry().GetPoint(i).X()-(*it)->GetGeometry().GetPoint(i+1).X()) +
                       ((*it)->GetGeometry().GetPoint(i).Y()-(*it)->GetGeometry().GetPoint(i+1).Y())*
                       ((*it)->GetGeometry().GetPoint(i).Y()-(*it)->GetGeometry().GetPoint(i+1).Y()));
        
                if(AuxEDim < EDim) EDim = AuxEDim;
            }
      
            ElementDimension[Elem_it] = EDim;
      
            TotalArea += (*it)->GetGeometry().Area();
      
            Elem_it += 1;
        }
    
        Elem_it = 0;
        
    
        if(mmesh_optimality_criteria == "Global_Error_Equidistribution")
        {
            for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
            {
                bool condition_is_active = true;
                if ((*it)->IsDefined(ACTIVE))
                {
                    condition_is_active = (*it)->Is(ACTIVE);
                }


                //std::cout << (*it)->Is(ACTIVE) << "   " << (*it)->Id() << std::endl;

                double damage = (*it)->GetValue(DAMAGE_ELEMENT);
                
                if (condition_is_active)
                {
                    pElementRefinementParameter[Elem_it] = pow(GlobalError/(mpermissible_error*GlobalStrainEnergy),1/m_coef)*pow(sqrt(mNElements)*ElementError[Elem_it]/GlobalError,1/q_coef);

                    if (damage > 0.0  && pElementRefinementParameter[Elem_it] < 1.0) // damaged elements remain with constant size
                    {
                        pElementRefinementParameter[Elem_it] = 1.2; // to maintain the crack width
                    }

                    if (damage > 0.965 && pElementRefinementParameter[Elem_it] < 4.0) // we refine the crack lips elements
                    {
                        pElementRefinementParameter[Elem_it] = 4.0; // to maintain the crack width
                    } 

                    // if (pElementRefinementParameter[Elem_it] < 1.0) // test
                    // {
                    //     pElementRefinementParameter[Elem_it] = 1.0;
                    // }

                    pNewElementDimension[Elem_it] =  ElementDimension[Elem_it] / pElementRefinementParameter[Elem_it];


                }
                else // Inactive
                {
                    pElementRefinementParameter[Elem_it] = 10.0;
                    //pNewElementDimension[Elem_it] = ElementDimension[Elem_it]/4;
                    //pNewElementDimension[Elem_it] = ElementDimension[Elem_it] / pElementRefinementParameter[Elem_it];
					//pNewElementDimension[Elem_it] = ElementDimension[Elem_it] / pElementRefinementParameter[Elem_it];
					pNewElementDimension[Elem_it] =  ElementDimension[Elem_it] / pElementRefinementParameter[Elem_it];
                }
                    
                Elem_it += 1;
            }
        }
        else
        {
            for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
            {
                pElementRefinementParameter[Elem_it] = pow(ElementError[Elem_it]/(mpermissible_error*GlobalStrainEnergy)*
                                                  sqrt(TotalArea/((*it)->GetGeometry().Area())),1/m_coef);
                
                               
                pNewElementDimension[Elem_it] = ElementDimension[Elem_it] / pElementRefinementParameter[Elem_it];
                Elem_it += 1;
            }
        }
    
        //Dellocate vector:
        delete[] ElementDimension;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void WriteNewMeshFiles(const double* NewElementDimension)
    {
        std::string eletyp;

        if((*(mr_model_part.Elements().ptr_begin()))->GetGeometry().PointsNumber()==3) //Only one type of element (triangle or quadrilateral)
            eletyp = "Triangle";
        else
            eletyp = "Tetrahedra"; // TODO: this means Quadrilateral (bug)
  
        //Writing background mesh file .BGM
        std::fstream bgmfile;
        bgmfile.open((mproblem_name+".bgm").c_str(),std::fstream::out);
  
        bgmfile << "BackgroundMesh V 1.0" << std::endl;
        bgmfile << "MESH dimension " << 2 << " ElemType " << eletyp << " Nnode " << (*(mr_model_part.Elements().ptr_begin()))->GetGeometry().PointsNumber() << std::endl;
        bgmfile << "Coordinates" << std::endl;
        for(ModelPart::NodeIterator i = mr_model_part.NodesBegin(); i != mr_model_part.NodesEnd(); ++i)
            bgmfile << i->Id() << " " << i->X0() << " " << i->Y0() << std::endl;
            
        bgmfile << "end coordinates" << std::endl << std::endl;
  
        bgmfile << "Elements" << std::endl;
  
        if(eletyp == "Triangle")
        {
            for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
            {
                bgmfile << (*it)->Id() << " " << (*it)->GetGeometry().GetPoint(0).Id() << " " << (*it)->GetGeometry().GetPoint(1).Id() << " " 
                        << (*it)->GetGeometry().GetPoint(2).Id() << std::endl;
            }
        }
        else
        {
            for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
            {
                bgmfile << (*it)->Id() << " " << (*it)->GetGeometry().GetPoint(0).Id() << " " << (*it)->GetGeometry().GetPoint(1).Id() << " " 
                        << (*it)->GetGeometry().GetPoint(2).Id() << " " << (*it)->GetGeometry().GetPoint(3).Id() << std::endl;
            }
        }
  
        bgmfile << "end elements" << std::endl << std::endl;
  
        bgmfile << "DesiredSize Elements" << std::endl;
        
        double GlobalDimension = 0.0;
        
        for(unsigned int i = 0; i < mNElements; i++)
        {
            bgmfile << i + 1 << " " << NewElementDimension[i] << std::endl;
            GlobalDimension += NewElementDimension[i];
        }
        GlobalDimension = GlobalDimension / mNElements;
        bgmfile << "End DesiredSize" << std::endl;
  
        bgmfile.close();
  
        // Writing batch file for GiD
        std::fstream batchfile;
      
        batchfile.open((mproblem_name+".bch").c_str(),std::fstream::out);
      
		batchfile << "Mescape Files Read" << std::endl;
        batchfile << mproblem_path << std::endl;
        batchfile << "Mescape Meshing AssignSizes BackgMesh" << std::endl;
        batchfile << "Yes" << std::endl;
        batchfile << "escape" << std::endl;
		batchfile << "Mescape Meshing AssignSizes BackgMesh" << std::endl;
		batchfile << mproblem_path + "\\" + mproblem_name + ".bgm" << std::endl;
        batchfile << "Mescape Meshing Generate" << std::endl;
        batchfile << "Yes" << std::endl;
        batchfile << GlobalDimension << " MeshingParametersFrom=Preferences" << std::endl;
        batchfile << "Mescape Meshing MeshView" << std::endl;
        batchfile << "Mescape Utilities Calculate" << std::endl; //added
		batchfile << "Mescape Files WriteCalcFile" << std::endl;
		//batchfile << mproblem_path + "\\" + mproblem_name + ".dat" << std::endl;
        batchfile << "escape escape escape escape escape Quit" << std::endl;
        batchfile << "No" << std::endl;
      
        batchfile.close();
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void WriteLastMeshInfo(const NodeStresses* NodeStressesVector,const double* ElementRefinementParameter,const double& GlobalError,const double& GlobalStrainEnergy)
    {
        // AMR PostProcess Files 
        std::string eletyp;

        if((*(mr_model_part.Elements().ptr_begin()))->GetGeometry().PointsNumber()==3)//Only one type of element (triangle or quadrilateral)
            eletyp = "Triangle";
        else
            eletyp = "Quadrilateral";
        
        // Writing Post Mesh File
        std::fstream meshfile;
        meshfile.open((mproblem_name + "_AMR_parameters.post.msh").c_str(),std::fstream::out);
    
        meshfile << "MESH \"AMR_Mesh\" dimension " << 2 << " ElemType " << eletyp << " Nnode " << (*(mr_model_part.Elements().ptr_begin()))->GetGeometry().PointsNumber() << std::endl;
        meshfile << "Coordinates" << std::endl;
        
        for(ModelPart::NodeIterator i = mr_model_part.NodesBegin(); i != mr_model_part.NodesEnd(); ++i)
            meshfile << i->Id() << " " << i->X0() << " " << i->Y0() << std::endl;
        
        meshfile << "End Coordinates" << std::endl << std::endl;
    
        meshfile << "Elements" << std::endl;
    
        if(eletyp == "Triangle")
        {
            for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
            {
                meshfile << (*it)->Id() << " " << (*it)->GetGeometry().GetPoint(0).Id() << " " << (*it)->GetGeometry().GetPoint(1).Id() << " " 
                         << (*it)->GetGeometry().GetPoint(2).Id() << std::endl;
            }
        }
        else
        {
            for(ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
            {
                meshfile << (*it)->Id() << " " << (*it)->GetGeometry().GetPoint(0).Id() << " " << (*it)->GetGeometry().GetPoint(1).Id() << " " 
                         << (*it)->GetGeometry().GetPoint(2).Id() << " " << (*it)->GetGeometry().GetPoint(3).Id() << std::endl;
            }
        }
    
        meshfile << "End Elements" << std::endl << std::endl;
    
        meshfile.close();
    
        // Writing Results File
        std::fstream resfile;
        resfile.open((mproblem_name+"_AMR_parameters.post.res").c_str(),std::fstream::out);
    
        resfile << "GiD Post Results File 1.0" << std::endl;
    
        resfile << "GaussPoints \"AMR_PostProcess\" Elemtype " << eletyp << std::endl;
        resfile << "Number of Gauss Points: " << 1 << std::endl;
        resfile << "Natural Coordinates: Internal" << std::endl;
        resfile << "End GaussPoints" << std::endl << std::endl;
        
        resfile << "Result \"Nodal_Effective_Stresses\" \"Kratos_AMR\" " << 1 << " Vector OnNodes" << std::endl;
        resfile << "ComponentNames \"Sxx\", \"Syy\", \"Sxy\"" << std::endl;
        resfile << "Values" << std::endl;
        for(unsigned int i = 0; i < mNNodes; i++)
        {
          resfile << i+1 << " " << NodeStressesVector[i].EffectiveStressVector[0] << " " << NodeStressesVector[i].EffectiveStressVector[1] << " " << NodeStressesVector[i].EffectiveStressVector[2] << std::endl;
        }
        resfile << "End Values" << std::endl << std::endl;
        
        resfile << "Result \"Refinement_Parameter\" \"Kratos_AMR\" " << 1 << " Scalar OnGaussPoints AMR_PostProcess" << std::endl;
        resfile << "Values" << std::endl;
        for(unsigned int i = 0; i < mNElements; i++)
        {
          resfile << i+1 << " " << ElementRefinementParameter[i] << std::endl;
        }
        resfile << "End Values" << std::endl;

        resfile.close();
        
        //------------------------------------------------------------------------------------
        
        //AMR Info
        
        double GlobalRefinementParameter = GlobalError/(mpermissible_error*GlobalStrainEnergy);
        
        //Writing AMR info file
        std::fstream AMR_info;
        AMR_info.open("AMR_info.txt", std::fstream::out | std::fstream::app);
    
        if(mnumber_of_refinements == 0)
            AMR_info << "################ Adaptive Mesh Refinement Info ################" << std::endl;
            
        AMR_info << "Mesh: " << mnumber_of_refinements << std::endl;
        AMR_info << "Number of elements: " << mNElements << std::endl;
        AMR_info << "Global Refinement Parameter: " << GlobalRefinementParameter << std::endl;
        AMR_info << "------------------------------------------------------------------" << std::endl;
    
        AMR_info.close();
        
        //------------------------------------------------------------------------------------
        
        std::cout << "################ Adaptive Mesh Refinement Info ################" << std::endl;
        std::cout << "Number of elements: " << mNElements << std::endl;
        std::cout << "Global Refinement Parameter: " << GlobalRefinementParameter << std::endl;
        std::cout << "###############################################################" << std::endl;
    }

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_ADAPTIVE_MESH_REFINEMENT_PROCESS defined */
