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
*   Date:                $Date: 2009-01-14 17:15:24 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/

#if !defined(KRATOS_RESTART_UTILITY_INCLUDED )
#define  KRATOS_RESTART_UTILITY_INCLUDED
//System includes
#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>

#if !defined(isnan)
#define isnan(x) ((x)!=(x))
#endif

//External includes
#include "boost/smart_ptr.hpp"

//Project includes
#include "includes/model_part.h"
#include "includes/variables.h"
#include "integration/integration_point.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{
	class RestartUtility
	{
		public:
			typedef ModelPart::ElementsContainerType ElementsArrayType;
			typedef Geometry<Node<3> >::IntegrationPointsArrayType IntegrationPointsArrayType;
			/** 
			 * Constructor.
			 */
			RestartUtility(const std::string& rDatafilename)
            {
				mDataFileName= rDatafilename;
				std::cout << "RestartUtility created" << std::endl;
			}
			
			/** 
			 * Destructor.
			 */
			virtual ~RestartUtility()
			{}
			
			void ChangeFileName( const std::string& rDatafilename)
			{
				mDataFileName= rDatafilename;
			}
			/**
			 */
			void StoreNodalVariables( ModelPart& source_model_part )
			{
				std::string nodesfile= mDataFileName+".nodes.restart";
				std::ofstream outfile;
				outfile.open((char*)nodesfile.c_str());
				outfile.precision(14);
				outfile<<"TIME["<<source_model_part.GetProcessInfo()[TIME]<<"]\n";
				for( ModelPart::NodeIterator it = source_model_part.NodesBegin() ;
								 it != source_model_part.NodesEnd(); it++ )
				{
					outfile<<"NODE["<<(*it).Id()<<"]\n";
					if(it->HasDofFor(DISPLACEMENT_X)
						|| it->HasDofFor(DISPLACEMENT_Y)
						|| it->HasDofFor(DISPLACEMENT_Z))
					{
						outfile<<"\tDISPLACEMENT_NULL["<<it->GetSolutionStepValue(DISPLACEMENT_NULL)[0]<<","<<it->GetSolutionStepValue(DISPLACEMENT_NULL)[1]<<","<<it->GetSolutionStepValue(DISPLACEMENT_NULL)[2]<<"]\n";
						outfile<<"\tDISPLACEMENT_NULL_DT["<<it->GetSolutionStepValue(DISPLACEMENT_NULL_DT)[0]<<","<<it->GetSolutionStepValue(DISPLACEMENT_NULL_DT)[1]<<","<<it->GetSolutionStepValue(DISPLACEMENT_NULL_DT)[2]<<"]\n";
						outfile<<"\tACCELERATION_NULL["<<it->GetSolutionStepValue(ACCELERATION_NULL)[0]<<","<<it->GetSolutionStepValue(ACCELERATION_NULL)[1]<<","<<it->GetSolutionStepValue(ACCELERATION_NULL)[2]<<"]\n";
						outfile<<"\tDISPLACEMENT_OLD["<<it->GetSolutionStepValue(DISPLACEMENT_OLD)[0]<<","<<it->GetSolutionStepValue(DISPLACEMENT_OLD)[1]<<","<<it->GetSolutionStepValue(DISPLACEMENT_OLD)[2]<<"]\n";
					}
					if(it->HasDofFor(WATER_PRESSURE))
					{
						outfile<<"\tWATER_PRESSURE_NULL["<<it->GetSolutionStepValue(WATER_PRESSURE_NULL)<<"]]\n";
						outfile<<"\tWATER_PRESSURE_NULL_DT["<<it->GetSolutionStepValue(WATER_PRESSURE_NULL_DT)<<"]\n";
						outfile<<"\tWATER_PRESSURE_NULL_ACCELERATION["<<it->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION)<<"]\n";
					}
					if(it->HasDofFor(AIR_PRESSURE))
					{
						outfile<<"\tAIR_PRESSURE_NULL["<<it->GetSolutionStepValue(AIR_PRESSURE_NULL)<<"]\n";
						outfile<<"\tAIR_PRESSURE_NULL_DT["<<it->GetSolutionStepValue(AIR_PRESSURE_NULL_DT)<<"]\n";
						outfile<<"\tAIR_PRESSURE_NULL_ACCELERATION["<<it->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION)<<"]\n";
					}
				}
				outfile.close();
			}

			/**
			 */
			double WriteNodalVariables( ModelPart& model_part )
			{
				double time = 0.0;
				unsigned int id;
				char c;
				std::string line;
				std::string out;
				std::string nodesfile= mDataFileName+".nodes.restart";
				std::ifstream infile;
				infile.precision(14);
				infile.open((char*)nodesfile.c_str());
				infile.seekg(0,std::ios::beg);

				while (!infile.eof()) 
				{
					out="";
					infile.get(c);
					if(c!= '\t')
					{
						out+=c;
						while(c!='\n')
						{
							infile.get(c);
							if(c=='[') break;
							out+=c;
						}
						if(out=="TIME")
						{
							out="";
							while(c!='\n')
							{
								infile.get(c);
								if(c==']') break;
								out+=c;
							}

							time= std::atof((char*)out.c_str());
							break;
						}
  					}
					std::getline(infile, line);
				}
				for( ModelPart::NodeIterator it = model_part.NodesBegin() ;
								 it != model_part.NodesEnd(); it++ )
				{
					bool node_found= false;
					infile.seekg(0,std::ios::beg);

					while (!infile.eof() ) 
					{
							out="";
							infile.get(c);
							if(c!= '\t')
							{
								out+=c;
								while(c!='\n')
								{
									infile.get(c);
									if(c=='[') break;
									out+=c;
								}
							}
							else
							{
								std::getline(infile, line);
								continue;
							}

							if(out=="NODE")
							{
								out="";
								while(c!='\n')
								{
									infile.get(c);
									if(c==']') break;
									out+=c;
								}
								id= std::atoi((char*)out.c_str());
								std::getline(infile, line);
							}
							else
							{
								std::getline(infile, line);
								continue;
							}
							if(id== it->Id())
							{
								node_found= true;
								if(it->HasDofFor(DISPLACEMENT_X)
									|| it->HasDofFor(DISPLACEMENT_Y)
									|| it->HasDofFor(DISPLACEMENT_Z))
								{
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									for(int k=0; k<3; k++)
									{
										out="";
										while(c!='\n')
										{
											infile.get(c);
											if(c==']'|| c== ',') break;
											out+=c;
										}
										double value= std::atof((char*)out.c_str());
										if(isnan(value))
										{
											std::cout<<"##### NaN FOUND ("<<it->Id()<<") ####"<<std::endl;
											continue;
										}
										it->GetSolutionStepValue(DISPLACEMENT_NULL)[k] = value;
									}
									std::getline(infile, line);
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									for(int k=0; k<3; k++)
									{
										out="";
										while(c!='\n')
										{
											infile.get(c);
											if(c==']'|| c== ',') break;
											out+=c;
										}
										double value= std::atof((char*)out.c_str());
										if(isnan(value))
										{
											std::cout<<"##### NaN FOUND ("<<it->Id()<<") ####"<<std::endl;
											continue;
										}
										it->GetSolutionStepValue(DISPLACEMENT_NULL_DT)[k] =value;
									}
									std::getline(infile, line);
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									for(int k=0; k<3; k++)
									{
										out="";
										while(c!='\n')
										{
											infile.get(c);
											if(c==']'|| c== ',') break;
											out+=c;
										}
										double value= std::atof((char*)out.c_str());
										if(isnan(value))
										{
											std::cout<<"##### NaN FOUND ("<<it->Id()<<") ####"<<std::endl;
											continue;
										}
										it->GetSolutionStepValue(ACCELERATION_NULL)[k] = value;
									}
									std::getline(infile, line);
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									for(int k=0; k<3; k++)
									{
										out="";
										while(c!='\n')
										{
											infile.get(c);
											if(c==']'|| c== ',') break;
											out+=c;
										}
										double value= std::atof((char*)out.c_str());
										if(isnan(value))
										{
											std::cout<<"##### NaN FOUND ("<<it->Id()<<") ####"<<std::endl;
											continue;
										}
										it->GetSolutionStepValue(DISPLACEMENT_OLD)[k] = value;
									}
									std::getline(infile, line);
								}
								if(it->HasDofFor(WATER_PRESSURE))
								{
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									out="";
									while(c!='\n')
									{
										infile.get(c);
										if(c==']'|| c== ',') break;
										out+=c;
									}
									double value= std::atof((char*)out.c_str());
									if(isnan(value))
									{
										std::cout<<"##### NaN FOUND ("<<it->Id()<<") ####"<<std::endl;
									}
									else
										it->GetSolutionStepValue(WATER_PRESSURE_NULL) = value;
									std::getline(infile, line);
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									out="";
									while(c!='\n')
									{
										infile.get(c);
										if(c==']'|| c== ',') break;
										out+=c;
									}
									value= std::atof((char*)out.c_str());
									if(isnan(value))
									{
										std::cout<<"##### NaN FOUND ("<<it->Id()<<") ####"<<std::endl;
									}
									else
										it->GetSolutionStepValue(WATER_PRESSURE_NULL_DT) = value;
									std::getline(infile, line);
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									out="";
									while(c!='\n')
									{
										infile.get(c);
										if(c==']'|| c== ',') break;
										out+=c;
									}
									value= std::atof((char*)out.c_str());
									if(isnan(value))
									{
										std::cout<<"##### NaN FOUND ("<<it->Id()<<") ####"<<std::endl;
									}
									else
										it->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION) = value;
									std::getline(infile, line);
								}
								if(it->HasDofFor(AIR_PRESSURE))
								{
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									out="";
									while(c!='\n')
									{
										infile.get(c);
										if(c==']'|| c== ',') break;
										out+=c;
									}
									double value= std::atof((char*)out.c_str());
									if(isnan(value))
									{
										std::cout<<"##### NaN FOUND ("<<it->Id()<<") ####"<<std::endl;
									}
									else
										it->GetSolutionStepValue(AIR_PRESSURE_NULL) = value;
									std::getline(infile, line);
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									out="";
									while(c!='\n')
									{
										infile.get(c);
										if(c==']'|| c== ',') break;
										out+=c;
									}
									value= std::atof((char*)out.c_str());
									if(isnan(value))
									{
										std::cout<<"##### NaN FOUND ("<<it->Id()<<") ####"<<std::endl;
									}
									else
										it->GetSolutionStepValue(AIR_PRESSURE_NULL_DT) = value;
									std::getline(infile, line);
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									out="";
									while(c!='\n')
									{
										infile.get(c);
										if(c==']'|| c== ',') break;
										out+=c;
									}
									value= std::atof((char*)out.c_str());
									if(isnan(value))
									{
										std::cout<<"##### NaN FOUND ("<<it->Id()<<") ####"<<std::endl;
									}
									else
										it->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION) = value;
									std::getline(infile, line);
								}
								break;
							}
					}

					if(!node_found)
						std::cout<<"##### NO CORRESPONDING NODE ("<<it->Id()<<") IN INPUT FILE ####"<<std::endl;
				}
				infile.close();
				return time;
			}
			/**
			 */
			void StoreInSituStress( ModelPart& model_part)
			{
				ElementsArrayType& ElementsArray= model_part.Elements();

				std::string nodesfile= mDataFileName+".insitu_stress.restart";
				std::ofstream outfile;
				outfile.open((char*)nodesfile.c_str());
				outfile.precision(14);
				outfile<<"TIME["<<model_part.GetProcessInfo()[TIME]<<"]\n";
			
				for( ElementsArrayType::ptr_iterator it = ElementsArray.ptr_begin(); 
				 	it != ElementsArray.ptr_end(); ++it )
				{
					const IntegrationPointsArrayType& integration_points 
						= (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());
					outfile<<"ELEMENT["<<(*it)->Id()<<"]\n";

					std::vector<Vector> ValuesOnIntPoint(integration_points.size());
					
					(*it)->GetValueOnIntegrationPoints(INSITU_STRESS, ValuesOnIntPoint, model_part.GetProcessInfo());

					if(ValuesOnIntPoint[0].size()!=6)
					{
						outfile<<"\t[NOVALUE]\n";
						continue;
					}

					for(unsigned int point=0; point< integration_points.size(); point++)
					{
						outfile<<"\t["<<ValuesOnIntPoint[point](0)<<","<<ValuesOnIntPoint[point](1)<<","<<ValuesOnIntPoint[point](2)<<","<<ValuesOnIntPoint[point](3)<<","<<ValuesOnIntPoint[point](4)<<","<<ValuesOnIntPoint[point](5)<<"]\n";
					}
				}

				outfile.close();
			}

			/**
			 */
			void WriteInSituStress( ModelPart& model_part)
			{
				unsigned int id;
				char c;
				std::string line;
				std::string out;
				std::string nodesfile= mDataFileName+".insitu_stress.restart";
				std::ifstream infile;
				infile.open((char*)nodesfile.c_str());
				infile.precision(14);
				infile.seekg(0,std::ios::beg);

				ElementsArrayType& ElementsArray= model_part.Elements();

				for( ElementsArrayType::ptr_iterator it = ElementsArray.ptr_begin(); 
				 	it != ElementsArray.ptr_end(); ++it )
				{
					bool element_found= false;

					const IntegrationPointsArrayType& integration_points 
						= (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

					std::vector<Vector> ValuesOnIntPoint(integration_points.size());

					infile.seekg(0,std::ios::beg);

					while (!infile.eof() ) 
					{
							out="";
							infile.get(c);
							if(c!= '\t')
							{
								out +=c;
								while(c!='\n')
								{
									infile.get(c);
									if(c=='[') break;
									out +=c;
								}
							}
							else
							{
								std::getline(infile, line);
								continue;
							}

							if(out=="ELEMENT")
							{
								out="";
								while(c!='\n')
								{
									infile.get(c);
									if(c==']') break;
									out+=c;
								}
								id= std::atoi((char*)out.c_str());

								std::getline(infile, line);
							}
							else
							{
								std::getline(infile, line);
								continue;
							}
							if(id== (*it)->Id())
							{
								element_found=true;
								for(unsigned int point=0; point< integration_points.size(); point++)
								{
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									if(ValuesOnIntPoint[point].size()!=6)
										ValuesOnIntPoint[point].resize(6);
									for(unsigned int k=0; k<6; k++)
									{
										out="";
										while(c!='\n')
										{
											infile.get(c);
											if(c==']'|| c== ',') break;
											out+=c;
										}
										if(out=="NOVALUE")
										{	
										goto end_of_outer_loop;
										}
										double value= std::atof((char*)out.c_str());
										if(isnan(value))
										{
											std::cout<<"##### NaN FOUND ("<<(*it)->Id()<<") ####"<<std::endl;
											goto end_of_outer_loop;
										}
										ValuesOnIntPoint[point](k)= value;
									}
									std::getline(infile, line);
								}

								(*it)->SetValueOnIntegrationPoints( INSITU_STRESS, ValuesOnIntPoint,
									model_part.GetProcessInfo());
								break;
							}
					}
					end_of_outer_loop://labeling outer loops end

					if(!element_found)
						std::cout<<"##### NO CORRESPONDING ELEMENT ("<<(*it)->Id()<<") IN INPUT FILE ####"<<std::endl;

				}
			}

			/**
			 */
			void StoreConstitutiveLawVariables( ModelPart& model_part )
			{

				ElementsArrayType& ElementsArray= model_part.Elements();

				std::string nodesfile= mDataFileName+".constitutive_law.restart";
				std::ofstream outfile;
				outfile.open((char*)nodesfile.c_str());
				outfile.precision(14);
				outfile<<"TIME["<<model_part.GetProcessInfo()[TIME]<<"]\n";
			
				for( ElementsArrayType::ptr_iterator it = ElementsArray.ptr_begin(); 
				 	it != ElementsArray.ptr_end(); ++it )
				{
					const IntegrationPointsArrayType& integration_points 
						= (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());
					outfile<<"ELEMENT["<<(*it)->Id()<<"]\n";

					std::vector<Matrix> ValuesOnIntPoint(integration_points.size());
					
					(*it)->GetValueOnIntegrationPoints(ELASTIC_LEFT_CAUCHY_GREEN_OLD, ValuesOnIntPoint, model_part.GetProcessInfo());

					if(ValuesOnIntPoint[0].size1()!=3 || ValuesOnIntPoint[0].size2()!=3 )
					{
						outfile<<"\t[NOVALUE]\n";
						continue;
					}

					for(unsigned int point=0; point< integration_points.size(); point++)
					{
						outfile<<"\t["<<ValuesOnIntPoint[point](0,0)<<","<<ValuesOnIntPoint[point](0,1)<<","<<ValuesOnIntPoint[point](0,2)<<","<<ValuesOnIntPoint[point](1,0)<<","<<ValuesOnIntPoint[point](1,1)<<","<<ValuesOnIntPoint[point](1,2)<<","<<ValuesOnIntPoint[point](2,0)<<","<<ValuesOnIntPoint[point](2,1)<<","<<ValuesOnIntPoint[point](2,2)<<"]\n";
					}
				}

				outfile.close();
			}

			/**
			 */
			void WriteConstitutiveLawVariables( ModelPart& model_part )
			{
				unsigned int id;
				char c;
				std::string line;
				std::string out;
				std::string nodesfile= mDataFileName+".constitutive_law.restart";
				std::ifstream infile;
				infile.open((char*)nodesfile.c_str());
				infile.precision(14);
				infile.seekg(0,std::ios::beg);

				ElementsArrayType& ElementsArray= model_part.Elements();

				for( ElementsArrayType::ptr_iterator it = ElementsArray.ptr_begin(); 
				 	it != ElementsArray.ptr_end(); ++it )
				{
					bool element_found= false;

					infile.seekg(0,std::ios::beg);

					const IntegrationPointsArrayType& integration_points 
						= (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

					std::vector<Matrix> ValuesOnIntPoint(integration_points.size());

					while (!infile.eof()) 
					{
							out="";
							infile.get(c);
							if(c!= '\t')
							{
								out+=c;
								while(c!='\n')
								{
									infile.get(c);
									if(c=='[') break;
									out+=c;
								}
							}
							else
							{
								std::getline(infile, line);
								continue;
							}
							if(out=="ELEMENT")
							{
								out="";
								while(c!='\n')
								{
									infile.get(c);
									if(c==']') break;
									out+=c;
								}
								id= std::atoi((char*)out.c_str());
								std::getline(infile, line);
							}
							else
							{
								std::getline(infile, line);
								continue;
							}
							if(id== (*it)->Id())
							{
								element_found= true;
								for(unsigned int point=0; point< integration_points.size(); point++)
								{
									while(c!='\n')
									{
										infile.get(c);
										if(c=='[') break;
									}
									if(ValuesOnIntPoint[point].size1()!=3|| ValuesOnIntPoint[point].size2()!=3)
									ValuesOnIntPoint[point].resize(3,3);
									for(int k=0; k<3; k++)
									{
										for(int l=0; l<3; l++)
										{
											out="";
											while(c!='\n')
											{
												infile.get(c);
												if(c==']'|| c== ',') break;
												out+=c;
											}
											if(out=="NOVALUE")
											{	
												goto end_of_outer_loop;
											}
											double value= std::atof((char*)out.c_str());
											if(isnan(value))
											{
												std::cout<<"##### NaN FOUND ("<<(*it)->Id()<<") ####"<<std::endl;
												goto end_of_outer_loop;
											}
											ValuesOnIntPoint[point](k,l)= value;
										}
									}
									std::getline(infile, line);
								}
								(*it)->SetValueOnIntegrationPoints( ELASTIC_LEFT_CAUCHY_GREEN_OLD, ValuesOnIntPoint,
									model_part.GetProcessInfo());
								break;
							}
					}
					end_of_outer_loop:

					if(!element_found)
						std::cout<<"##### NO CORRESPONDING ELEMENT ("<<(*it)->Id()<<") IN INPUT FILE ####"<<std::endl;//labeling outer loops end
				}
			}
		
			private:
				std::string mDataFileName;

	};//Class Scheme
}//namespace Kratos.

#endif /* KRATOS_VARIABLE_TRANSFER_UTILITY  defined */
