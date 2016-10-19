#ifndef PRE_UTILITES_H
#define PRE_UTILITES_H

// Project includes
#include "utilities/timer.h"
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "cluster_information.h"

namespace Kratos
{

class PreUtilities
    {
    public:

    typedef ModelPart::ElementsContainerType                         ElementsArrayType;
    typedef ModelPart::NodesContainerType::ContainerType             NodesContainerType;
    typedef WeakPointerVector<Element>                               ParticleWeakVectorType;
    typedef WeakPointerVector<Element>::iterator                     ParticleWeakIteratorType;

    KRATOS_CLASS_POINTER_DEFINITION(PreUtilities);

    PreUtilities() {}
    
    /// Default constructor
    PreUtilities(ModelPart& rModelPart)
    {
        //mInitialCenterOfMassAndMass = CalculateCenterOfMass(rModelPart);
        //mInitialMass                = CalculateTotalMass(rModelPart);
    }
    
    /// Destructor
    virtual ~PreUtilities() {}

    void SetClusterInformationInProperties(std::string const& name,
                                           boost::python::list& list_of_coordinates, 
                                           boost::python::list& list_of_radii, 
                                           double size, 
                                           double volume, 
                                           boost::python::list& inertias, 
                                           Properties::Pointer& p_properties) {
        ClusterInformation cl_info;

        cl_info.mName = name;

        array_1d<double,3> coords(3,0.0);
        
        for (int i = 0; i < boost::python::len(list_of_coordinates); i++) {
            boost::python::list list(list_of_coordinates[i]);
            coords[0] =  boost::python::extract<double>(list[0]);
            coords[1] =  boost::python::extract<double>(list[1]);
            coords[2] =  boost::python::extract<double>(list[2]);
            cl_info.mListOfCoordinates.push_back(coords);
        }
        for (int i = 0; i < boost::python::len(list_of_radii); i++) {
            cl_info.mListOfRadii.push_back(boost::python::extract<double>(list_of_radii[i]));
        }
        //TODO: check the sizes (should be the same)
        cl_info.mSize = size;
        cl_info.mVolume = volume;
        cl_info.mInertias[0] = boost::python::extract<double>(inertias[0]);
        cl_info.mInertias[1] = boost::python::extract<double>(inertias[1]);
        cl_info.mInertias[2] = boost::python::extract<double>(inertias[2]);

        p_properties->SetValue(CLUSTER_INFORMATION, cl_info);
    }        

    void CreateCartesianSpecimenMdpa(std::string filename) {
        
        // We have a prismatic specimen of dimensions 1m x 1m x 2m
        const double side = 0.15;
        int divisions;
        std::cout << "\nEnter the number of divisions: ";
        std::cin >> divisions;
        if (!divisions) {
            std::cout << "\nCannot divide by zero. Program stopped.\n\n";
            exit(EXIT_FAILURE);
        }
        const double radius = 0.5 * side / divisions;
        int node_counter = 0;
        std::vector<int> skin_nodes;
        std::vector<int> top_nodes;
        std::vector<int> bottom_nodes;
        filename += "DEM.mdpa";
        
        //
        
        std::ifstream infile(filename);
        if(infile.good()) {
            while(1){
                std::cout << "\nThe file already exists. Do you want to overwrite it? (y/n) ";
                char yn;
                std::cin >> yn;
                if(yn == 'n') {
                    std::cout << "\nStopped.\n\n";
                    exit(EXIT_FAILURE);
                }
                if(yn=='y') break;
            }
        }
        
        std::cout << "\nGenerating mesh...\n\n";
        
        clock_t initial_time, final_time;
        initial_time = clock();
        std::ofstream outputfile(filename, std::ios_base::out);
        outputfile << "Begin ModelPartData\nEnd ModelPartData\n\n";        
        outputfile << "Begin Properties 1\n";
        outputfile << "PARTICLE_DENSITY 2550.0\n";
        outputfile << "YOUNG_MODULUS 35e9\n";
        outputfile << "POISSON_RATIO 0.20\n";
        outputfile << "PARTICLE_FRICTION 0.5773502691896257\n";
        outputfile << "PARTICLE_COHESION 0.0\n";
        outputfile << "COEFFICIENT_OF_RESTITUTION 0.2\n";
        outputfile << "PARTICLE_MATERIAL 1\n";
        outputfile << "ROLLING_FRICTION 0.01\n";
        outputfile << "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME DEM_Dempack\n";
        outputfile << "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME DEM_D_Linear_viscous_Coulomb\n";
        outputfile << "SLOPE_LIMIT_COEFF_C1 24\n";
        outputfile << "SLOPE_LIMIT_COEFF_C2 28\n";
        outputfile << "SLOPE_LIMIT_COEFF_C3 1\n";
        outputfile << "SLOPE_FRACTION_N1 1\n";
        outputfile << "SLOPE_FRACTION_N2 1\n";
        outputfile << "SLOPE_FRACTION_N3 35e9\n";
        outputfile << "YOUNG_MODULUS_PLASTIC 1000\n";
        outputfile << "PLASTIC_YIELD_STRESS 0.2\n";
        outputfile << "DAMAGE_FACTOR 1\n";
        outputfile << "SHEAR_ENERGY_COEF 1\n";
        outputfile << "CONTACT_TAU_ZERO 5\n";
        outputfile << "CONTACT_SIGMA_MIN 1\n";
        outputfile << "CONTACT_INTERNAL_FRICC 20\n";
        outputfile << "End Properties\n";
        
        outputfile << "\nBegin Nodes\n";

        // Relative sizes according to axes:
        int ai=1;
        int aj=2;
        int ak=1;
        
        //Generation of the samble
        for (int k = 0; k < ai*divisions; k++) {
            for (int j = 0; j < aj* divisions; j++) {
                for (int i = 0; i < ak*divisions; i++) {
                    outputfile << ++node_counter << " " << (1 + 2 * i) * radius - 0.5*side << " " << (1 + 2 * j) * radius << " " << (1 + 2 * k) * radius - 0.5*side << '\n';
                    if ((i == 0) || (j == 0) || (k == 0) || (i == ai* divisions - 1) || (j == aj*divisions - 1) || (k ==  ak*divisions - 1)) skin_nodes.push_back(node_counter);
                    if (k == 0) bottom_nodes.push_back(node_counter);
                    if (k == 2 * divisions - 1) top_nodes.push_back(node_counter);
                }
            }
        }
        //
        outputfile << "End Nodes\n";
        outputfile << "\nBegin Elements SphericContinuumParticle3D\n";
        for (int i = 1; i <= node_counter; i++) outputfile << i << " 1 " << i << '\n';
        outputfile << "End Elements\n";
        outputfile <<  "\nBegin NodalData RADIUS\n";
        for (int i = 1; i <= node_counter; i++) outputfile << i << " 0 " << radius << '\n';
        outputfile << "End NodalData\n";
        outputfile << "\nBegin NodalData COHESIVE_GROUP // whole specimen\n";
        for (int i = 1; i <= node_counter; i++) outputfile << i << " 0 1\n";
        outputfile << "End NodalData\n";
        //outputfile << "\nBegin NodalData COHESIVE_GROUP // bottom nodes\n";
        //for (std::vector<int>::iterator it_bottom = bottom_nodes.begin(); it_bottom != bottom_nodes.end(); it_bottom++) outputfile << *it_bottom << " 0 1\n";
        //outputfile << "End NodalData\n\nBegin NodalData COHESIVE_GROUP // top nodes\n";
        //for (std::vector<int>::iterator it_top = top_nodes.begin(); it_top != top_nodes.end(); it_top++) outputfile << *it_top << " 0 1\n";
        //outputfile << "End NodalData\n";
        outputfile << "\nBegin NodalData SKIN_SPHERE\n";
        for (std::vector<int>::iterator it_skin = skin_nodes.begin(); it_skin != skin_nodes.end(); it_skin++) outputfile << *it_skin << " 0 1\n";
        outputfile << "End NodalData\n\n";
        /*outputfile << "Begin Mesh 1 // bottom nodes\n  Begin MeshData\n  VELOCITY_START_TIME 0.0\n";
        outputfile << "  FORCE_INTEGRATION_GROUP 0\n  VELOCITY_STOP_TIME 100.0\n  TOP 0\n";
        outputfile << "  IMPOSED_VELOCITY_Z_VALUE 0.0005\n  BOTTOM 0\n  End MeshData\n  Begin MeshNodes\n";
        for (std::vector<int>::iterator it_bottom = bottom_nodes.begin(); it_bottom != bottom_nodes.end(); it_bottom++) outputfile << "  " << *it_bottom << '\n';
        outputfile << "  End MeshNodes\nEnd Mesh\n\n";
        outputfile << "Begin Mesh 2 // top nodes\n  Begin MeshData\n  VELOCITY_START_TIME 0.0\n";
        outputfile << "  FORCE_INTEGRATION_GROUP 0\n  VELOCITY_STOP_TIME 100.0\n  TOP 0\n";
        outputfile << "  IMPOSED_VELOCITY_Z_VALUE -0.0005\n  BOTTOM 0\n  End MeshData\n  Begin MeshNodes\n";
        for (std::vector<int>::iterator it_top = top_nodes.begin(); it_top != top_nodes.end(); it_top++) outputfile << "  " << *it_top << '\n';
        outputfile << "  End MeshNodes\nEnd Mesh\n";*/
        outputfile.close();
        final_time = clock();
        double elapsed_time = (double(final_time) - double(initial_time)) / CLOCKS_PER_SEC;
        std::cout << "\nfinished!\n\n";
        std::cout << "\nTotal number of elements: " << node_counter << '\n';
        std::cout << "\nTime required to create the mdpa file: " << elapsed_time << " seconds\n\n";
    } 

    void MeasureTopHeight(ModelPart& rModelPart, double& subtotal, double& weight)
    {
        /*
        ElementsArrayType& pElements        = rModelPart.Elements();

        for (ElementsArrayType::iterator it= pElements.begin(); it!=pElements.end(); ++it)
        {

            if( it->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID) == 1 )
            {
                ParticleWeakVectorType& mrNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);            

                for(ParticleWeakIteratorType ineighbour = mrNeighbours.begin();  
                ineighbour != mrNeighbours.end(); ineighbour++)
                {
                    if( ineighbour->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID) != 1 )
                    {
                        subtotal += it->GetGeometry()[0].Coordinates()[1]*it->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                        weight += it->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);

                        break;
                    }
                }
            }
        }
        */
    }

    void MeasureBotHeight(ModelPart& rModelPart, double& subtotal, double& weight)
    {
        /*
            ElementsArrayType& pElements        = rModelPart.Elements();

            for (ElementsArrayType::iterator it= pElements.begin(); it!=pElements.end(); ++it)
            {
                if( it->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID) == 2 )
                {
                    ParticleWeakVectorType& mrNeighbours = it->GetValue(NEIGHBOUR_ELEMENTS);

                    for(ParticleWeakIteratorType ineighbour = mrNeighbours.begin();
                    ineighbour != mrNeighbours.end(); ineighbour++)
                    {
                        if( ineighbour->GetGeometry()[0].FastGetSolutionStepValue(GROUP_ID) != 2 )
                        {
                            subtotal += it->GetGeometry()[0].Coordinates()[1]*it->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                            weight += it->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);

                            break;
                        }
                    }
                }
            }
       */
    }

    array_1d<double, 3> GetInitialCenterOfMass()
    {
        return mInitialCenterOfMassAndMass;
    }

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.

    virtual std::string Info() const
    {
        return "";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    vector<unsigned int>& GetElementPartition() {return (mElementPartition);};

    protected:

        vector<unsigned int> mElementPartition;

    private:

        array_1d<double, 3> mInitialCenterOfMassAndMass;
        double mInitialMass;

        /// Assignment operator
        PreUtilities & operator=(PreUtilities const& rOther);

    }; // Class PreUtilities

    /// output stream function
    // 	template<std::size_t TDim>
    // 	inline std::ostream& operator << (std::ostream& rOStream)
    // 	{
    // 		rThis.PrintInfo(rOStream);
    // 		rOStream << std::endl;
    // 		rThis.PrintData(rOStream);
    //
    // 		return rOStream;
    // 	}

} // namespace Kratos

#endif // PRE_UTILITES_H
