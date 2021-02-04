#ifndef PRE_UTILITES_H
#define PRE_UTILITES_H

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
#include "utilities/timer.h"
#include "includes/variables.h"
#include "utilities/openmp_utils.h"
#include "cluster_information.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos
{

class PreUtilities
    {
    public:

    typedef ModelPart::ElementsContainerType                         ElementsArrayType;
    typedef ModelPart::NodesContainerType::ContainerType             NodesContainerType;
    typedef GlobalPointersVector<Element>                               ParticleWeakVectorType;
    typedef GlobalPointersVector<Element>::iterator                     ParticleWeakIteratorType;

    KRATOS_CLASS_POINTER_DEFINITION(PreUtilities);

    /// Default constructor
    PreUtilities() {}

    PreUtilities(ModelPart& rModelPart)
    {
        //mInitialCenterOfMassAndMass = CalculateCenterOfMass(rModelPart);
        //mInitialMass                = CalculateTotalMass(rModelPart);
    }

    /// Destructor
    virtual ~PreUtilities() {}

    void SetClusterInformationInProperties(std::string const& name,
                                           pybind11::list& list_of_coordinates,
                                           pybind11::list& list_of_radii,
                                           double size,
                                           double volume,
                                           pybind11::list& inertias,
                                           Properties::Pointer& p_properties) {
        ClusterInformation cl_info;

        cl_info.mName = name;

        array_1d<double,3> coords(3,0.0);

        for (int i = 0; i < (int)pybind11::len(list_of_coordinates); i++) {
            pybind11::list list(list_of_coordinates[i]);
            coords[0] =  pybind11::cast<double>(list[0]);
            coords[1] =  pybind11::cast<double>(list[1]);
            coords[2] =  pybind11::cast<double>(list[2]);
            cl_info.mListOfCoordinates.push_back(coords);
        }
        for (int i = 0; i < (int)pybind11::len(list_of_radii); i++) {
            cl_info.mListOfRadii.push_back(pybind11::cast<double>(list_of_radii[i]));
        }
        //TODO: check the sizes (should be the same)
        cl_info.mSize = size;
        cl_info.mVolume = volume;
        cl_info.mInertias[0] = pybind11::cast<double>(inertias[0]);
        cl_info.mInertias[1] = pybind11::cast<double>(inertias[1]);
        cl_info.mInertias[2] = pybind11::cast<double>(inertias[2]);

        p_properties->SetValue(CLUSTER_INFORMATION, cl_info);
    }

    void PrintNumberOfNeighboursHistogram(const ModelPart& rSpheresModelPart, std::string const& filename) {
        std::vector<int> number_of_spheres_with_i_neighbours;
        number_of_spheres_with_i_neighbours.resize(20);
        for(int i=0; i<(int)number_of_spheres_with_i_neighbours.size(); i++) {number_of_spheres_with_i_neighbours[i] = 0;}

        const ElementsArrayType& pElements = rSpheresModelPart.GetCommunicator().LocalMesh().Elements();
        ElementsArrayType::ptr_const_iterator begin = pElements.ptr_begin();

        for(int i=0; i<(int)pElements.size(); i++) {
            ElementsArrayType::ptr_const_iterator it = begin + i;
            const Element& el = **it;
            const SphericContinuumParticle* p_cont_sphere = dynamic_cast<const SphericContinuumParticle*>(&el);
            if(p_cont_sphere) {
                unsigned int size = p_cont_sphere->mContinuumInitialNeighborsSize;
                if(size > number_of_spheres_with_i_neighbours.size() - 1) size = number_of_spheres_with_i_neighbours.size() - 1;
                number_of_spheres_with_i_neighbours[size] += 1;
            } else {
                const SphericParticle* p_sphere = dynamic_cast<const SphericParticle*>(&el);
                unsigned int size = p_sphere->mNeighbourElements.size();
                if(size > number_of_spheres_with_i_neighbours.size() - 1) size = number_of_spheres_with_i_neighbours.size() - 1;
                number_of_spheres_with_i_neighbours[size] += 1;
            }
        }
        std::ofstream outputfile(filename, std::ios_base::out | std::ios_base::app);
        outputfile << "number_of_neighbours   percentage_of_spheres_with_that_number_of_neighbours    number_of_spheres_with_that_number_of_neighbours\n";
        for(int i=0; i<(int)number_of_spheres_with_i_neighbours.size(); i++) {
            const double percentage = (double)(number_of_spheres_with_i_neighbours[i]) / (double)(rSpheresModelPart.NumberOfElements(0)) * 100.0;
            outputfile <<i<<"        "<<percentage<<"        "<<number_of_spheres_with_i_neighbours[i]<<"\n";
        }

    }


    void FillAnalyticSubModelPartUtility(ModelPart& rSpheresModelPart, ModelPart& rAnalyticSpheresModelPart){
        ElementsArrayType& pElements = rSpheresModelPart.GetCommunicator().LocalMesh().Elements();
        std::vector<std::vector<std::size_t> > thread_vectors_of_ids;
        int mNumberOfThreads = ParallelUtilities::GetNumThreads();
        thread_vectors_of_ids.resize(mNumberOfThreads);

        #pragma omp parallel for
        for (int k = 0; k < (int)pElements.size(); k++) {
            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            int analytic_particle_id = it->Id();
            thread_vectors_of_ids[OpenMPUtils::ThisThread()].push_back(analytic_particle_id);
        }
        std::vector<std::size_t> vector_of_ids;
        for (int i = 0; i < mNumberOfThreads; i++) {
            vector_of_ids.insert(vector_of_ids.end(), thread_vectors_of_ids[i].begin(), thread_vectors_of_ids[i].end());
        }
        rAnalyticSpheresModelPart.AddElements(vector_of_ids);
    }

//    non-OMP version
//    void FillAnalyticSubModelPartUtility(ModelPart& rSpheresModelPart, ModelPart& rAnalyticSpheresModelPart){
//        ElementsArrayType& pElements = rSpheresModelPart.GetCommunicator().LocalMesh().Elements();
//        std::vector<long unsigned int> vector_of_ids;
//        for (int k = 0; k < (int)pElements.size(); k++) {
//            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
//            int analytic_particle_id = it->Id();
//            vector_of_ids.push_back(analytic_particle_id);
//        }
//        rAnalyticSpheresModelPart.AddElements(vector_of_ids);
//    }

    void ResetSkinParticles(ModelPart& r_model_part) {
        auto& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();
        #pragma omp parallel for
        for (int k = 0; k < (int)pNodes.size(); k++) {
            auto it = pNodes.begin() + k;
            it->FastGetSolutionStepValue(SKIN_SPHERE) = 0.0;
        }
    }

    void SetSkinParticlesInnerCircularBoundary(ModelPart& r_model_part, const double inner_radius, const double detection_radius) {
        auto& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();

        #pragma omp parallel for
        for (int k = 0; k < (int)pNodes.size(); k++) {
            auto it = pNodes.begin() + k;
            const array_1d<double, 3>& coords = it->Coordinates();
            array_1d<double, 3> vector_distance_to_center;
            noalias(vector_distance_to_center) = coords;
            const double distance_to_center = MathUtils<double>::Norm3(vector_distance_to_center);
            if(distance_to_center < inner_radius + detection_radius) {
                it->FastGetSolutionStepValue(SKIN_SPHERE) = 1.0;
            }
        }
    }

    void SetSkinParticlesOuterCircularBoundary(ModelPart& r_model_part, const double outer_radius, const double detection_radius) {
        auto& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();

        #pragma omp parallel for
        for (int k = 0; k < (int)pNodes.size(); k++) {
            auto it = pNodes.begin() + k;
            const array_1d<double, 3>& coords = it->Coordinates();
            array_1d<double, 3> vector_distance_to_center;
            noalias(vector_distance_to_center) = coords;
            const double distance_to_center = MathUtils<double>::Norm3(vector_distance_to_center);
            const double radius = it->FastGetSolutionStepValue(RADIUS);
            if (distance_to_center + radius > outer_radius - detection_radius) {
                it->FastGetSolutionStepValue(SKIN_SPHERE) = 1.0;
            }
        }
    }

    void SetSkinParticlesOuterSquaredBoundary(ModelPart& r_model_part, const double outer_radius, const array_1d<double, 3>& center, const double detection_radius) {

        auto& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();

        #pragma omp parallel for
        for (int k = 0; k < (int)pNodes.size(); k++) {
            auto it = pNodes.begin() + k;
            const array_1d<double, 3>& coords = it->Coordinates();
            array_1d<double, 3> vector_distance_to_center;
            noalias(vector_distance_to_center) = coords - center;
            const double total_x_distance = fabs(vector_distance_to_center[0]);
            const double total_y_distance = fabs(vector_distance_to_center[1]);
            const double radius = it->FastGetSolutionStepValue(RADIUS);

            if ((total_x_distance + radius > outer_radius - detection_radius) || (total_y_distance + radius > outer_radius - detection_radius)) {
                it->FastGetSolutionStepValue(SKIN_SPHERE) = 1.0;
            }
        }
    }

    void BreakBondUtility(ModelPart& rSpheresModelPart) {

        ElementsArrayType& pElements = rSpheresModelPart.GetCommunicator().LocalMesh().Elements();
        #pragma omp parallel for
        for (int k = 0; k < (int)pElements.size(); k++) {

            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            Element* p_element = &(*it);
            SphericContinuumParticle* p_sphere = dynamic_cast<SphericContinuumParticle*>(p_element);

            if (p_sphere->mNeighbourElements[k] == NULL) continue;

            double x_node = p_sphere->GetGeometry()[0].Coordinates()[0];
            double y_node = p_sphere->GetGeometry()[0].Coordinates()[1];
            double z_node = p_sphere->GetGeometry()[0].Coordinates()[2];
            double radius = 0.0225; // radi

            if ((x_node*x_node + z_node*z_node >= radius*radius && y_node < 0.01) || (x_node*x_node + z_node*z_node >= radius*radius && y_node > 0.07)) {   // 1- geometry condition
                unsigned int number_of_neighbors = p_sphere->mContinuumInitialNeighborsSize;
                for (unsigned int i = 0; i < number_of_neighbors; i++)
                {
                    SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(p_sphere->mNeighbourElements[i]);
                    double x_node_it = neighbour_iterator->GetGeometry()[0].Coordinates()[0];
                    double z_node_it = neighbour_iterator->GetGeometry()[0].Coordinates()[2];
                    double radius_it = 0.0225; // radi de la entalla en el shear test.
                    if (x_node_it*x_node_it + z_node_it*z_node_it < radius_it*radius_it) {   // 2- geometry condition

                        //int& failure_type = p_sphere->mIniNeighbourFailureId[i];
                        //failure_type = 1;
                        p_sphere->Set(TO_ERASE, true);
                        neighbour_iterator->Set(TO_ERASE, true);

                        //noalias(other_to_me_vector)         = p_sphere->GetGeometry()[0].Coordinates() - p_sphere->mNeighbourElements[i]->GetGeometry()[0].Coordinates();
                        //noalias(initial_other_to_me_vector) = p_sphere->GetGeometry()[0].GetInitialPosition() - p_sphere->mNeighbourElements[i]->GetGeometry()[0].GetInitialPosition();
                    }
                }
            } else if ((x_node*x_node + z_node*z_node < radius*radius && y_node < 0.01) || (x_node*x_node + z_node*z_node < radius*radius && y_node > 0.07)) {
                unsigned int number_of_neighbors = p_sphere->mContinuumInitialNeighborsSize;
                for (unsigned int i = 0; i < number_of_neighbors; i++)
                {
                    SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(p_sphere->mNeighbourElements[i]);
                    double x_node_it = neighbour_iterator->GetGeometry()[0].Coordinates()[0];
                    double z_node_it = neighbour_iterator->GetGeometry()[0].Coordinates()[2];
                    double radius_it = 0.0225; // radi de la entalla en el shear test.
                    if (x_node_it*x_node_it + z_node_it*z_node_it > radius_it*radius_it) {   // 2- geometry condition

                        //int& failure_type = p_sphere->mIniNeighbourFailureId[i];
                        //failure_type = 1;
                        p_sphere->Set(TO_ERASE, true);
                        neighbour_iterator->Set(TO_ERASE, true);
                    }
                }
            }
        }
    }

    void CreateCartesianSpecimenMdpa(std::string filename) {

        // We have a prismatic specimen of dimensions 1m x 1m x 2m
        const double side = 0.15;
        int divisions;
        KRATOS_WARNING("DEM") << "\nEnter the number of divisions: ";
        std::cin >> divisions;
        if (!divisions) {
            KRATOS_WARNING("DEM") << "\nCannot divide by zero. Program stopped.\n\n";
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
                KRATOS_WARNING("DEM") << "\nThe file already exists. Do you want to overwrite it? (y/n) ";
                char yn;
                std::cin >> yn;
                if(yn == 'n') {
                    KRATOS_WARNING("DEM") << "\nStopped.\n\n";
                    exit(EXIT_FAILURE);
                }
                if(yn=='y') break;
            }
        }

        KRATOS_INFO("DEM") << "\nGenerating mesh...\n\n";

        clock_t initial_time, end_time;
        initial_time = clock();
        std::ofstream outputfile(filename, std::ios_base::out);
        outputfile << "Begin ModelPartData\nEnd ModelPartData\n\n";
        outputfile << "Begin Properties 1\n";
        outputfile << "PARTICLE_DENSITY 2550.0\n";
        outputfile << "YOUNG_MODULUS 35e9\n";
        outputfile << "POISSON_RATIO 0.20\n";
        outputfile << "STATIC_FRICTION 0.5773502691896257\n";
        outputfile << "DYNAMIC_FRICTION 0.5773502691896257\n";
        outputfile << "PARTICLE_COHESION 0.0\n";
        outputfile << "COEFFICIENT_OF_RESTITUTION 0.2\n";
        outputfile << "PARTICLE_MATERIAL 1\n";
        outputfile << "ROLLING_FRICTION 0.01\n";
        outputfile << "ROLLING_FRICTION_WITH_WALLS 0.01\n";
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
        end_time = clock();
        double elapsed_time = (double(end_time) - double(initial_time)) / CLOCKS_PER_SEC;
        KRATOS_INFO("DEM") << "\nfinished!\n\n";
        KRATOS_INFO("DEM") << "\nTotal number of elements: " << node_counter << '\n';
        KRATOS_INFO("DEM") << "\nTime required to create the mdpa file: " << elapsed_time << " seconds\n\n";
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

    void MarkToEraseParticlesOutsideRadius(ModelPart& r_model_part, const double max_radius, const array_1d<double, 3>& center, const double tolerance_for_erasing) {
        auto& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();

        #pragma omp parallel for
        for (int k = 0; k < (int)pNodes.size(); k++) {
            auto it = pNodes.begin() + k;
            const array_1d<double, 3>& coords = it->Coordinates();
            array_1d<double, 3> vector_distance_to_center;
            noalias(vector_distance_to_center) = coords - center;
            const double distance_to_center = MathUtils<double>::Norm3(vector_distance_to_center);
            const double radius = it->FastGetSolutionStepValue(RADIUS);
            if(distance_to_center + radius > max_radius + tolerance_for_erasing) {
                it->Set(TO_ERASE, true);
            }
        }
    }

    void ApplyConcentricForceOnParticles(ModelPart& r_model_part, const array_1d<double, 3>& center, const double density_for_artificial_gravity) {
        auto& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        #pragma omp parallel for
        for (int k = 0; k < (int)pElements.size(); k++) {
            auto it = pElements.begin() + k;
            auto& node = it->GetGeometry()[0];
            const array_1d<double, 3>& coords = node.Coordinates();
            array_1d<double, 3> vector_particle_to_center;
            noalias(vector_particle_to_center) = center - coords;
            const double distance_to_center = MathUtils<double>::Norm3(vector_particle_to_center);
            const double inv_dist = 1.0 / distance_to_center;
            array_1d<double, 3> force;
            SphericParticle* spheric_p_particle = dynamic_cast<SphericParticle*> (&*it);
            const double volume = spheric_p_particle->CalculateVolume();
            noalias(force) = inv_dist * vector_particle_to_center * volume * density_for_artificial_gravity;
            node.FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) = force;
        }
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

    std::vector<unsigned int>& GetElementPartition() {return (mElementPartition);};

    protected:

        std::vector<unsigned int> mElementPartition;

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
