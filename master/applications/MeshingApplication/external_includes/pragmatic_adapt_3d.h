#if !defined(KRATOS_PRAGMATIC_ADAPTOR_H_INCLUDED )
#define  KRATOS_PRAGMATIC_ADAPTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cstdlib>
#include <boost/timer.hpp>



#include "tetgen.h" // Defined tetgenio, tetrahedralize().

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


#include "include/Mesh.h"
//#include "include/VTKTools.h"
#include "include/MetricField.h"

#include "include/Coarsen.h"
#include "include/Refine.h"
#include "include/Smooth.h"
#include "include/Swapping.h"
#include "include/ticker.h"

#include <mpi.h>


namespace Kratos
{

class PragmaticAdaptor
{
public:
//     void cout_quality(const Mesh<double> *mesh, std::string operation)
//     {
//         double qmean = mesh->get_qmean();
//         double qmin = mesh->get_qmin();
// 
//         int rank;
//         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// 
//         if(rank==0)
//             std::cout<<operation<<": step in quality (mean, min): ("<<qmean<<", "<<qmin<<")"<<std::endl;
//     }

    void AdaptMesh(ModelPart& rmodel_part)
    {

//         int required_thread_support=MPI_THREAD_SINGLE;
//         int provided_thread_support;
//         MPI_Init_thread(&argc, &argv, required_thread_support, &provided_thread_support);
//         assert(required_thread_support==provided_thread_support);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        bool verbose = true;

        // Benchmark times.
        double time_coarsen=0, time_refine=0, time_swap=0, time_smooth=0, time_adapt=0, tic;

        //**************************************
        //here i must import mesh from Kratos

        mesh->create_boundary();


        //**************************************
        //here i shall assign the target metric field

        MetricField<double,3> metric_field(*mesh);

        size_t NNodes = mesh->get_number_nodes();
        double eta=0.05;

        for(size_t i=0; i<NNodes; i++)
        {
            double target_h = 1.0; //TODO: read target H from the nodes!!...
            double d = 1/sqrt(h);
            double m[] = {d,0.0,0.0,d, 0.0, d};

//             double x = 2*mesh->get_coords(i)[0] - 1;
//             double y = 2*mesh->get_coords(i)[1];
// 
//             double m[] = {0.2*(-8*x + 4*sin(5*y))/pow(pow(2*x - sin(5*y), 2) + 0.01, 2) - 250.0*sin(50*x),  2.0*(2*x - sin(5*y))*cos(5*y)/pow(pow(2*x - sin(5*y), 2) + 0.01, 2),                                                        0,
//                           -5.0*(2*x - sin(5*y))*pow(cos(5*y), 2)/pow(pow(2*x - sin(5*y), 2) + 0.01, 2) + 2.5*sin(5*y)/(pow(2*x - sin(5*y), 2) + 0.01), 0,
//                           0
//                          };

/*
            for(int j=0; j<5; j++)
                m[j]/=eta;
            m[5] = 1.0;*/

            metric_field.set_metric(m, i);
        }
        metric_field.apply_max_aspect_ratio(10);
        metric_field.update_mesh();

//  VTKTools<double>::export_vtu("../data/test_adapt_3d-initial", mesh);

//         cout_quality(mesh, "Initial quality");

        // See Eqn 7; X Li et al, Comp Methods Appl Mech Engrg 194 (2005) 4915-4950
        double L_up = sqrt(2.0);
        double L_low = L_up/2;

        Coarsen<double, 3> coarsen(*mesh);
        Smooth<double, 3> smooth(*mesh);
        Refine<double, 3> refine(*mesh);
        Swapping<double, 3> swapping(*mesh);

        time_adapt = get_wtime();

        double L_max = mesh->maximal_edge_length();
        double alpha = sqrt(2.0)/2;
        double L_ref = std::max(alpha*L_max, L_up);

        if(verbose)
            std::cout<<"Phase I\n";

        for(size_t i=0; i<10; i++)
        {
            // Coarsen
            tic = get_wtime();
            coarsen.coarsen(L_low, L_ref, true);
            time_coarsen += get_wtime() - tic;

//             if(verbose)
//                 cout_quality(mesh, "Coarsen");

            // Refine
            tic = get_wtime();
            refine.refine(L_ref);
            time_refine += get_wtime() - tic;

//             if(verbose)
//                 cout_quality(mesh, "refine");

            // Swap
            tic = get_wtime();
            swapping.swap(0.1);
            time_swap += get_wtime() - tic;

//             if(verbose)
//                 cout_quality(mesh, "Swap");

            // Smooth
            tic = get_wtime();
            smooth.smooth(1);
            time_smooth += get_wtime()-tic;

//             if(verbose)
//                 cout_quality(mesh, "Smooth");

            alpha = (1.0-1e-2*i*i)*sqrt(2.0)/2;
            L_max = mesh->maximal_edge_length();
            L_ref = std::max(alpha*L_max, L_up);


            if(L_max>1.0 and (L_max-L_up)<0.01)
                break;
        }

        if(verbose)
            std::cout<<"Phase II\n";

        for(size_t i=0; i<5; i++)
        {
            tic = get_wtime();
            coarsen.coarsen(L_up, L_up, true);
            time_coarsen += get_wtime() - tic;
//             if(verbose)
//                 cout_quality(mesh, "coarsen");

            tic = get_wtime();
            swapping.swap(0.7);
//             if(verbose)
//                 cout_quality(mesh, "Swap");
            time_swap += get_wtime() - tic;

            tic = get_wtime();
            smooth.smooth(1);
//             if(verbose)
//                 cout_quality(mesh, "Smooth");
            time_smooth += get_wtime()-tic;
        }

        double time_defrag = get_wtime();
        mesh->defragment();
        time_defrag = get_wtime()-time_defrag;

        if(verbose)
        {
            mesh->verify();

            VTKTools<double>::export_vtu("../data/test_adapt_3d-basic", mesh);
        }

        tic = get_wtime();
        smooth.smooth(20);
        time_smooth += get_wtime()-tic;

//         if(verbose)
//             cout_quality(mesh, "Final smooth");

        time_adapt = get_wtime()-time_adapt;

        if(verbose)
        {
            if(rank==0)
                std::cout<<"After optimisation based smoothing:\n";
            mesh->verify();
        }

        VTKTools<double>::export_vtu("../data/test_adapt_3d", mesh);

        double qmean = mesh->get_qmean();
        double qmin = mesh->get_qmin();

        long double volume = mesh->calculate_volume();
        long double area = mesh->calculate_area();

        delete mesh;

        if(rank==0)
        {
            std::cout<<"BENCHMARK: time_coarsen time_refine time_swap time_smooth time_defrag time_adapt time_other\n";
            double time_other = (time_adapt-(time_coarsen+time_refine+time_swap+time_smooth+time_defrag));
            std::cout<<"BENCHMARK: "
                     <<std::setw(12)<<time_coarsen<<" "
                     <<std::setw(11)<<time_refine<<" "
                     <<std::setw(9)<<time_swap<<" "
                     <<std::setw(11)<<time_smooth<<" "
                     <<std::setw(11)<<time_defrag<<" "
                     <<std::setw(10)<<time_adapt<<" "
                     <<std::setw(10)<<time_other<<"\n";

            std::cout<<"Expecting qmean>0.7, qmin>0.2: ";
            if((qmean>0.8)&&(qmin>0.2))
                std::cout<<"pass"<<std::endl;
            else
                std::cout<<"fail (qmean="<<qmean<<", qmin="<<qmin<<")"<<std::endl;

//             std::cout<<"Expecting volume == 1: ";
//             if(fabs(volume-1)<DBL_EPSILON)
//                 std::cout<<"pass"<<std::endl;
//             else
//                 std::cout<<"fail (volume="<<volume<<")"<<std::endl;
// 
//             std::cout<<"Expecting area == 6: ";
//             if(fabs(area-6)<DBL_EPSILON)
//                 std::cout<<"pass"<<std::endl;
//             else
//                 std::cout<<"fail (area="<<area<<")"<<std::endl;
        }

        return 0;
    }
};

}

#endif // KRATOS_PRAGMATIC_ADAPTOR_H_INCLUDED  defined 