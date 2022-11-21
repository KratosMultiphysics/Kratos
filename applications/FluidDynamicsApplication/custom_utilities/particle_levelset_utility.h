//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include "utilities/binbased_fast_point_locator.h"

// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// implements the Particle Levelset method as reviewed in https://doi.org/10.1016/j.cma.2015.05.017.
/** This classes uses particles to improve tracking of the levelset position
*/
template<unsigned int TDim>
class ParticleLevelsetUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParticleLevelsetUtility
    KRATOS_CLASS_POINTER_DEFINITION(ParticleLevelsetUtility<TDim>);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParticleLevelsetUtility(
        ModelPart& rVolumeModelPart,
        ModelPart& rParticlesModelPart,
        const Variable<double>& rDistanceVar,
        const double DistanceAroundZero
    ):
        mrVolumeModelPart(rVolumeModelPart),
        mrParticlesModelPart(rParticlesModelPart),
        mrDistanceVar(rDistanceVar),
        mDistanceAroundZero(DistanceAroundZero)
    {
    }

    /// Destructor.
    virtual ~ParticleLevelsetUtility() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**this function is to be called before convecting the distance field by eulerian methods.
    It seeds particles around the zero of the levelset*/
    void Prepare()
    {
        //erase all nodes from rParticlesModelPart (note that we assume that they are created by this utility, so
        //we will keep track of them internally to avoid generating unnecessary constructions)
        //note also that we take shortcuts instead of using the RemoveNodes functions of the modelpart
        mrParticlesModelPart.Nodes() = ModelPart::NodesContainerType();

        //identify where particles should be seeded within rVolumeModelPart
        std::vector<Element*> active_elements;
        active_elements.reserve(10000);
        for(auto& rElem : mrVolumeModelPart.Elements())
        {
            for(const auto& rNode : rElem.GetGeometry())
                if(rNode.FastGetSolutionStepValue(mrDistanceVar) < mDistanceAroundZero)
                {
                    active_elements.push_back(&rElem);
                    break;
                }
        }

        //ensure that enough nodes exist, otherwise create new
        constexpr unsigned int particles_per_element = ParticlesPerElement<TDim>();
        const unsigned int particles_needed = particles_per_element*active_elements.size();
        GenerateParticles(particles_needed);

        //assign to all of the particles a correct position and initial distance
        BoundedMatrix<double, particles_per_element, 3> positions;
        BoundedMatrix<double, particles_per_element, TDim+1> N;
        auto it = mParticles.begin();
        for(auto pelem : active_elements)
        {

            const auto& geom = pelem->GetGeometry();
            ReseedPositions(geom,positions,N);

            for(unsigned int i=0; i<particles_per_element; ++i)
            {
                auto& pparticle = (*it);

                pparticle->Coordinates() = row(positions,i);
                pparticle->GetInitialPosition().Coordinates() = pparticle->Coordinates();
                (pparticle->FastGetSolutionStepValue(DISPLACEMENT)).clear();
                const double d = InterpolateHistoricalVariable(geom, mrDistanceVar, row(N,i));
                AssignParticleRadius(pparticle, d);
                it++; //go to the next particle
            }
        }

        //assign nodes to rParticlesModelPart
        for(unsigned int i=0; i<particles_needed; ++i)
        {
            auto& pparticle = mParticles[i];
            mrParticlesModelPart.Nodes().push_back(pparticle);
        }

    }

    /**this function moves the particles in rParticlesModelPart and use their values to correct the
    distance field within the rVolumeModelPart*/
    void MoveAndCorrect(BinBasedFastPointLocator<TDim>& rLocator)
    {
        //move the particles by RK4

        //generate a list of "escaped particles"
        std::vector< std::pair<Node<3>*, Element::Pointer> > escaped_particles;
        unsigned int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        double SearchTolerance =1e-6;

        for(auto& rParticle : mrParticlesModelPart.Nodes())
        {
            Vector shape_functions;
            Element::Pointer p_element;
            const bool is_found = rLocator.FindPointOnMesh(rParticle.Coordinates(), shape_functions, p_element, results.begin(), results.size(), SearchTolerance);

            if(is_found)
            {
                auto& r_geom = p_element->GetGeometry();
                const double eulerian_distance = InterpolateHistoricalVariable(r_geom, mrDistanceVar, shape_functions);

                const double d = rParticle.FastGetSolutionStepValue(mrDistanceVar);
                const double rp = std::abs(d);

                bool same_side = (std::signbit(eulerian_distance)==std::signbit(d));
                bool error_is_large = (rp < std::abs(eulerian_distance));

                if( !same_side && error_is_large ) //particle on the wrong side of the interface
                {
                    std::pair<Node<3>*, Element::Pointer> ppair = std::make_pair(&rParticle, p_element);
                    escaped_particles.push_back(ppair);
                }
            }
            else
            {
                std::cout << "particle not found !! " << rParticle << std::endl;
                KRATOS_WATCH(rParticle.Coordinates());
                KRATOS_WATCH(rParticle.Id())
            }
        }

        KRATOS_WATCH(escaped_particles.size())

        //initialize the phi_plus and phi_minus for all nodes around escaped particles
        std::map<Node<3>::IndexType, double> phi_plus_map;
        std::map<Node<3>::IndexType, double> phi_minus_map;
        for(auto& item : escaped_particles)
        {
            auto& pelem = item.second;
            for(const auto& rNode : pelem->GetGeometry())
            {
                double d = rNode.FastGetSolutionStepValue(DISTANCE);
                phi_plus_map[rNode.Id()]  = d;
                phi_minus_map[rNode.Id()] = d;
//std::cout << rNode.Id() << " " << phi_plus_map[rNode.Id()] << std::endl;
            }
        }


        //correct the levelset
        for(auto& item : escaped_particles)
        {
            auto& pparticle = item.first;
            auto& pelem = item.second;

            PerformLocalCorrection(pparticle,pelem,phi_plus_map,phi_minus_map);

        }

        //perform final correction
        auto it_plus = phi_plus_map.begin();
        auto it_minus = phi_minus_map.begin();
        for(unsigned int i=0; i<phi_plus_map.size(); ++i)
        {
            double phi_plus = it_plus->second;
            double phi_minus = it_minus->second;
            auto& r_node = mrVolumeModelPart.Nodes()[it_plus->first];

            if(std::abs(phi_plus) > std::abs(phi_minus))
                r_node.FastGetSolutionStepValue(DISTANCE) = phi_minus;
            else
                r_node.FastGetSolutionStepValue(DISTANCE) = phi_plus;

            it_plus++;
            it_minus++;
        }

    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    // /// Turn back information as a string.
    // virtual std::string Info() const {};

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ParticleLevelSetUtility" << std::endl;
    };

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {

    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{
    //*****************************************************************
    void GenerateParticles(unsigned int particles_needed)
    {
        if(particles_needed>mParticles.size())
        {
            //ensure that ids are numbered consecutively in mParticles
            for(unsigned int i=0; i<mParticles.size(); ++i)
                mParticles[i]->SetId(i+1);

            mParticles.reserve(particles_needed);

            auto buffer_size  = mrParticlesModelPart.GetBufferSize();
            auto pvar_list = mrParticlesModelPart.pGetNodalSolutionStepVariablesList();

            for(unsigned int i=mParticles.size(); i<particles_needed; ++i)
            {
                //create a new node
                auto p_new_node = Kratos::make_intrusive< Node<3> >( i+1, 0.0,0.0,0.0 );

                // Giving model part's variables list to the node
                p_new_node->SetSolutionStepVariablesList(pvar_list);

                //set buffer size
                p_new_node->SetBufferSize(buffer_size);

                mParticles.push_back(p_new_node);
            }
        }
KRATOS_WATCH(mParticles.size())
    }

    //*****************************************************************
    void AssignParticleRadius(Node<3>::Pointer pparticle, const double d)
    {


        double rp = std::abs(d);
        double rmax = 0.3*mDistanceAroundZero;
        double rmin = 0.2*rmax; //0.2
        if(rp > rmax) rp = rmax;
        else if(rp < rmin) rp = rmin;
        double sp = (d>0) ? 1.0 : -1.0;
        double assigned_d = sp*rp;
        pparticle->FastGetSolutionStepValue(mrDistanceVar) = assigned_d;
    }


    //*****************************************************************
    void PerformLocalCorrection
        (Node<3>* pParticle,
        Element::Pointer pElement,
        std::map<Node<3>::IndexType, double>& rPhiPlus,
        std::map<Node<3>::IndexType, double>& rPhiMinus)
    {
        double d = pParticle->FastGetSolutionStepValue(mrDistanceVar);
        double sp = (d>0) ? 1.0 : -1.0;
        double rp = std::abs(d);
        auto& particle_coords = pParticle->Coordinates();

            for(auto& rNode : pElement->GetGeometry())
            {
                double phi_p = sp*(rp - norm_2(rNode.Coordinates() - particle_coords));
                if(d>0)
                {
KRATOS_WATCH(rPhiPlus.find(rNode.Id())->second)
                    double& phi_plus = rPhiPlus.find(rNode.Id())->second;
                    phi_plus = std::max(phi_plus, phi_p);
KRATOS_WATCH(rPhiPlus.find(rNode.Id())->second)
                }
                else
                {
                    double& phi_minus = rPhiMinus.find(rNode.Id())->second;
double old_value = phi_minus;

                    phi_minus = std::min(phi_minus, phi_p);
if(phi_minus < old_value)
{
    KRATOS_WATCH(rNode.Id())
    KRATOS_WATCH(old_value)
    KRATOS_WATCH(phi_minus)
    KRATOS_WATCH(rPhiMinus.find(rNode.Id())->second)
}
                }
            }
    }

    inline double InterpolateHistoricalVariable(
        const Geometry<Node<3>>& rGeom,
        const Variable<double>& rVar,
        const Vector& rN) noexcept
    {
        double res = rN[0]*rGeom[0].FastGetSolutionStepValue(rVar);
        for(unsigned int i=1; i<rN.size(); ++i)
            res += rN[i]*rGeom[i].FastGetSolutionStepValue(rVar);
        return res;
    }

    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrVolumeModelPart;
    ModelPart& mrParticlesModelPart;
    const Variable<double>& mrDistanceVar;
    const double mDistanceAroundZero;
    std::vector<Node<3>::Pointer> mParticles;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    //*****************************************************************
    template<unsigned int Dim>
    static constexpr unsigned int ParticlesPerElement()
    {
        if constexpr (Dim==2)
            return 20;
        if constexpr (Dim==3)
            return 9;
    }

    //*****************************************************************
    //2D
    void ReseedPositions(const Geometry< Node < 3 > >& rGeom,
                         BoundedMatrix<double, 20, 3 > & rPositions,
                         BoundedMatrix<double, 20, 3 > & rN)
    {
        rN( 0 ,0)= 0.20922092720866203 ; rN( 0 ,1)= 0.058036486618220806 ;
        rN( 1 ,0)= 0.1294716987758875 ; rN( 1 ,1)= 0.4986436367034912 ;
        rN( 2 ,0)= 0.5541051775217056 ; rN( 2 ,1)= 0.42883730959147215 ;
        rN( 3 ,0)= 0.5994162689894438 ; rN( 3 ,1)= 0.11212435364723206 ;
        rN( 4 ,0)= 0.36456844210624695 ; rN( 4 ,1)= 0.28130291029810905 ;
        rN( 5 ,0)= 0.2890447061508894 ; rN( 5 ,1)= 0.22459513042122126 ;
        rN( 6 ,0)= 0.10543740168213844 ; rN( 6 ,1)= 0.8317829854786396 ;
        rN( 7 ,0)= 0.056022411212325096 ; rN( 7 ,1)= 0.6467874338850379 ;
        rN( 8 ,0)= 0.07221280224621296 ; rN( 8 ,1)= 0.46329878829419613 ;
        rN( 9 ,0)= 0.02670334279537201 ; rN( 9 ,1)= 0.030866234563291073 ;
        rN( 10 ,0)= 0.6647166814655066 ; rN( 10 ,1)= 0.08413012884557247 ;
        rN( 11 ,0)= 0.2558545283973217 ; rN( 11 ,1)= 0.6056225309148431 ;
        rN( 12 ,0)= 0.47911561094224453 ; rN( 12 ,1)= 0.18851773906499147 ;
        rN( 13 ,0)= 0.4298989810049534 ; rN( 13 ,1)= 0.2548649590462446 ;
        rN( 14 ,0)= 0.23850986175239086 ; rN( 14 ,1)= 0.6745376391336322 ;
        rN( 15 ,0)= 0.8373343050479889 ; rN( 15 ,1)= 0.13351291231811047 ;
        rN( 16 ,0)= 0.12786749750375748 ; rN( 16 ,1)= 0.37275556940585375 ;
        rN( 17 ,0)= 0.20339027605950832 ; rN( 17 ,1)= 0.17994011007249355 ;
        rN( 18 ,0)= 0.597377423197031 ; rN( 18 ,1)= 0.236325872130692 ;
        rN( 19 ,0)= 0.547961475327611 ; rN( 19 ,1)= 0.3008382674306631 ;

        for(unsigned int i=0; i<rN.size1(); ++i)
        {
            rN(i,2) = 1.0 - rN(i,0) - rN(i,1);
            rPositions(i,0) = rN(i,0)*rGeom[0].X() + rN(i,1)*rGeom[1].X() + rN(i,2)*rGeom[2].X();
            rPositions(i,1) = rN(i,0)*rGeom[0].Y() + rN(i,1)*rGeom[1].Y() + rN(i,2)*rGeom[2].Y();
            rPositions(i,2) = 0.0;
        }
    }

    void ReseedPositions(const Geometry< Node < 3 > >& rGeom,
                         BoundedMatrix<double, 7, 3 > & rPositions,
                         BoundedMatrix<double, 7, 3 > & rN)
    {
        double one_third = 1.0 / 3.0;
        double one_eight = 0.12; //1.0 / 6.0;
        double three_quarters = 0.76; //2.0 * one_third;

        rN(0, 0) = one_eight;
        rN(0, 1) = one_eight;
        rN(0, 2) = three_quarters;

        rN(1, 0) = three_quarters;
        rN(1, 1) = one_eight;
        rN(1, 2) = one_eight;

        rN(2, 0) = one_eight;
        rN(2, 1) = three_quarters;
        rN(2, 2) = one_eight;

        rN(3, 0) = one_third;
        rN(3, 1) = one_third;
        rN(3, 2) = one_third;

        rN(4, 0) = one_eight;
        rN(4, 1) = 0.44;
        rN(4, 2) = 0.44;

        rN(5, 0) = 0.44;
        rN(5, 1) = one_eight;
        rN(5, 2) = 0.44;

        rN(6, 0) = 0.44;
        rN(6, 1) = 0.44;
        rN(6, 2) = one_eight;


        //first
        rPositions(0, 0) = one_eight * rGeom[0].X() + one_eight * rGeom[1].X() + three_quarters * rGeom[2].X();
        rPositions(0, 1) = one_eight * rGeom[0].Y() + one_eight * rGeom[1].Y() + three_quarters * rGeom[2].Y();
        rPositions(0, 2) = one_eight * rGeom[0].Z() + one_eight * rGeom[1].Z() + three_quarters * rGeom[2].Z();

        //second
        rPositions(1, 0) = three_quarters * rGeom[0].X() + one_eight * rGeom[1].X() + one_eight * rGeom[2].X();
        rPositions(1, 1) = three_quarters * rGeom[0].Y() + one_eight * rGeom[1].Y() + one_eight * rGeom[2].Y();
        rPositions(1, 2) = three_quarters * rGeom[0].Z() + one_eight * rGeom[1].Z() + one_eight * rGeom[2].Z();

        //third
        rPositions(2, 0) = one_eight * rGeom[0].X() + three_quarters * rGeom[1].X() + one_eight * rGeom[2].X();
        rPositions(2, 1) = one_eight * rGeom[0].Y() + three_quarters * rGeom[1].Y() + one_eight * rGeom[2].Y();
        rPositions(2, 2) = one_eight * rGeom[0].Z() + three_quarters * rGeom[1].Z() + one_eight * rGeom[2].Z();

        //fourth
        rPositions(3, 0) = one_third * rGeom[0].X() + one_third * rGeom[1].X() + one_third * rGeom[2].X();
        rPositions(3, 1) = one_third * rGeom[0].Y() + one_third * rGeom[1].Y() + one_third * rGeom[2].Y();
        rPositions(3, 2) = one_third * rGeom[0].Z() + one_third * rGeom[1].Z() + one_third * rGeom[2].Z();

        //fifth
        rPositions(4, 0) = one_eight * rGeom[0].X() + 0.44 * rGeom[1].X() + 0.44 * rGeom[2].X();
        rPositions(4, 1) = one_eight * rGeom[0].Y() + 0.44 * rGeom[1].Y() + 0.44 * rGeom[2].Y();
        rPositions(4, 2) = one_eight * rGeom[0].Z() + 0.44 * rGeom[1].Z() + 0.44 * rGeom[2].Z();

        //sixth
        rPositions(5, 0) = 0.44 * rGeom[0].X() + one_eight * rGeom[1].X() + 0.44 * rGeom[2].X();
        rPositions(5, 1) = 0.44 * rGeom[0].Y() + one_eight * rGeom[1].Y() + 0.44 * rGeom[2].Y();
        rPositions(5, 2) = 0.44 * rGeom[0].Z() + one_eight * rGeom[1].Z() + 0.44 * rGeom[2].Z();

        //seventh
        rPositions(6, 0) = 0.44 * rGeom[0].X() + 0.44 * rGeom[1].X() + one_eight * rGeom[2].X();
        rPositions(6, 1) = 0.44 * rGeom[0].Y() + 0.44 * rGeom[1].Y() + one_eight * rGeom[2].Y();
        rPositions(6, 2) = 0.44 * rGeom[0].Z() + 0.44 * rGeom[1].Z() + one_eight * rGeom[2].Z();
    }

//*****************************************************************
    //3D
    void ReseedPositions(
        const Geometry< Node < 3 > >& rGeom,
        BoundedMatrix<double, 9, 3 > & rPositions,
        BoundedMatrix<double, 9, 4 > & rN)
    {
        double one_quarter = 0.25;
        double small_fraction = 0.1; //1.0 / 6.0;
        double big_fraction = 0.7; //2.0 * one_third;
        double mid_fraction = 0.3; //2.0 * one_third;

        rN(0, 0) = big_fraction;
        rN(0, 1) = small_fraction;
        rN(0, 2) = small_fraction;
        rN(0, 3) = small_fraction;

        rN(1, 0) = small_fraction;
        rN(1, 1) = big_fraction;
        rN(1, 2) = small_fraction;
        rN(1, 3) = small_fraction;

        rN(2, 0) = small_fraction;
        rN(2, 1) = small_fraction;
        rN(2, 2) = big_fraction;
        rN(2, 3) = small_fraction;

        rN(3, 0) = small_fraction;
        rN(3, 1) = small_fraction;
        rN(3, 2) = small_fraction;
        rN(3, 3) = big_fraction;

        rN(4, 0) = one_quarter;
        rN(4, 1) = one_quarter;
        rN(4, 2) = one_quarter;
        rN(4, 3) = one_quarter;

        rN(5, 0) = small_fraction;
        rN(5, 1) = mid_fraction;
        rN(5, 2) = mid_fraction;
        rN(5, 3) = mid_fraction;

        rN(6, 0) = mid_fraction;
        rN(6, 1) = small_fraction;
        rN(6, 2) = mid_fraction;
        rN(6, 3) = mid_fraction;

        rN(7, 0) = mid_fraction;
        rN(7, 1) = mid_fraction;
        rN(7, 2) = small_fraction;
        rN(7, 3) = mid_fraction;

        rN(8, 0) = mid_fraction;
        rN(8, 1) = mid_fraction;
        rN(8, 2) = mid_fraction;
        rN(8, 3) = small_fraction;

        rPositions=ZeroMatrix(9,3);
        for (unsigned int i=0; i!=4; i++) //going through the 4 nodes
        {
            const auto& coordinates = rGeom[i].Coordinates();
            for (unsigned int j=0; j!=9; j++) //going through the 9 particles
            {
                for (unsigned int k=0; k!=3; k++) //x,y,z
                    rPositions(j,k) += rN(j,i) * coordinates[k];
            }
        }
    }



    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ParticleLevelsetUtility& operator=(ParticleLevelsetUtility const& rOther) = delete;

    /// Copy constructor.
    ParticleLevelsetUtility(ParticleLevelsetUtility const& rOther) = delete;


    ///@}

}; // Class ParticleLevelsetUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim>
inline std::istream& operator >> (std::istream& rIStream,
                                  ParticleLevelsetUtility<TDim>& rThis)
                                  {
                                    return rIStream;
                                  }

/// output stream function
template< unsigned int TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ParticleLevelsetUtility<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.
