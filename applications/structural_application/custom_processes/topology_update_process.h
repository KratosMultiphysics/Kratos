/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 31 Mar 2016 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_TOPOLOGY_UPDATE_PROCESS_H_INCLUDED )
#define  KRATOS_TOPOLOGY_UPDATE_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
//#include "spatial_containers/bins_static.h"
#include "spatial_containers/bins_dynamic.h"
#include "structural_application.h"


//#define DEBUG_WITH_MPI


#ifdef DEBUG_WITH_MPI
#include "mpi.h"
#endif


namespace Kratos
{

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

/// Short class definition.
/** Detail class definition.
*/
class TopologyUpdateProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    typedef Element::GeometryType GeometryType;
    typedef GeometryType::PointType NodeType;
    typedef GeometryType::PointType::PointType PointType;
    typedef ModelPart::ElementsContainerType ElementsContainerType;

    /// Pointer definition of TopologyUpdateProcess
    KRATOS_CLASS_POINTER_DEFINITION(TopologyUpdateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// This constructor will take all elements of the model_part for topology optimization
    TopologyUpdateProcess(ModelPart& r_model_part, double volfrac, double rmin, int ft)
    : mr_model_part(r_model_part), mvolfrac(volfrac), mrmin(rmin), mtotal_volume(0.0), mft(ft), mbinsize(100)
    {
        mrElements = r_model_part.Elements();
    }

    /// Destructor.
    virtual ~TopologyUpdateProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void SetBinSize(std::size_t size)
    {
        mbinsize = size;
    }

    std::size_t GetBinSize() const
    {
        return mbinsize;
    }

    virtual void Execute()
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {
        std::cout << "At TopologyUpdateProcess::" << __FUNCTION__ << std::endl;
        double start = OpenMPUtils::GetCurrentTime();
//        KRATOS_WATCH(mvolfrac)
//        KRATOS_WATCH(mrmin)
//        KRATOS_WATCH(mft)

        /* find the neighbour elements and initilize filter structure */
        // firstly compute the center of all elements and store it to list
        typedef std::vector<NodeType::Pointer> NodesContainerType;
        NodesContainerType Centers;
        for(ElementsContainerType::iterator i_element = mrElements.begin(); i_element != mrElements.end(); ++i_element)
        {
            NodeType::Pointer p_new_node = NodeType::Pointer(new NodeType(i_element->Id(), i_element->GetGeometry().Center()));
            Centers.push_back(p_new_node);
        }

        // create the spatial bin containing all the centers
        typedef BinsDynamic<3, PointType, NodesContainerType> BinsType;
        BinsType SpatialBin(Centers.begin(), Centers.end());

        // loop through each center and find out the one in range
        // TODO parallize the neighbour search (look for bins_dynamic.h, line 482)
        std::size_t max_number_of_results = mbinsize;
        NodesContainerType Results(max_number_of_results);
        std::vector<double> Distances(max_number_of_results);
        for(NodesContainerType::iterator it = Centers.begin(); it != Centers.end(); ++it)
        {
            NodeType::Pointer& ThisCenter = *it;

            std::size_t num_results = SpatialBin.SearchInRadius(*ThisCenter, mrmin, Results.begin(), Distances.begin(), max_number_of_results); // TODO parameterize this

            if(num_results >= max_number_of_results)
            {
                KRATOS_WATCH((*it)->Id())
                KRATOS_THROW_ERROR(std::runtime_error, "The number of founded neighbour elements is larger than maximum number of results.", "Try to increase the maximum number of results.")
            }

            if(num_results > 0)
            {
//                mrElements[(*it)->Id()].GetValue(NEIGHBOUR_ELEMENTS).reserve(num_results); // TODO do we need this?
                for(std::size_t i = 0; i < num_results; ++i)
                {
                    Element::WeakPointer temp = mr_model_part.pGetElement(Results[i]->Id());
                    AddUniqueWeakPointer(mrElements[(*it)->Id()].GetValue(NEIGHBOUR_ELEMENTS), temp);
                }

                Vector weights(num_results);
                for(std::size_t i = 0; i < num_results; ++i)
                    weights(i) = mrmin - sqrt(Distances[i]);
                mrElements[(*it)->Id()].GetValue(NEIGHBOUR_WEIGHTS) = weights;
//                KRATOS_WATCH(weights)
            }
        }

        /* Initialize the material density and young modulus for each element */
        for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
        {
            i_element->SetValue(MATERIAL_DENSITY, mvolfrac);
            i_element->SetValue(MATERIAL_DENSITY_FILTERED, mvolfrac);

            const GeometryType::IntegrationPointsArrayType& integration_points = i_element->GetGeometry().IntegrationPoints( i_element->GetIntegrationMethod() );
            double young = GetModulus(i_element->GetProperties(), i_element->GetValue(MATERIAL_DENSITY_FILTERED));
            std::vector<double> Modulus(integration_points.size(), young);
            i_element->SetValueOnIntegrationPoints(YOUNG_MODULUS, Modulus, mr_model_part.GetProcessInfo());
        }

        /* Compute the total volume of the system */
        mtotal_volume = 0.0;
        for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
            mtotal_volume += i_element->GetValue(GEOMETRICAL_DOMAIN_SIZE);
        mtotal_volume = this->SumAll(mtotal_volume);

        double stop = OpenMPUtils::GetCurrentTime();
        std::cout << "TopologyUpdateProcess::" << __FUNCTION__ << " completed: " << (stop - start) << std::endl;
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep()
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
        std::cout << "At TopologyUpdateProcess::" << __FUNCTION__ << std::endl;

        #ifdef DEBUG_WITH_MPI
        int mpi_rank, mpi_size;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        #endif

        /*** OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS ***/
        double Compliance = 0.0;
        for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
        {
            // get the integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = i_element->GetGeometry().IntegrationPoints( i_element->GetIntegrationMethod() );

            // compute strain energy (i.e compliance) at the integration points
            std::vector<double> StrainEnergyAtIntegrationPoints;
            i_element->GetValueOnIntegrationPoints(STRAIN_ENERGY, StrainEnergyAtIntegrationPoints, mr_model_part.GetProcessInfo());

            // get the Jacobian at integration points
            std::vector<double> Jacobian;
            i_element->GetValueOnIntegrationPoints(JACOBIAN_0, Jacobian, mr_model_part.GetProcessInfo());

            // compute strain energy at the element
            double StrainEnergy = 0.0;
            for(unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
            {
                StrainEnergy += integration_points[PointNumber].Weight() * StrainEnergyAtIntegrationPoints[PointNumber] * Jacobian[PointNumber];
            }

            // add to system compliance
//            double aux = GetModulus(i_element->GetProperties(), i_element->GetValue(MATERIAL_DENSITY_FILTERED));
//            KRATOS_WATCH(aux)
//            KRATOS_WATCH(2.0 * StrainEnergy / aux)
//            KRATOS_WATCH(2.0 * StrainEnergy)
            Compliance += 2.0 * StrainEnergy;

            // compute the derivative of the strain energy w.r.t physical material density
            // TODO check why -1
            double E = GetModulus(i_element->GetProperties(), i_element->GetValue(MATERIAL_DENSITY_FILTERED));
            double dE = GetDerivativeModulus(i_element->GetProperties(), i_element->GetValue(MATERIAL_DENSITY_FILTERED));
            i_element->SetValue(ELEMENT_DC, (dE / E) * 2.0 * StrainEnergy);

            // compute the derivative of the volume w.r.t physical material density
            i_element->SetValue(ELEMENT_DV, 1.0);
        }
        mobjective = Compliance;
//        KRATOS_WATCH(Compliance)

        #ifdef DEBUG_WITH_MPI
        if(mpi_rank == 0)
        {
            std::cout << "BaseType::mrElements (proc 0 before synchronize ELEMENT_DC and ELEMENT_DV):";
            for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
                std::cout << " " << i_element->Id();
            std::cout << std::endl;
        }
        #endif

        // synchronize
        this->Synchronize(ELEMENT_DC);
        this->Synchronize(ELEMENT_DV);

        #ifdef DEBUG_WITH_MPI
        if(mpi_rank == 0)
        {
            std::cout << "BaseType::mrElements (proc 0 after synchronize ELEMENT_DC and ELEMENT_DV):";
            for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
                std::cout << " " << i_element->Id();
            std::cout << std::endl;
        }
        #endif

        /*** FILTERING/MODIFICATION OF SENSITIVITIES ***/
        if(mft == 1)
        {
            double tol = 1.0e-3; // TODO parameterize this
            for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
            {
                double coeff = std::max(tol, i_element->GetValue(MATERIAL_DENSITY));
		        double nom = 0.0, denom = 0.0;
		        WeakPointerVector<Element>& rE = i_element->GetValue(NEIGHBOUR_ELEMENTS);
                const Vector& weights = i_element->GetValue(NEIGHBOUR_WEIGHTS);
		        for(std::size_t i = 0; i < rE.size(); ++i)
                {
                    nom += weights[i] * rE[i].GetValue(MATERIAL_DENSITY) * rE[i].GetValue(ELEMENT_DC);
                    denom += weights[i];
                }
                denom *= coeff;
                #ifdef DEBUG_WITH_MPI
                if(mpi_rank == 0)
                    std::cout << i_element->Id() << ": nom: " << nom << ", denom: " << denom << ", coeff: " << coeff << std::endl;
                #endif
                i_element->SetValue(ELEMENT_DC_FILTERED, nom/denom);

                i_element->SetValue(ELEMENT_DV_FILTERED, i_element->GetValue(ELEMENT_DV));
            }
        }
        else if(mft == 2)
        {
            for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
            {
                double nom = 0.0, denom = 0.0;
                WeakPointerVector<Element>& rE = i_element->GetValue(NEIGHBOUR_ELEMENTS);
                const Vector& weights = i_element->GetValue(NEIGHBOUR_WEIGHTS);
                for(std::size_t i = 0; i < rE.size(); ++i)
                {
                    nom += weights[i] * rE[i].GetValue(ELEMENT_DC);
                    denom += weights[i];
                }
                i_element->SetValue(ELEMENT_DC_FILTERED, nom/denom);

                nom = 0.0; denom = 0.0;
                for(std::size_t i = 0; i < rE.size(); ++i)
                {
                    nom += weights[i] * rE[i].GetValue(ELEMENT_DV);
                    denom += weights[i];
                }
                i_element->SetValue(ELEMENT_DV_FILTERED, nom/denom);
            }
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid type of filter:", mft)

        // synchronize
        this->Synchronize(ELEMENT_DC_FILTERED);
        this->Synchronize(ELEMENT_DV_FILTERED);

        /*** OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES ***/
        double l1 = 0.0;
        double l2 = 1.0e9;
        double lmid;
        double move = 0.2; // TODO parameterize these options
        double tol = 1.0e-3;
//        unsigned int cnt = 0;

        while( (l2 - l1) / (l1 + l2) > tol )
        {
            lmid = 0.5 * (l1 + l2);

            // for each element, compute xnew
            for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
            {
                double x = i_element->GetValue(MATERIAL_DENSITY);
                double dc = i_element->GetValue(ELEMENT_DC_FILTERED);
                double dv = i_element->GetValue(ELEMENT_DV_FILTERED);
                double xnew = std::max(0.0, std::max(x-move, std::min(1.0, std::min(x+move, x*sqrt((dc/dv)/lmid)))));
                i_element->SetValue(MATERIAL_DENSITY_NEW, xnew);
            }

            // synchronize
            this->Synchronize(MATERIAL_DENSITY_NEW);

            // for each element, filter the density
            if(mft == 1)
            {
                for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
                    i_element->SetValue(MATERIAL_DENSITY_FILTERED, i_element->GetValue(MATERIAL_DENSITY_NEW));
            }
            else if(mft == 2)
            {
                for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
                {
                    WeakPointerVector<Element>& rE = i_element->GetValue(NEIGHBOUR_ELEMENTS);
                    const Vector& weights = i_element->GetValue(NEIGHBOUR_WEIGHTS);
                    double nom = 0.0, denom = 0.0;
		            for(std::size_t i = 0; i < rE.size(); ++i)
                    {
                        nom += weights[i] * rE[i].GetValue(MATERIAL_DENSITY_NEW);
                        denom += weights[i];
                    }
                    i_element->SetValue(MATERIAL_DENSITY_FILTERED, nom/denom);
                }
            }

            // synchronize
            this->Synchronize(MATERIAL_DENSITY_FILTERED);

            // update the bisection
            double total_density_volume = 0.0;
            for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
            {
                total_density_volume += i_element->GetValue(MATERIAL_DENSITY_FILTERED) * i_element->GetValue(GEOMETRICAL_DOMAIN_SIZE);
//                total_density_volume += i_element->GetValue(MATERIAL_DENSITY_FILTERED);
//                KRATOS_WATCH(i_element->GetValue(GEOMETRICAL_DOMAIN_SIZE))
            }
            total_density_volume = this->SumAll(total_density_volume);
            if(total_density_volume > mvolfrac * mtotal_volume)
            {
                l1 = lmid;
            }
            else
            {
                l2 = lmid;
            }
//            ++cnt;
//            KRATOS_WATCH(cnt)
//            KRATOS_WATCH(total_density_volume)
//            KRATOS_WATCH(mtotal_volume)
//            KRATOS_WATCH(mvolfrac)
//            KRATOS_WATCH(l1)
//            KRATOS_WATCH(l2)
//            KRATOS_WATCH(lmid)
        }

        // compute the change in total density
        mchange = 0.0;
        for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
        {
            double aux = fabs(i_element->GetValue(MATERIAL_DENSITY_NEW) - i_element->GetValue(MATERIAL_DENSITY));
            if(mchange < aux)
                mchange = aux;
        }
        mchange = this->MaxAll(mchange);

        // update the material density
        for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
        {
            i_element->SetValue(MATERIAL_DENSITY, i_element->GetValue(MATERIAL_DENSITY_NEW));
        }

        // synchronize
        this->Synchronize(MATERIAL_DENSITY);

        // update the element YOUNG_MODULUS
        for(ElementsContainerType::iterator i_element = mrElements.begin() ; i_element != mrElements.end(); ++i_element)
        {
            const GeometryType::IntegrationPointsArrayType& integration_points = i_element->GetGeometry().IntegrationPoints( i_element->GetIntegrationMethod() );
            double young = GetModulus(i_element->GetProperties(), i_element->GetValue(MATERIAL_DENSITY_FILTERED));
            std::vector<double> Modulus(integration_points.size(), young);
            i_element->SetValueOnIntegrationPoints(YOUNG_MODULUS, Modulus, mr_model_part.GetProcessInfo());
        }

        std::cout << "TopologyUpdateProcess::" << __FUNCTION__ << " completed" << std::endl;
    }

    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
    }

    void ClearNeighbours()
    {
        for(ElementsContainerType::iterator i_element = mrElements.begin(); i_element != mrElements.end(); ++i_element)
        {
            WeakPointerVector<Element>& rE = i_element->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(), rE.end());
        }
    }

    double GetTopologyChange() const {return mchange;}
    double GetObjective() const {return mobjective;}

    /// Reduce-Sum over all processes
    virtual double SumAll(double value) const
    {
        return value;
    }

    /// Reduce-Max over all processes
    virtual double MaxAll(double value) const
    {
        return value;
    }

    /// Synchronize the elemental values across processes
    virtual void Synchronize(const Variable<int>& rVariable)
    {
        // DO NOTHING
    }

    virtual void Synchronize(const Variable<double>& rVariable)
    {
        // DO NOTHING
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

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "TopologyUpdateProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TopologyUpdateProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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

    ModelPart& mr_model_part;
    ElementsContainerType mrElements;
    double mvolfrac;        // volume fraction
    double mrmin;           // radius for neighbour search
    double mtotal_volume;   // total volume of the domain
    int mft;                // filter type: 1: sensitivity filter, 2: density filter
    std::size_t mbinsize;   // for neighbour search using Bin

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    //******************************************************************************************
    //******************************************************************************************
    template<class TDataType>
    void AddUniqueWeakPointer(WeakPointerVector<TDataType>& v, const typename TDataType::WeakPointer candidate)
    {
        typename WeakPointerVector< TDataType >::iterator i = v.begin();
        typename WeakPointerVector< TDataType >::iterator endit = v.end();
        while ( i != endit && (i)->Id() != (candidate.lock())->Id())
        {
            i++;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }

    }

    double GetModulus(const Properties& props, double xe) const
    {
        double E0 = props[YOUNG_MODULUS];
        double Emin = props[YOUNG_MODULUS_MIN];
        double penal = props[PENALIZATION_FACTOR];
        return Emin + pow(xe, penal) * (E0 - Emin);
    }

    double GetDerivativeModulus(const Properties& props, double xe) const
    {
        double E0 = props[YOUNG_MODULUS];
        double Emin = props[YOUNG_MODULUS_MIN];
        double penal = props[PENALIZATION_FACTOR];
        return penal * pow(xe, penal - 1) * (E0 - Emin);
    }


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

    double mchange;     // convergence criteria
    double mobjective;  // value of objective function

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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
    TopologyUpdateProcess& operator=(TopologyUpdateProcess const& rOther);

    /// Copy constructor.
    //TopologyUpdateProcess(FindConditionsNeighboursProcess const& rOther);


    ///@}

}; // Class FindConditionsNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TopologyUpdateProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TopologyUpdateProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#undef DEBUG_WITH_MPI

#endif // KRATOS_TOPOLOGY_UPDATE_PROCESS_H_INCLUDED  defined
