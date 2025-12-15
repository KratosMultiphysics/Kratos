//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nunez, based on A. Geiser, M. Fusseder, I. Lopez and R. Rossi work
//
#include "compressible_potential_flow_application_variables.h"
#include "incompressible_potential_flow_element.h"
#include "incompressible_perturbation_potential_flow_element.h"
#include "transonic_perturbation_potential_flow_element.h"
#include "compressible_potential_flow_element.h"
#include "embedded_incompressible_potential_flow_element.h"
#include "embedded_compressible_potential_flow_element.h"
#include "adjoint_base_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"

namespace Kratos
{
    template <class TPrimalElement>
    Element::IntegrationMethod AdjointBasePotentialFlowElement<TPrimalElement>::GetIntegrationMethod() const
    {
        return mpPrimalElement->GetIntegrationMethod();
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->Initialize(rCurrentProcessInfo);
        FindUpwindElement(rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        auto adjoint_geometry_data = this->GetGeometry().GetData();
        mpPrimalElement->GetGeometry().SetData(adjoint_geometry_data);
        mpPrimalElement->Set(Flags(*this));
        mpPrimalElement->InitializeSolutionStep(rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {

    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrix,
                                              rRightHandSideVector,
                                              rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       const ProcessInfo& rCurrentProcessInfo)
    {
        MatrixType tmp;
        mpPrimalElement->CalculateLeftHandSide(tmp, rCurrentProcessInfo);
        rLeftHandSideMatrix = trans(tmp);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        const ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVector,
                                                rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateOnIntegrationPoints(const Variable<int>& rVariable,
            std::vector<int>& rValues,
            const ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->CalculateOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::GetValuesVector(Vector& rValues, int Step) const
    {
        KRATOS_TRY
        const AdjointBasePotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);
        const auto& r_geometry = GetGeometry();

        if (wake == 1) // wake element
        {
            if(rValues.size() != 2*NumNodes)
                rValues.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances = PotentialFlowUtilities::GetWakeDistances<Dim, NumNodes>(r_this);
            GetValuesOnSplitElement(rValues,distances);

        }else{ // normal element
            if(rValues.size() != NumNodes)
                rValues.resize(NumNodes, false);

            if(kutta == 0){
                for(unsigned int i=0; i<NumNodes; i++)
                    rValues[i] = r_geometry[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            }else{
                for(unsigned int i=0; i<NumNodes; i++){
                    if (!r_geometry[i].GetValue(TRAILING_EDGE))
                        rValues[i] = r_geometry[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rValues[i] = r_geometry[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling CalculateSensitivityMatrix from adjoint potential flow base element." << std::endl;

        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling CalculateSensitivityMatrix from adjoint potential flow base element." << std::endl;

        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
    {
        const AdjointBasePotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);
        const auto& r_geometry = GetGeometry();

        if(wake == 0)//normal element
        {
            if (rResult.size() != NumNodes+1)
                rResult.resize(NumNodes+1, false);

            if(kutta == 0){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rResult[i] = r_geometry[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (!r_geometry[i].GetValue(TRAILING_EDGE))
                        rResult[i] = r_geometry[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                    else
                        rResult[i] = r_geometry[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL).EquationId();
                }
            }
            if (r_this.IsNot(INLET))
                AddUpwindEquationId(rResult);
        }
        else//wake element
        {
            if (rResult.size() != 2*NumNodes)
                rResult.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances = PotentialFlowUtilities::GetWakeDistances<Dim, NumNodes>(r_this);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rResult[i] = r_geometry[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                else
                    rResult[i] = r_geometry[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0).EquationId();
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rResult[NumNodes+i] = r_geometry[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                else
                    rResult[NumNodes+i] = r_geometry[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0).EquationId();
            }
        }


    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const
    {
        const AdjointBasePotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);
        const auto& r_geometry = GetGeometry();

        if(wake == 0) //normal element
        {
            if (rElementalDofList.size() != NumNodes)
                rElementalDofList.resize(NumNodes);

            if(kutta == 0){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rElementalDofList[i] = r_geometry[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (!r_geometry[i].GetValue(TRAILING_EDGE))
                        rElementalDofList[i] = r_geometry[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rElementalDofList[i] = r_geometry[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        else//wake element
        {
            if (rElementalDofList.size() != 2*NumNodes)
                rElementalDofList.resize(2*NumNodes);

            array_1d<double,NumNodes> distances = PotentialFlowUtilities::GetWakeDistances<Dim, NumNodes>(r_this);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rElementalDofList[i] = r_geometry[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                else
                    rElementalDofList[i] = r_geometry[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rElementalDofList[NumNodes+i] = r_geometry[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                else
                    rElementalDofList[NumNodes+i] = r_geometry[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
            }
        }
    }

    template <class TPrimalElement>
    int AdjointBasePotentialFlowElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
    {

        KRATOS_TRY

        const auto& const_elem_ref = *mpPrimalElement;
        int Check = const_elem_ref.Check(rCurrentProcessInfo);

        if (Check != 0)
        {
            return Check;
        }
        else
        {
            for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
            {
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_VELOCITY_POTENTIAL,
                                                this->GetGeometry()[i]);
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,
                                                    this->GetGeometry()[i]);

                return Check;
            }
        }

        return 0;

        KRATOS_CATCH("");
    }


    /// Turn back information as a string.
    template <class TPrimalElement>
    std::string AdjointBasePotentialFlowElement<TPrimalElement>::Info() const
    {
        std::stringstream buffer;
        buffer << "AdjointBasePotentialFlowElement #" << Id();
        return buffer.str();
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AdjointBasePotentialFlowElement #" << Id();
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::PrintData(std::ostream& rOStream) const
    {
        pGetGeometry()->PrintData(rOStream);
    }

    template <class TPrimalElement>
    Element::Pointer AdjointBasePotentialFlowElement<TPrimalElement>::pGetPrimalElement()
    {
        return mpPrimalElement;
    }

    /*PROTECTED*/

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances ) const
    {

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] > 0)
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            else
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }

        //negative part - sign is opposite to the previous case
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] < 0)
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            else
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }
    }

    /*PRIVATE*/

    template <class TPrimalElement>
    inline GlobalPointer<Element> AdjointBasePotentialFlowElement<TPrimalElement>::pGetUpwindElement() const
    {
        KRATOS_ERROR_IF(mpUpwindElement.get() == nullptr)
            << "No upwind element found for element #" << this->Id() << std::endl;
        return mpUpwindElement;
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::pSetUpwindElement(GlobalPointer<Element> pUpwindElement)
    {
        mpUpwindElement = pUpwindElement;
    }

    template <class TPrimalElement>
    bool AdjointBasePotentialFlowElement<TPrimalElement>::CheckUpwindElement()
    {
        return mpUpwindElement.get() == nullptr;
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::AddUpwindEquationId(
        EquationIdVectorType& rResult) const
    {
        const int additional_upwind_node_index = GetAdditionalUpwindNodeIndex();
        const auto& r_upstream_element = *pGetUpwindElement();
        const auto& r_upwind_geomtery = r_upstream_element.GetGeometry();
        const int upstream_kutta = r_upstream_element.GetValue(KUTTA);
        if (upstream_kutta == 0) { // upwind element is not kutta
            // TODO special treatment for upwind wake elements
            rResult[NumNodes] = r_upwind_geomtery[additional_upwind_node_index].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
        } else { // upwind element is kutta
            if (!r_upwind_geomtery[additional_upwind_node_index].GetValue(TRAILING_EDGE)) {
                // upwind node is not trailing edge
                rResult[NumNodes] = r_upwind_geomtery[additional_upwind_node_index].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
            } else {
                // upwind node is trailing edge
                rResult[NumNodes] = r_upwind_geomtery[additional_upwind_node_index].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL).EquationId();
            }
        }
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::FindUpwindElement(const ProcessInfo& rCurrentProcessInfo)
    {
        GeometryType upwind_element_boundary;
        FindUpwindEdge(upwind_element_boundary, rCurrentProcessInfo);
        std::vector<size_t> upwind_element_nodes;
        PotentialFlowUtilities::GetSortedIds<Dim, NumNodes>(upwind_element_nodes, upwind_element_boundary);

        GlobalPointersVector<Element> upwind_element_candidates;
        PotentialFlowUtilities::GetNodeNeighborElementCandidates<Dim, NumNodes>(upwind_element_candidates, upwind_element_boundary);
        SelectUpwindElement(upwind_element_nodes, upwind_element_candidates);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::FindUpwindEdge(GeometryType& rUpwindEdge,
        const ProcessInfo& rCurrentProcessInfo)
    {

        GeometriesArrayType element_boundary_geometry;
        GetElementGeometryBoundary(element_boundary_geometry);

        // free stream values
        const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];

        double minimum_edge_flux = 0.0;
        for (SizeType i = 0; i < element_boundary_geometry.size(); i++)
        {
            const auto edge_normal = GetEdgeNormal(element_boundary_geometry[i]);

            const double edge_flux = inner_prod(edge_normal, free_stream_velocity);

            if(edge_flux < minimum_edge_flux)
            {
                minimum_edge_flux = edge_flux;
                rUpwindEdge = element_boundary_geometry[i];
            }
        }
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::GetElementGeometryBoundary(GeometriesArrayType& rElementGeometryBoundary)
    {
        const AdjointBasePotentialFlowElement& r_this = *this;

        // current element geometry
        const GeometryType& r_geom = r_this.GetGeometry();

        // get element edges or faces depending on dimension of the problem
        if constexpr (Dim == 2)
        {
            // current element edges
            rElementGeometryBoundary = r_geom.GenerateEdges();
        }
        else if constexpr (Dim == 3)
        {
            // current element faces
            rElementGeometryBoundary = r_geom.GenerateFaces();
        }
    }

    template <class TPrimalElement>
    array_1d<double, 3> AdjointBasePotentialFlowElement<TPrimalElement>::GetEdgeNormal(const GeometryType& rEdge)
    {
        // get local coordinates of edge center
        array_1d<double, 3> edge_center_coordinates;
        rEdge.PointLocalCoordinates(edge_center_coordinates, rEdge.Center());

        // outward pointing normals of each edge
        return rEdge.Normal(edge_center_coordinates);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::SelectUpwindElement(
        std::vector<IndexType>& rUpwindElementNodeIds,
        GlobalPointersVector<Element>& rUpwindElementCandidates)
    {
        for (SizeType i = 0; i < rUpwindElementCandidates.size(); i++)
        {
            // get sorted node ids of neighbording elements
            std::vector<size_t> neighbor_element_ids;
            PotentialFlowUtilities::GetSortedIds<Dim, NumNodes>(neighbor_element_ids, rUpwindElementCandidates[i].GetGeometry());

            // find element which shares the upwind element nodes with current element
            // but is not the current element
            if(std::includes(neighbor_element_ids.begin(), neighbor_element_ids.end(),
                rUpwindElementNodeIds.begin(), rUpwindElementNodeIds.end())
                && rUpwindElementCandidates[i].Id() != this->Id())
            {
                mpUpwindElement = rUpwindElementCandidates(i);
                break;
            }
        }

        // If no upwind element is found, the element is an INLET element and the
        // upwind element pointer points to itself
        if (this->CheckUpwindElement())
        {
            this->pSetUpwindElement(this);
            this->SetFlags(INLET);
        }
    }

    template <class TPrimalElement>
    int AdjointBasePotentialFlowElement<TPrimalElement>::GetAdditionalUpwindNodeIndex() const
    {
        // current and upwind element geometry
        const GeometryType& r_geom = this->GetGeometry();
        const GeometryType& r_upwind_geom = pGetUpwindElement()->GetGeometry();
        std::vector<size_t> element_nodes_ids;
        PotentialFlowUtilities::GetSortedIds<Dim, NumNodes>(element_nodes_ids, r_geom);

        // Search for the Id of the upwind element node that
        // is not contained in the current element
        bool upstream_element_id_found = false;
        // loop over upwind element nodes
        for (unsigned int i = 0; i < NumNodes; i++) {
            if( std::find(element_nodes_ids.begin(), element_nodes_ids.end(),
                r_upwind_geom[i].Id()) == element_nodes_ids.end() )  {
                    upstream_element_id_found = true;
                    return i;
                }
        }

        KRATOS_ERROR_IF(!upstream_element_id_found) << "No upstream element id found for element #"
                << this->Id() << std::endl;

        return -1;
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        rSerializer.save("mpPrimalElement", mpPrimalElement);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
        rSerializer.load("mpPrimalElement", mpPrimalElement);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Template class instantiation

    template class AdjointBasePotentialFlowElement<IncompressiblePotentialFlowElement<2,3>>;
    template class AdjointBasePotentialFlowElement<IncompressiblePerturbationPotentialFlowElement<2,3>>;
    template class AdjointBasePotentialFlowElement<IncompressiblePerturbationPotentialFlowElement<3,4>>;
    template class AdjointBasePotentialFlowElement<CompressiblePotentialFlowElement<2,3>>;
    template class AdjointBasePotentialFlowElement<TransonicPerturbationPotentialFlowElement<2,3>>;
    template class AdjointBasePotentialFlowElement<TransonicPerturbationPotentialFlowElement<3,4>>;
    template class AdjointBasePotentialFlowElement<EmbeddedIncompressiblePotentialFlowElement<2,3>>;
    template class AdjointBasePotentialFlowElement<EmbeddedCompressiblePotentialFlowElement<2,3>>;
} // namespace Kratos.
