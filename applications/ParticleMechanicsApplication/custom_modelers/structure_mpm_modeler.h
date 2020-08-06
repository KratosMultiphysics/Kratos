//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//


#if !defined(KRATOS_STRUCTURE_MPM_MODELER_H_INCLUDED )
#define  KRATOS_STRUCTURE_MPM_MODELER_H_INCLUDED


// System includes

// External includes

// Project includes
#include "modeler/modeler.h"

#include "geometries/geometry_data.h"
#include "geometries/coupling_geometry.h"
#include "utilities/quadrature_points_utility.h"
#include "utilities/binbased_fast_point_locator.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Finds realtion between structure FEM and MPM to make Mortar Mapping.
/* Detail class definition.
*/
class KRATOS_API(MAPPING_APPLICATION) StructureMpmModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(StructureMpmModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StructureMpmModeler()
        : Modeler()
    {
        mpModelStructure = nullptr;
        mpModelMpm = nullptr;
    }

    /// Constructor.
    StructureMpmModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
    {
        mpModelStructure = &rModel;
        mpModelMpm = nullptr;
    }

    /// Destructor.
    virtual ~StructureMpmModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<StructureMpmModeler>(rModel, ModelParameters);
    }

    /// Adds the second model part to the modeler.
    void GenerateNodes(ModelPart& ThisModelPart) override
    {
        mpModelMpm = &ThisModelPart.GetModel();
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupGeometryModel() override;

    void UpdateGeometryModel();

    ///@}
    ///@name Model Part generation
    ///@{

    template<
        class TLineGeometriesList,
        class TQuadraturePointGeometriesList = TLineGeometriesList>
    void CreateStructureQuadraturePointGeometries(
        TLineGeometriesList& rInputLineGeometries,
        TQuadraturePointGeometriesList& rOuputQuadraturePointGeometries,
        GeometryData::IntegrationMethod ThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4
        )
    {
        for (auto line : rInputLineGeometries)
        {
            std::vector<GeometryPointerType> qudrature_point_geometries = CreateQuadraturePointsUtility<Node<3>>::Create(line, ThisIntegrationMethod, 2);
            for (IndexType i = 0; i < qudrature_point_geometries.size(); ++i) {
                rOuputQuadraturePointGeometries.push_back(qudrature_point_geometries[i]);
            }
        }
    }

    template<SizeType TDimension, class TQuadraturePointGeometriesList>
    void CreateMpmQuadraturePointGeometries(
        const TQuadraturePointGeometriesList& rInputQuadraturePointGeometries,
        TQuadraturePointGeometriesList& rOuputQuadraturePointGeometries,
        ModelPart& rBackgroundGridModelPart) {
        if (rOuputQuadraturePointGeometries.size() != rInputQuadraturePointGeometries.size()) {
            rOuputQuadraturePointGeometries.resize(rInputQuadraturePointGeometries.size());
        }

        BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
        SearchStructure.UpdateSearchDatabase();
        typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(100);

        // Loop over the submodelpart of rInitialModelPart
        for (IndexType i; i < rInputQuadraturePointGeometries.size(); ++i)
        {
            typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

            array_1d<double, 3> coordinates = rInputQuadraturePointGeometries[i]->Center();

            Element::Pointer p_elem;
            Vector N;

            // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
            bool is_found = SearchStructure.FindPointOnMesh(coordinates, N, p_elem, result_begin);

            if (is_found) {
                rOuputQuadraturePointGeometries.push_back(CreateQuadraturePointsUtility<NodeType>::CreateFromCoordinates(
                    p_elem->pGetGeometry(),
                    coordinates,
                    rInputQuadraturePointGeometries[i]->IntegrationPoints()[0].Weight()));
            }
        }
    }

    template<SizeType TDimension,
        class TConditionsList>
    void UpdateMpmQuadraturePointGeometries(
        const TConditionsList& rInputConditions,
        ModelPart& rBackgroundGridModelPart) {
        BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
        SearchStructure.UpdateSearchDatabase();
        typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(100);

        // Loop over the submodelpart of rInitialModelPart
        for (IndexType i; i < rInputConditions.size(); ++i)
        {
            typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

            array_1d<double, 3> coordinates = rInputConditions[i]->GetGeometry().pGetGeometryPart(0)->Center();

            Element::Pointer p_elem;
            Vector N;

            // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
            bool is_found = SearchStructure.FindPointOnMesh(coordinates, N, p_elem, result_begin);
            array_1d<double, 3> local_coordinates;
            p_elem->GetGeometry().PointLocalCoordinates(local_coordinates, coordinates);

            if (is_found) {
                CreateQuadraturePointsUtility<NodeType>::UpdateFromLocalCoordinates(
                    rInputConditions[i]->GetGeometry().pGetGeometryPart(1),
                    local_coordinates,
                    rInputConditions[i]->pGetGeometry()->IntegrationPoints()[0].Weight(),
                    p_elem->GetGeometry());
            }
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "StructureMpmModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream & rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream & rOStream) const override
    {
    }

    ///@}

private:
    Model* mpModelStructure;

    Model* mpModelMpm;

    void CopySubModelPart(ModelPart& rDestinationMP, ModelPart& rReferenceMP)
    {
        rDestinationMP.SetNodes(rReferenceMP.pNodes());
        ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        rDestinationMP.SetConditions(coupling_conditions.pConditions());
    }

    void CreateInterfaceLineCouplingConditions(ModelPart& rInterfaceModelPart);

    void CheckParameters();

    ///@}
    ///@name Serializer
    ///@{

    friend class Serializer;

    ///@}
}; // Class StructureMpmModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    StructureMpmModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const StructureMpmModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_STRUCTURE_MPM_MODELER_H_INCLUDED  defined
