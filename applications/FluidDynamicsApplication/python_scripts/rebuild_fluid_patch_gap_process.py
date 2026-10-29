#     |  /           |
#     ' /   __| _` | __|  _ \   __|
#     . \  |   (   | |   (   |\__ `
#    _|\_\_|  \__,_|\__|\___/ ____/
#                    Multi-Physics
#
#   Fluid Dynamics Application
#
#   License:         BSD License
#                    Kratos default license: kratos/license.txt
#
#   Main author:     Nicolò Antonelli
#

"""Rebuild positive-gap entities after prescribed patch motion.

The caller must clear its strategy before Execute and recreate it afterwards.
Nodes, DOFs, bulk elements, nodal histories and materials are retained. History
initialization is deliberately left to the caller (exact fields in the MMS).
"""
import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as CFD


class RebuildFluidPatchGapProcess(KM.Process):
    def __init__(self, model, coupling_settings, geometry_name="Interfaces",
                 computing_name="Computing", tolerance=1e-10):
        super().__init__()
        self.model = model
        self.settings = [settings.Clone() for settings in coupling_settings]
        if not self.settings:
            raise ValueError("Provide all patch pairs sharing the background.")
        self.background = model.GetModelPart(self.settings[0]["background_model_part_name"].GetString())
        self.root = self.background.GetRootModelPart()
        if any(s["background_model_part_name"].GetString() != self.background.FullName()
               for s in self.settings):
            raise ValueError("All pairs must share one background.")
        self.inner = [model.GetModelPart(s["inner_model_part_name"].GetString()) for s in self.settings]
        self.geometry_name = geometry_name
        self.computing_name = computing_name
        self.tolerance = tolerance

    def GetPatchGeometry(self, index):
        """Return the geometry view corresponding to an inner patch index."""
        if index < 0 or index >= len(self.inner):
            raise IndexError('Inner patch index is out of range.')
        return self.root.GetSubModelPart(self.geometry_name).GetSubModelPart(f'Patch_{index + 1}')

    def Execute(self):
        # Delete only entities owned by the previous gap construction. Bulk
        # inactive elements and their nodes are never removed from the root.
        if self.root.HasSubModelPart(self.computing_name):
            self.root.RemoveSubModelPart(self.computing_name)
        for settings in self.settings:
            name = settings["interface_model_part_name"].GetString()
            if self.model.HasModelPart(name):
                part = self.model.GetModelPart(name)
                for entity_id in [e.Id for e in part.Elements]:
                    self.root.RemoveElementFromAllLevels(entity_id)
                for entity_id in [c.Id for c in part.Conditions]:
                    self.root.RemoveConditionFromAllLevels(entity_id)
                part.GetParentModelPart().RemoveSubModelPart(part.Name)
        if self.root.HasSubModelPart(self.geometry_name):
            part = self.root.GetSubModelPart(self.geometry_name)
            for geometry_id in [g.Id for g in part.Geometries]:
                self.root.RemoveGeometryFromAllLevels(geometry_id)
            self.root.RemoveSubModelPart(self.geometry_name)

        owners = KM.FemMultipatchGapUtilities.UpdateActivation(self.background, self.inner, self.tolerance)
        views = self.root.CreateSubModelPart(self.geometry_name)
        for index, (inner, settings) in enumerate(zip(self.inner, self.settings), 1):
            # Refresh the modeler's retained per-hole membership as well.
            name = f"Inactive_Patch_{index}"
            if self.background.HasSubModelPart(name):
                self.background.RemoveSubModelPart(name)
            inactive = self.background.CreateSubModelPart(name)
            inactive.AddElements([entity_id for entity_id, owner in owners.items() if owner == index])
            inactive.AddNodes(sorted({n.Id for e in inactive.Elements for n in e.GetGeometry()}))
            KM.FemMultipatchGapUtilities.CreateGeometry(
                self.background, inner, views.CreateSubModelPart(f"Patch_{index}"), self.tolerance)
            if settings.Has("allow_transient"):
                settings["allow_transient"].SetBool(True)
            else:
                settings.AddBool("allow_transient", True)
            CFD.CreateFluidPatchGapCouplingConditionsProcess(self.model, settings).Execute()

        self.computing = self.root.CreateSubModelPart(self.computing_name)
        self.computing.AddElements([e.Id for patch in [self.background, *self.inner]
                                   for e in patch.Elements if e.Is(KM.ACTIVE)])
        for settings in self.settings:
            part = self.model.GetModelPart(settings["interface_model_part_name"].GetString())
            self.computing.AddElements([e.Id for e in part.Elements])
            self.computing.AddConditions([c.Id for c in part.Conditions])
        for patch in self.inner:
            if patch.HasSubModelPart('CylinderBoundary'):
                self.computing.AddConditions([c.Id for c in patch.GetSubModelPart('CylinderBoundary').Conditions])
        if self.root.HasSubModelPart('CylinderGapVolumes'):
            self.computing.AddElements([e.Id for e in self.root.GetSubModelPart('CylinderGapVolumes').Elements])
            self.computing.AddConditions([c.Id for c in self.root.GetSubModelPart('CylinderGapInterfaces').Conditions])
        ids = {n.Id for e in self.computing.Elements for n in e.GetGeometry()}
        ids.update(n.Id for c in self.computing.Conditions for n in c.GetGeometry())
        self.computing.AddNodes(sorted(ids))
