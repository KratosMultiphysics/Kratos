import time
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.kratos_utilities as kratos_utilities
from chimera_analysis_base_test import ChimeraAnalysisBaseTest

import KratosMultiphysics.KratosUnittest as UnitTest

class FlowOverCylinderMonolithic(ChimeraAnalysisBaseTest):
    def test_MonolithicFlowOverCylinder(self):
        start = time.clock()
        work_folder = "flow_over_cylinder_monolithic"
        settings_file_name = "flow_over_cylinder_monolithic.json"
        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._run_test(settings_file_name)

        end = time.clock()
        print("Time taken for monolithic chimera simulation",end-start)

class FlowOverCylinderFractionalStep(ChimeraAnalysisBaseTest):
    def test_FractionalStepFlowOverCylinder(self):
        start = time.clock()
        work_folder = "flow_over_cylinder_fractionalstep"
        settings_file_name = "flow_over_cylinder_fractionalstep.json"
        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._run_test(settings_file_name)

        end = time.clock()
        print("Time taken for fractional step chimera simulation",end-start)


class FlowOverCrossMonolithic(ChimeraAnalysisBaseTest):
    def test_FlowOverCrossMonolithic(self):
        start = time.clock()
        work_folder = "flow_over_cross_monolithic"
        settings_file_name = "flow_over_cross_monolithic.json"
        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._run_test(settings_file_name)

        end = time.clock()
        print("Time taken for flow over cross using monolithic solver ",end-start)

class FlowOverCrossFractionalStep(ChimeraAnalysisBaseTest):
    def test_FlowOverCrossFractionalStep(self):
        start = time.clock()
        work_folder = "flow_over_cross_fractional_step"
        settings_file_name = "flow_over_cross_fractional_step.json"
        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._run_test(settings_file_name)

        end = time.clock()
        print("Time taken for flow over cross using fractionalstep solver ",end-start)