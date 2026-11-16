"""Workflow fixtures for the lightweight Dam initial-state tests.

The helpers build a small thermo-mechanical DamAnalysis and execute it through
a child process so the production solver and process callbacks are used.
"""
import json
from pathlib import Path
import subprocess
import tempfile

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as Dam
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_dam_initial_stress_handling as mechanical_tests
import test_dam_thermo_mechanic_initial_output as thermal_tests

_MDPA = """Begin ModelPartData
End ModelPartData
Begin Properties 0
End Properties
Begin Properties 1
DENSITY 2.0
YOUNG_MODULUS 10000.0
POISSON_RATIO 0.2
THICKNESS 1.0
CONSTITUTIVE_LAW_NAME LinearElasticPlaneStrain2DLaw
End Properties
Begin Nodes
1 0 0 0
2 1 0 0
3 0 1 0
4 1 1 0
5 0 2 0
6 1 2 0
End Nodes
Begin Elements SmallDisplacementElement2D4N
1 1 1 2 4 3
2 1 3 4 6 5
End Elements
Begin NodalData DISPLACEMENT_X
1 1 0.0
2 1 0.0
End NodalData
Begin NodalData DISPLACEMENT_Y
1 1 0.0
2 1 0.0
End NodalData
Begin SubModelPart Parts_Concrete
Begin SubModelPartNodes
1
2
3
4
5
6
End SubModelPartNodes
Begin SubModelPartElements
1
2
End SubModelPartElements
End SubModelPart
"""


def snapshot(part):
    residual = [0.] * 12
    external = [0.] * 12
    stress, strain, thermal_strain = [], [], []
    for element in part.Elements:
        stress.extend([list(v) for v in element.CalculateOnIntegrationPoints(
            KM.CAUCHY_STRESS_VECTOR, part.ProcessInfo)])
        strain.extend([list(v) for v in element.CalculateOnIntegrationPoints(
            KM.GREEN_LAGRANGE_STRAIN_VECTOR, part.ProcessInfo)])
        if part.HasNodalSolutionStepVariable(Dam.NODAL_REFERENCE_TEMPERATURE):
            thermal_strain.extend([list(v) for v in element.CalculateOnIntegrationPoints(
                KM.KratosGlobals.GetVariable("THERMAL_STRAIN_VECTOR"), part.ProcessInfo)])
        lhs, rhs = KM.Matrix(), KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, part.ProcessInfo)
        # Two unit-square Q4s with constant density and body acceleration.
        weight = element.Properties[KM.DENSITY] / 4.
        for i, node in enumerate(element.GetGeometry()):
            for j in range(2):
                index = 2 * (node.Id - 1) + j
                residual[index] += rhs[2*i+j]
                external[index] += weight * node.GetSolutionStepValue(KM.VOLUME_ACCELERATION)[j]
    result = {
        'time': part.ProcessInfo[KM.TIME], 'step': part.ProcessInfo[KM.STEP],
        'restarted': part.ProcessInfo[KM.IS_RESTARTED],
        'stress': stress, 'strain': strain, 'thermal_strain': thermal_strain,
        'residual': residual, 'external': external,
        'internal': [e-r for e, r in zip(external, residual)],
        'coordinates': [[n.X, n.Y, n.Z] for n in part.Nodes],
        'fixity': [[n.IsFixed(KM.DISPLACEMENT_X), n.IsFixed(KM.DISPLACEMENT_Y)] for n in part.Nodes],
    }
    for name, variable in (('displacement', KM.DISPLACEMENT), ('reaction', KM.REACTION),
                           ('velocity', KM.VELOCITY), ('acceleration', KM.ACCELERATION),
                           ('body_force', KM.VOLUME_ACCELERATION)):
        if part.HasNodalSolutionStepVariable(variable):
            result[name] = [list(n.GetSolutionStepValue(variable)) for n in part.Nodes]
        else:
            result[name] = [[0., 0., 0.]] * len(part.Nodes)
    result['previous_displacement'] = [list(n.GetSolutionStepValue(KM.DISPLACEMENT, 1)) for n in part.Nodes]
    if part.HasNodalSolutionStepVariable(Dam.NODAL_REFERENCE_TEMPERATURE):
        result['reference_temperature'] = [n.GetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE) for n in part.Nodes]
        result['temperature'] = [n.GetSolutionStepValue(KM.TEMPERATURE) for n in part.Nodes]
    return result


class TestDamInitialStateWorkflow(KratosUnittest.TestCase):
    @staticmethod
    def _load(value, start=0.):
        return {'python_module': 'assign_vector_variable_process',
                'kratos_module': 'KratosMultiphysics',
                'Parameters': {'model_part_name': 'MainModelPart.Parts_Concrete',
                               'variable_name': 'VOLUME_ACCELERATION', 'constrained': False,
                               'interval': [start, 'End'], 'value': value}}

    def _parameters(self, thermal=False, increment='horizontal',
                    temperature_offset=0., reference_shift=0.):
        settings = mechanical_tests.TestDamInitialStressHandling._base_parameters('result', False, 3.)
        settings['solver_settings']['processes_sub_model_part_list'] = []
        settings['constraints_process_list'] = []
        settings['transfer_results_process']['add_previous_results'] = False
        settings['loads_process_list'] = [self._load([0., -10., 0.])]
        if increment == 'horizontal':
            settings['loads_process_list'].append(self._load([2., -10., 0.], 2.5))
        if thermal:
            thermal_settings = thermal_tests.TestDamThermoMechanicInitialOutput()._base_parameters('unused')
            settings['solver_settings']['solver_type'] = 'dam_thermo_mechanic_solver'
            thermal_solver = thermal_settings['solver_settings']['thermal_solver_settings']
            thermal_solver['move_mesh_flag'] = False
            settings['solver_settings']['thermal_solver_settings'] = thermal_solver
            settings['constraints_process_list'] = [{
                'python_module': 'impose_nodal_reference_temperature_process',
                'kratos_module': 'KratosMultiphysics.DamApplication',
                'Parameters': {'model_part_name': 'MainModelPart',
                               'variable_name': 'NODAL_REFERENCE_TEMPERATURE',
                               'initial_value': -999., 'input_file_name': 'reference_temperature.txt',
                               'interval': [0., 0.]}}]
            expression = '{}+2*y+{}'.format(20+reference_shift, temperature_offset)
            if increment == 'thermal':
                expression = expression + '+5 if t>2.5 else ' + expression
            # Top/bottom thermal Dirichlet data; the middle row is solved.
            settings['constraints_process_list'].append({
                'python_module': 'assign_scalar_variable_process', 'kratos_module': 'KratosMultiphysics',
                'Parameters': {'model_part_name': 'MainModelPart.Thermal_Boundary',
                               'variable_name': 'TEMPERATURE', 'constrained': True,
                               'value': expression, 'interval': [0., 'End']}})
            settings['output_configuration']['result_file_configuration']['nodal_results'] = [
                'DISPLACEMENT', 'TEMPERATURE', 'NODAL_REFERENCE_TEMPERATURE']
        return settings

    def _run(self, directory, settings, thermal=False, offset=0., fully_fixed=False,
             reference_shift=0., expected_error=None):
        directory.mkdir(parents=True, exist_ok=True)
        mdpa = _MDPA
        if thermal:
            mdpa = mdpa.replace('LinearElasticPlaneStrain2DLaw', 'ThermalLinearElastic2DPlaneStrain')
            mdpa = mdpa.replace('THICKNESS 1.0', 'THICKNESS 1.0\nTHERMAL_EXPANSION 0.00001\nCONDUCTIVITY 1.0\nSPECIFIC_HEAT 1.0')
            mdpa += '\nBegin Elements EulerianConvDiff2D4N\n3 1 1 2 4 3\n4 1 3 4 6 5\nEnd Elements\n'
            mdpa += '\nBegin SubModelPart Thermal_Concrete\nBegin SubModelPartNodes\n1\n2\n3\n4\n5\n6\nEnd SubModelPartNodes\nBegin SubModelPartElements\n3\n4\nEnd SubModelPartElements\nEnd SubModelPart\n'
            mdpa += '\nBegin SubModelPart Thermal_Boundary\nBegin SubModelPartNodes\n1\n2\n5\n6\nEnd SubModelPartNodes\nEnd SubModelPart\n'
            mdpa += '\nBegin NodalData TEMPERATURE\n'
            mdpa += ''.join('{} 0 {}\n'.format(i, 20+reference_shift+2*((i-1)//2)+offset) for i in range(1, 7))
            mdpa += 'End NodalData\n'
            for variable, value in (('DENSITY', 2.), ('SPECIFIC_HEAT', 1.), ('CONDUCTIVITY', 1.)):
                mdpa += '\nBegin NodalData '+variable+'\n'
                mdpa += ''.join('{} 0 {}\n'.format(i, value) for i in range(1, 7))
                mdpa += 'End NodalData\n'
        if fully_fixed:
            for variable in ('DISPLACEMENT_X', 'DISPLACEMENT_Y'):
                mdpa = mdpa.replace('Begin NodalData '+variable+'\n1 1 0.0\n2 1 0.0',
                                    'Begin NodalData '+variable+'\n'+''.join('{} 1 0.0\n'.format(i) for i in range(1, 7)).rstrip())
        (directory / 'mini.mdpa').write_text(mdpa)
        (directory / 'reference_temperature.txt').write_text(''.join(
            '{} {}\n'.format(i, 20+reference_shift+2*((i-1)//2)) for i in range(1, 7)))
        (directory / 'parameters.json').write_text(json.dumps(settings))
        tests = Path(__file__).resolve().parent
        import sys
        result = subprocess.run([sys.executable, str(tests / 'selective_initial_stress/run.py'),
                                 str(tests.parent / 'python_scripts')], cwd=directory,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
        (directory / 'run.log').write_text(result.stdout)
        if expected_error is not None:
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn(expected_error, result.stdout)
            return None
        self.assertEqual(result.returncode, 0, result.stdout)
        return json.loads((directory / 'snapshots.json').read_text())

    def _compare(self, actual, expected, atol=1.e-10):
        if isinstance(expected, list):
            self.assertEqual(len(actual), len(expected))
            return max((self._compare(a, e, atol) for a, e in zip(actual, expected)), default=0.)
        error = abs(actual-expected)
        self.assertLessEqual(error, atol + 1.e-10 * abs(expected))
        return error

    def test_existing_reference_temperature_process(self):
        with tempfile.TemporaryDirectory(prefix='dam_reference_temperature_') as tmp:
            settings = self._parameters(thermal=True, increment='none', temperature_offset=10.)
            settings['loads_process_list'] = [self._load([0., 0., 0.])]
            result = self._run(Path(tmp), settings, thermal=True, offset=10., fully_fixed=True)
            # Plane strain: alpha*(1+nu)*dT; C*[e,e,0] = [E/(1+nu)/(1-2nu)*e,...].
            for state in [result['initial']] + result['steps']:
                self._compare(state['reference_temperature'], [20., 20., 22., 22., 24., 24.])
                self._compare(state['temperature'], [30., 30., 32., 32., 34., 34.])
                self._compare(state['thermal_strain'], [[0.00012, 0.00012, 0.]] * 8)
                self._compare(state['stress'], [[-5./3., -5./3., 0.]] * 8)
                self._compare(state['displacement'], [[0., 0., 0.]] * 6)
            self._compare(result['initial']['reaction'], [[0., 0., 0.]] * 6)
            self.assertGreater(max(abs(v) for row in result['steps'][0]['reaction'] for v in row), 0.1)

    def test_use_process_based_local_equivalent_strain_is_registered(self):
        # Dam writes this variable to ProcessInfo. Without the C++ registration the
        # core serializer cannot resolve it when loading a real Dam checkpoint.
        variable = KM.KratosGlobals.GetVariable("USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN")
        self.assertEqual(variable.Name(), "USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN")
        with tempfile.TemporaryDirectory(prefix='dam_variable_registration_') as tmp:
            part = KM.Model().CreateModelPart("MainModelPart")
            part.ProcessInfo[Dam.USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN] = True
            path = str(Path(tmp) / "checkpoint")
            KM.FileSerializer(path, KM.SerializerTraceType.SERIALIZER_NO_TRACE).Save(
                "MainModelPart", part)
            loaded = KM.Model().CreateModelPart("MainModelPart")
            KM.FileSerializer(path, KM.SerializerTraceType.SERIALIZER_NO_TRACE).Load(
                "MainModelPart", loaded)
            self.assertTrue(
                loaded.ProcessInfo[Dam.USE_PROCESS_BASED_LOCAL_EQUIVALENT_STRAIN])


if __name__ == '__main__':
    KratosUnittest.main()
