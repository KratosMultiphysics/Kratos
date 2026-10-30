"""Selective stress initialization: real DamAnalysis, independent source Model.

The stress tests use the production DamApplication utility and process. The
source checkpoint is loaded into an isolated Model; the destination is a fresh
DamAnalysis that keeps its own geometry, time, step, loads and constraints.
"""
import json
import os
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


class TestDamSelectiveInitialization(KratosUnittest.TestCase):
    @staticmethod
    def _load(value, start=0.):
        return {'python_module': 'assign_vector_variable_process',
                'kratos_module': 'KratosMultiphysics',
                'Parameters': {'model_part_name': 'MainModelPart.Parts_Concrete',
                               'variable_name': 'VOLUME_ACCELERATION', 'constrained': False,
                               'interval': [start, 'End'], 'value': value}}

    def _parameters(self, thermal=False, increment='horizontal', checkpoint=False, source=None,
                    temperature_offset=0., reference_shift=0.):
        settings = mechanical_tests.TestDamInitialStressHandling._base_parameters('result', False, 3.)
        settings['solver_settings']['processes_sub_model_part_list'] = []
        settings['constraints_process_list'] = []
        settings['transfer_results_process']['add_previous_results'] = False
        settings['loads_process_list'] = [self._load([0., -10., 0.])]
        if increment == 'horizontal':
            settings['loads_process_list'].append(self._load([2., -10., 0.], 2.5))
        if checkpoint:
            settings['output_processes'] = {'restart': [{
                'python_module': 'save_restart_process', 'kratos_module': 'KratosMultiphysics',
                'Parameters': {'model_part_name': 'MainModelPart', 'restart_save_frequency': 0.,
                               'output_path': 'checkpoints', 'save_restart_files_in_folder': True}}]}
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
        if source:
            settings['loads_process_list'].append({
                'python_module': 'selective_initial_stress_process',
                'kratos_module': 'KratosMultiphysics.DamApplication',
                'Parameters': {
                    'model_part_name': 'MainModelPart.Parts_Concrete',
                    'source_model_part_name': 'Parts_Concrete',
                    'source_root_model_part_name': 'MainModelPart',
                    'restart_settings': {
                        'input_filename': str(source / 'MainModelPart'),
                        'input_output_path': 'checkpoints',
                        'restart_load_file_label': '1.0'}}})
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

    def _selective(self, thermal=False):
        with tempfile.TemporaryDirectory(prefix='dam_selective_') as tmp:
            root = Path(tmp)
            increment = 'thermal' if thermal else 'horizontal'
            reference = self._run(root / 'reference', self._parameters(thermal, increment, True), thermal)
            destination = self._run(root / 'new', self._parameters(
                thermal, increment, source=root / 'reference', reference_shift=50.),
                thermal, reference_shift=50.)
            base = reference['steps'][0]
            initial = destination['initial']
            self.assertEqual(initial['time'], 0.)
            self.assertFalse(initial['restarted'])
            self.assertEqual(initial['step'], 0)
            self.assertEqual(base['step'], 41)
            for key in ('displacement', 'reaction', 'velocity', 'acceleration', 'previous_displacement'):
                self._compare(initial[key], [[0., 0., 0.]] * 6, atol=0.)
            self._compare(initial['coordinates'], reference['initial']['coordinates'], atol=0.)
            self._compare(initial['residual'][4:], [0.] * 8)
            self.assertAlmostEqual(sum(base['external'][1::2]), -40., delta=1.e-10)
            self.assertAlmostEqual(sum(r[1] for r in base['reaction']), 40., delta=1.e-10)
            for key in ('stress', 'internal', 'external', 'residual'):
                error = self._compare(initial[key], base[key])
                print('selective {} initial {} error {:.6e}'.format(increment, key, error))
            for actual, expected in zip(destination['steps'], reference['steps']):
                expected_increment = [[a-b for a, b in zip(u, u0)]
                                      for u, u0 in zip(expected['displacement'], base['displacement'])]
                error = self._compare(actual['displacement'], expected_increment)
                stress_error = self._compare(actual['stress'], expected['stress'])
                reaction_error = self._compare(actual['reaction'], expected['reaction'])
                print('  stress error {:.6e}; reaction error {:.6e}'.format(stress_error, reaction_error))
                self._compare(actual['residual'][4:], [0.] * 8)
                print('selective {} t={} displacement increment error {:.6e}'.format(increment, actual['time'], error))
            self._compare(destination['steps'][0]['displacement'], [[0., 0., 0.]] * 6)
            self.assertGreater(max(abs(v) for row in base['displacement'] for v in row), 1.e-4)
            self.assertGreater(max(abs(v) for row in destination['steps'][-1]['displacement'] for v in row), 1.e-6)
            output = (root / 'new/result.post.res').read_text()
            header = 'Result "CAUCHY_STRESS_TENSOR" "Kratos" 0 Matrix OnGaussPoints'
            self.assertIn(header, output)
            displacement_frame = output.split('Result "DISPLACEMENT" "Kratos" 0 Vector OnNodes', 1)[1]
            displacement_rows = displacement_frame.split('Values\n', 1)[1].split('End Values', 1)[0].splitlines()
            self.assertEqual(len(displacement_rows), 6)
            for row in displacement_rows:
                self._compare([float(v) for v in row.split()[1:]], [0., 0., 0.], atol=0.)
            lines = output[output.index(header):].splitlines()[2:]
            tensors = []
            for line in lines:
                if line.strip() == 'End Values':
                    break
                values = line.split()
                tensors.append([float(v) for v in (values[1:] if len(values) == 7 else values)])
            self.assertEqual(len(tensors), 8)
            for tensor, stress in zip(tensors, base['stress']):
                self._compare([tensor[0], tensor[1], tensor[3]], stress, atol=5.e-5)
            if thermal:
                for state in [initial] + destination['steps']:
                    self._compare(state['reference_temperature'], [70., 70., 72., 72., 74., 74.])
                # The thermal solver computes the middle row, not the prescribed boundary process.
                temperature = destination['steps'][-1]['temperature']
                self.assertGreater(temperature[2], 72.)
                self.assertLess(temperature[2], 77.)
            # Optional retained numerical evidence, outside the source/checkpoint directories.
            if os.environ.get('DAM_SELECTIVE_RESULTS'):
                evidence = Path(os.environ['DAM_SELECTIVE_RESULTS'])
                evidence.mkdir(parents=True, exist_ok=True)
                (evidence / (increment+'.json')).write_text(json.dumps(
                    {'reference': reference, 'destination': destination}, indent=2))

    def test_selective_gravity_stresses(self):
        self._selective()

    def test_selective_gravity_and_temperature_increment(self):
        self._selective(thermal=True)

    def test_rejects_unsupported_states_and_requires_gravity(self):
        import test_dam_elastic_initial_state_transfer as elastic_tests
        import KratosMultiphysics.StructuralMechanicsApplication as SMA
        fixture = elastic_tests.TestDamElasticInitialStateTransfer()
        source_model, source, source_strategy = fixture._body(quadrilaterals=True)
        source_strategy.Solve()
        for defect, message in (('mesh', 'coordinates differ'),
                                ('material', 'material mismatch'),
                                ('rule', 'integration rule'),
                                ('law', 'unsupported constitutive law'),
                                ('displacement', 'nonzero displacement')):
            with self.subTest(defect=defect):
                model, target, strategy = fixture._body(quadrilaterals=True)
                if defect == 'mesh':
                    target.Nodes[6].X0 += 0.1
                elif defect == 'material':
                    target.Properties[1][KM.YOUNG_MODULUS] *= 2.
                elif defect == 'displacement':
                    target.Nodes[6].SetSolutionStepValue(KM.DISPLACEMENT, [0.1, 0., 0.])
                else:
                    if defect == 'rule':
                        target.Properties[1][KM.KratosGlobals.GetVariable("INTEGRATION_ORDER")] = 1
                    else:
                        target.Properties[1][KM.CONSTITUTIVE_LAW] = SMA.LinearElasticPlaneStress2DLaw()
                    for element in target.Elements:
                        element.Initialize(target.ProcessInfo)
                with self.assertRaisesRegex(RuntimeError, message):
                    Dam.SelectiveInitialStressUtility.AssignInitialStresses(source, target)
                if defect != 'displacement':
                    for element in fixture._values(target, KM.CAUCHY_STRESS_VECTOR):
                        for point in element:
                            self._compare(point, [0., 0., 0.])
        model, target, strategy = fixture._body(quadrilaterals=True)
        Dam.SelectiveInitialStressUtility.AssignInitialStresses(source, target)
        with self.assertRaisesRegex(RuntimeError, 'already has an initial state'):
            Dam.SelectiveInitialStressUtility.AssignInitialStresses(source, target)
        for node in target.Nodes:
            node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, [0., 0., 0.])
        self.assertGreater(max(abs(v) for v in fixture._residual(target)[4:]), 1.)

    def test_rejects_combination_with_previous_results_transfer(self):
        with tempfile.TemporaryDirectory(prefix='dam_selective_conflict_') as tmp:
            settings = self._parameters(source=Path(tmp) / 'unused_source')
            settings['transfer_results_process']['add_previous_results'] = True
            self._run(Path(tmp) / 'run', settings,
                      expected_error='cannot be combined with')

    def test_does_not_modify_source_state(self):
        import test_dam_elastic_initial_state_transfer as elastic_tests
        fixture = elastic_tests.TestDamElasticInitialStateTransfer()
        source_model, source, source_strategy = fixture._body(quadrilaterals=True)
        source_strategy.Solve()
        before = snapshot(source)
        model, target, strategy = fixture._body(quadrilaterals=True)
        Dam.SelectiveInitialStressUtility.AssignInitialStresses(source, target)
        self.assertEqual(snapshot(source), before)

    def test_destination_keeps_its_own_state_fields(self):
        import test_dam_elastic_initial_state_transfer as elastic_tests
        fixture = elastic_tests.TestDamElasticInitialStateTransfer()
        source_model, source, source_strategy = fixture._body(quadrilaterals=True)
        source_strategy.Solve()
        for node in source.Nodes:
            node.SetSolutionStepValue(KM.VELOCITY, [13., 14., 15.])
            node.SetSolutionStepValue(KM.ACCELERATION, [16., 17., 18.])
        source.ProcessInfo[KM.TIME] = 99.
        source.ProcessInfo[KM.STEP] = 11

        model, target, strategy = fixture._body(quadrilaterals=True)
        for node in target.Nodes:
            node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, [1., 2., 3.])
            node.SetSolutionStepValue(KM.REACTION, [4., 5., 6.])
            node.SetSolutionStepValue(KM.VELOCITY, [7., 8., 9.])
            node.SetSolutionStepValue(KM.ACCELERATION, [10., 11., 12.])
        target.ProcessInfo[KM.TIME] = 42.
        target.ProcessInfo[KM.STEP] = 7

        Dam.SelectiveInitialStressUtility.AssignInitialStresses(source, target)

        for node in target.Nodes:
            self._compare(node.GetSolutionStepValue(KM.VOLUME_ACCELERATION), [1., 2., 3.])
            self._compare(node.GetSolutionStepValue(KM.REACTION), [4., 5., 6.])
            self._compare(node.GetSolutionStepValue(KM.VELOCITY), [7., 8., 9.])
            self._compare(node.GetSolutionStepValue(KM.ACCELERATION), [10., 11., 12.])
            self.assertEqual(node.GetSolutionStepValue(KM.TEMPERATURE), 30.)
            self.assertEqual(node.GetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE), 20.)
        self.assertEqual(target.ProcessInfo[KM.TIME], 42.)
        self.assertEqual(target.ProcessInfo[KM.STEP], 7)

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
