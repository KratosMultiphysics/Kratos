import math
import matplotlib.pyplot as plt
import h5py

do_print_max_error = True
mat_deriv_recovery_types = [1, 3, 4, 6, 7]
laplacian_recovery_types = [1, 3, 4, 6, 7]
field_identifiers = ['ethier', 'sines']
mesh_regularities = [True, False]
mesh_tags = ['Altair', 'Kratos']

meshes = []
recovery_cases = []
fields = []

class Field:
    def __init__(self, identifier):
        self.identifier = identifier
        if identifier == 'ethier':
            self.characteristic_length = 4
        else:
            self.characteristic_length = 2

class Mesh:
    def __init__(self, is_regular, tag):
        self.is_regular = is_regular
        self.tag = tag

class RecoveryCase:
    def __init__(self, m_or_l, method):
        self.m_or_l = m_or_l
        self.method = method

for identifier in field_identifiers:
    fields.append(Field(identifier))

for is_regular in mesh_regularities:
    if is_regular:
        for tag in mesh_tags:
            meshes.append(Mesh(is_regular, tag))
    else:
        meshes.append(Mesh(is_regular,''))

for method in mat_deriv_recovery_types:
    recovery_cases.append(RecoveryCase('M', method))

for method in laplacian_recovery_types:
    recovery_cases.append(RecoveryCase('L', method))

def GetCurveCharacteristics(mat_deriv_or_laplacian, method):
    if mat_deriv_or_laplacian == 'M':
        if method == 1:
            type_name = 'standard'
            color = 'r'
        elif method == 2:
            type_name = 'Zhang and Naga'
            color = 'k'
        elif method == 3:
            type_name = 'L2-lumped'
            color = 'b'
        elif method == 4:
            type_name = 'L2'
            color = 'g'
        elif method == 5:
            type_name = 'L2 only gradient'
            color = 'c'
        elif method == 6:
            type_name = 'Pouliot et al.'
            color = 'brown'
        elif method == 7:
            type_name = 'Zhang and Naga'
            color = 'm'
    elif mat_deriv_or_laplacian == 'L':
        if method == 1:
            type_name = 'standard'
            color = 'r'
        elif method == 2:
            type_name = 'Zhang and Naga'
            color = 'm'
        elif method == 3:
            type_name = 'L2-lumped' # L2-div. of gradient from L2-lumped
            color = 'b'
        elif method == 4:
            type_name = 'L2' # L2-div. of gradient from L2
            color = 'g'
        elif method == 6:
            type_name = 'L2 + Pouliot et al.' # L2-div. of gradient from L2
            color = 'brown'
        elif method == 7:
            type_name = 'Guo et al.'
            color = 'm'
    return color, type_name

class Plotter:
    def __init__(self):
        self.marker_size = 10
        self.line_width = 1
        self.curves = []
        self.figures = []

    @staticmethod
    def PlotSlope(figure):
            plt.semilogx()
            plt.semilogy()
            x_min, x_max = figure.GetMinMax('sizes')
            error_avg_min, error_avg_max = figure.GetMinMax('average_errors')
            error_max_min, error_max_max = figure.GetMinMax('max_errors')
            error_min = min(error_avg_min, error_max_min), max(error_avg_max, error_max_max)[0]
            plt.xlim([0.5 * x_min, 2 * x_max])
            sizes = [x_min * 2 ** i for i in range(3)]
            expected_order = figure.GetExpectedOrder()
            slope = [0.5 * error_min * (size / sizes[0]) ** expected_order for size in sizes]

            plt.plot(sizes, slope, color='k', linestyle='dashed', label='slope = ' + str(expected_order))

    def ReadCurveData(self, curve):
        if curve.m_or_l == 'M':
            recovery_tag = 'material derivative/method = ' + str(curve.method)
        else:
            recovery_tag = 'laplacian/method = ' + str(curve.method)
        if curve.is_regular_mesh:
            mesh_identifier = 'regular mesh (' + curve.mesh_tag + ')'
        else:
            mesh_identifier = 'irregular mesh'

        data_path = curve.field_id + '/' + recovery_tag + '/' + mesh_identifier

        with h5py.File('errors_recorded/recovery_errors.hdf5', 'r') as f:
            for dset in f[data_path].values():
                size = float(dset.attrs['mesh size'])
                curve.sizes.append(size / curve.field_length_scale)
                curve.average_errors.append(dset[0])
                curve.max_errors.append(dset[1])

        self.curves.append(curve)

    def MakePlots(self):
        self.ClassifyCurvesIntoDifferentFigures()

        for figure in self.figures:
            figure.fig = plt.figure(figure.number, figsize=(12,10))

            for curve in figure.curves:
                figure.CalculateSlopes()
                color, derivative_type = GetCurveCharacteristics(curve.m_or_l, curve.method)
                curve.sizes, curve.average_errors, curve.max_errors = \
                zip(*sorted(zip(curve.sizes, curve.average_errors, curve.max_errors)))
                plt.plot(curve.sizes,
                         curve.average_errors,
                         marker=curve.average_marker_type,
                         ms=self.marker_size,
                         color=color,
                         label=derivative_type + curve.slope_msg,
                         linewidth=self.line_width,
                         linestyle='solid',
                         markersize=20)

                if do_print_max_error:
                    color, derivative_type = GetCurveCharacteristics(curve.m_or_l, curve.method)
                    plt.plot(curve.sizes,
                             curve.max_errors,
                             marker=curve.average_marker_type,
                             markeredgecolor=color,
                             markerfacecolor='None',
                             ms=self.marker_size,
                             color=color,
                             label= '” ' + curve.max_slope_msg,
                             linewidth=self.line_width,
                             linestyle='solid',
                             markersize=20)

            Plotter.PlotSlope(figure)
            plt.xlabel('$h$', fontsize=30)
            plt.ylabel('$E_2$', fontsize=30)
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            legend = plt.legend(loc='lower right', prop={'size': 15})
            for handle in legend.legendHandles:
                handle._legmarker.set_markersize(6)

            plt.savefig(figure.title, format='pdf', bbox_inches='tight')

    def ClassifyCurvesIntoDifferentFigures(self):
        figure_correspondences = dict()
        i_new_figure = 0
        for curve in self.curves:
            combination = [curve.field_id, curve.m_or_l, curve.mesh_tag]
            if str(combination) not in figure_correspondences.keys():
                new_figure = Figure(combination)
                new_figure.Add(curve)
                self.figures.append(new_figure)
                figure_correspondences[str(combination)] = i_new_figure
                i_new_figure += 1
            else:
                i_figure = figure_correspondences[str(combination)]
                self.figures[i_figure].Add(curve)

class Figure:
    total_number_of_figures = 0

    @classmethod
    def CalculateLastSlopes(self, sizes, results):
        delta_result = math.log(results[-1] / results[-2])
        delta_size   = math.log(sizes[-1] / sizes[-2])
        slope = abs(delta_result/delta_size)
        return slope

    def __init__(self, combination):
        Figure.total_number_of_figures += 1
        self.combination = combination
        self.number = Figure.total_number_of_figures
        self.curves = []
        self.SetTitle()
        self.SetOrder()

    def Add(self, curve):
        self.curves.append(curve)

    def SetTitle(self):
        field_id, m_or_l, mesh_tag = self.combination

        if m_or_l == 'M':
            self.title = 'M_'
        else:
            self.title = 'L_'

        self.title += field_id

        if mesh_tag == '':
            self.title += '_irregular'
        else:
            self.title += '_regular_' + mesh_tag

        self.title += '.pdf'

    def SetOrder(self):
        m_or_l = self.combination[1]
        if m_or_l == 'M':
            self.expected_order = 2
        else:
            self.expected_order = 1

    def GetExpectedOrder(self):
        return self.expected_order

    def GetMinMax(self, what='sizes'):
        minimum = min(curve.GetMinMax(what)[0] for curve in self.curves)
        maximum = max(curve.GetMinMax(what)[1] for curve in self.curves)
        return minimum, maximum

    def CalculateSlopes(self):
        for curve in self.curves:
            try:
                curve.slope = Figure.CalculateLastSlopes(curve.sizes, curve.average_errors)
                curve.slope_msg = ' (' + str(round(curve.slope, 2)) + ')'
            except:
                curve.slope_msg = ''
            if do_print_max_error:
                try:
                    curve.max_slope = Figure.CalculateLastSlopes(curve.sizes, curve.max_errors)
                    curve.max_slope_msg = ' (' + str(round(curve.max_slope, 2)) + ')'
                except:
                    curve.max_slope_msg = ''

class Curve:
    def __init__(self, field, mesh, case):
        self.field_id = field.identifier
        self.field_length_scale = field.characteristic_length
        self.mesh_tag = mesh.tag
        self.is_regular_mesh = mesh.is_regular
        self.m_or_l = case.m_or_l
        self.method = case.method

        if self.m_or_l == 'M':
            self.average_marker_type = 'v'
            self.expected_order = 2
        else:
            self.average_marker_type = '^'
            self.expected_order = 1

        self.sizes = []
        self.average_errors = []
        self.max_errors = []
        self.values = dict()
        self.values['sizes'] = self.sizes
        self.values['average_errors'] = self.average_errors
        self.values['max_errors'] = self.max_errors

    def GetMinMax(self, what='sizes'):
        values = self.values[what]
        return min(values), max(values)

plotter = Plotter()

for field in fields:
    for mesh in meshes:
        for case in recovery_cases:
            curve = Curve(field, mesh, case)
            plotter.ReadCurveData(curve)

plotter.MakePlots()
