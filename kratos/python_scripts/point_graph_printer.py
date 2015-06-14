from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
CheckForPreviousImport()


class PrintGraphPrinter:

    def __init__(self, output_nodes_list,
                 model_part, variables_dictionary, domain_size):
        self.model_part = model_part

        self.outfile_list = []
        self.point_list = []
        self.var_list = []
        self.N = []
        self.elem_list = []
        self.is_scalar_var = []

        self.point_locator = PointLocation(model_part)

        for item in output_nodes_list:
            pos = item[0]
            var = variables_dictionary[item[1]]
            shape_functions = Vector(4)

            X = pos[0]
            Y = pos[1]
            Z = pos[2]

            # find the element where the node falls
            if domain_size == 2:
                Id = self.point_locator.Find2D(X, Y, shape_functions)
            else:
                Id = self.point_locator.Find3D(X, Y, Z, shape_functions)

            found_point = self.point_locator.found()

            if(found_point):
                outfile = open(item[2], 'w')
                elem = self.model_part.Elements[Id]

                if((isinstance(elem.GetNodes()[0].GetSolutionStepValue(var, 0), int)) or (isinstance(elem.GetNodes()[0].GetSolutionStepValue(var, 0), float))):
                    self.is_scalar_var.append(True)
                else:
                    self.is_scalar_var.append(False)

                aaa = str("#output for var: ") + str(item[1]) + str(
                    " on node with coord") + str(item[0]) + "\n"
                outfile.write(aaa)

                self.N.append(shape_functions)
                self.outfile_list.append(outfile)
                self.point_list.append(pos)
                self.var_list.append(var)
                self.elem_list.append(elem)

    def PrintGraphs(self, time):

        for i in range(0, len(self.outfile_list)):
            N = self.N[i]
            outfile = self.outfile_list[i]
            var = self.var_list[i]
            elem = self.elem_list[i]

            if(self.is_scalar_var[i]):
                result = 0.0
                for node, N in zip(elem.GetNodes(), N):
                    result += N * node.GetSolutionStepValue(var, 0)

                aaa = str(time) + str(" ") + str(result) + str(" \n")
                outfile.write(aaa)
            else:
                result = Vector(3)
                result[0] = 0.0
                result[1] = 0.0
                result[2] = 0.0

                for node, N in zip(elem.GetNodes(), N):
                    a = node.GetSolutionStepValue(var, 0)
                    result[0] += N * a[0]
                    result[1] += N * a[1]
                    result[2] += N * a[2]

                aaa = str(time) + str(" ") + str(result[0]) + " " + str(
                    result[1]) + " " + str(result[2]) + str(" \n")
                outfile.write(aaa)

            outfile.flush()

    #
    def Close():
        for item in self.outfile_list:
            item.close()

    #
    def identity(item):
        return item

    #
    def first(iterable, predicate=identity):
        for item in iterable:
            if predicate(item):
                return item
        raise ValueError('No satisfactory value found')
