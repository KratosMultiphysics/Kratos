
from typing import Any
import matplotlib.pyplot as plt

def set_default_graph_params(width=8, height=6):
    print("graph width = ", width)
    print("graph height = ", height)
    params = {
        'text.usetex': True,
        'font.size': 10,
        'text.latex.unicode': True,
        'figure.titlesize': 10,
        'figure.figsize': (width, height),
        'figure.dpi': 300,
        'axes.titlesize': 10,
        'axes.labelsize': 10,
        'axes.grid': 'True',
        'axes.grid.which': 'both',
        'axes.xmargin': 0.05,
        'axes.ymargin': 0.05,
        'lines.linewidth': 1,
        'lines.markersize': 3,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'ytick.minor.visible': 'true',
        'xtick.minor.visible': 'true',
        'grid.linestyle': '-',
        'grid.linewidth': 0.5,
        'grid.alpha': 0.3,
        'legend.fontsize': 10,
        'legend.loc': 'lower right',
        'savefig.dpi': 600,
        'savefig.format': 'png',
        'savefig.bbox': 'tight'
    }

    plt.rcParams.update(params)

def ReplaceEngineeringNotation(text: str) -> str:
    index = text.find("e-")
    while (index > -1):
        end_index = index
        for i in range(index + 2, len(text)):
            try:
                _temp = int(text[i])
            except:
                end_index = i
                break

        text = text[:index] + r"\times 10^{-" + str(
            int(text[index + 2:end_index])) + "}" + text[end_index:]
        index = text.find("e-")

    return text

class GraphSettings:
    def __init__(self, number_of_total_graphs: int, pre_defined_line_style: str = "", pre_defined_marker_style: str = "", pre_defined_color: str = "", order_of_change: str = "lmc", always_change_color=True):
        self.graph_index = 0
        self.number_of_total_graphs = number_of_total_graphs

        self.default_params_dict = {
            "linestyle": pre_defined_line_style,
            "marker": pre_defined_marker_style,
            "color": pre_defined_color
        }

        self.unique_order_of_change = self.__GetOrderedList(order_of_change)
        if always_change_color:
            self.max_possibilities_map = {"l": (int(number_of_total_graphs / 7) + 1) % 14, "c": 7, "m": 7}
            self.frequency_map = {"l": 1, "c": 1, "m": self.max_possibilities_map["l"]}
        else:
            self.max_possibilities_map = {"l": 13, "c": 7, "m": 7}
            if pre_defined_line_style != "":
                self.unique_order_of_change = self.unique_order_of_change.replace("l", "")

            if pre_defined_marker_style != "":
                self.unique_order_of_change = self.unique_order_of_change.replace("m", "")

            if pre_defined_color != "":
                self.unique_order_of_change = self.unique_order_of_change.replace("c", "")

            self.frequency_map = {"l": 1, "c": 1, "m": 1}
            for i, unique_char in enumerate(self.unique_order_of_change[1:]):
                previous_unique_char = self.unique_order_of_change[i]
                self.frequency_map[unique_char] = self.frequency_map[previous_unique_char] * self.max_possibilities_map[previous_unique_char]

    def GetGraphParameters(self) -> dict[str, Any]:
        colors = ["r", "g", "b", "c", "m", "y", "k"]
        markers = ["o", "v", "^", "s", "D", "<", ">"]

        params_dict = dict(self.default_params_dict)
        if params_dict["linestyle"] == "":
            params_dict["linestyle"] = self.__GetLineStyle(int(self.graph_index / self.frequency_map["l"]) % self.max_possibilities_map["l"])
        if params_dict["color"] == "":
            params_dict["color"] = colors[int(self.graph_index / self.frequency_map["c"]) % self.max_possibilities_map["c"]]
        if params_dict["marker"] == "":
            params_dict["marker"] = markers[int(self.graph_index / self.frequency_map["m"]) % self.max_possibilities_map["m"]]

        self.graph_index += 1
        return params_dict

    @staticmethod
    def __GetLineStyle(index: int):
        line_style_tuple = [
            ('solid'                , 'solid'),  # Same as (0, ()) or '-'           #0
            ('loosely dashed'       , (0, (5, 10))),                       #1
            ('dashed'               , (0, (5, 5))),                                #2
            ('densely dashed'       , (0, (5, 1))),                        #3
            ('dashdotted'           , (0, (3, 5, 1, 5))),                      #5
            ('densely dashdotted'   , (0, (3, 1, 1, 1))),              #6
            ('dotted'               , (0, (1, 1))),                                #8
            ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1))),      #12
            ('loosely dotted'       , (0, (1, 10))),                       #7
            ('loosely dashdotted'   , (0, (3, 10, 1, 10))),            #4
            ('densely dotted'       , (0, (1, 1))),                        #9
            ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),  #11
            ('dashdotdotted'        , (0, (3, 5, 1, 5, 1, 5))),             #10
        ]
        return line_style_tuple[index][1]

    @staticmethod
    def __GetOrderedList(input_str: str) -> str:
        result = ""
        for character in input_str:
            if result.find(character) == -1:
                result += character
        return result

if __name__ == "__main__":
    gs = GraphSettings(13, pre_defined_marker_style="o", pre_defined_color="r")
    # import numpy

    # x = numpy.arange(-1, 1, 0.1)
    # for i in range(13):
    #     plt.plot(x, x+i, **gs.GetGraphParameters(), label=f"{i}")
    # plt.legend(loc="upper right")
    # plt.show()

    # test_string = "lccmml"
    # print(GraphSettings._GetOrderedList(test_string))
