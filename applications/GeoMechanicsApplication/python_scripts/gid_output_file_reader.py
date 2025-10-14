import math


class GiDOutputFileReader:
    def __init__(self):
        self._reset_internal_state()

    def read_output_from(self, gid_output_file_path):
        self._reset_internal_state()

        with open(gid_output_file_path, "r") as result_file:
            for line in result_file:
                line = line.strip()
                if line.startswith("GaussPoints"):
                    self._process_begin_of_gauss_points(line)
                elif line == "End GaussPoints":
                    self._process_end_of_gauss_points(line)
                elif line.startswith("Result"):
                    self._process_result_header(line)
                elif line == "Values":
                    self._process_begin_of_block(line)
                elif line == "End Values":
                    self._process_end_of_block(line)
                elif self.current_block_name == "GaussPoints":
                    self._process_gauss_point_data(line)
                elif self.current_block_name == "Values":
                    self._process_value_data(line)

        return self.output_data

    def _reset_internal_state(self):
        self.output_data = {}
        self.current_block_name = None
        self.result_name = None
        self.result_type = None
        self.result_location = None
        self.gauss_points_name = None
        self.current_integration_point = 0

    def _process_begin_of_gauss_points(self, line):
        self._process_begin_of_block(line)
        self.gauss_points_name = self._strip_off_quotes(line.split()[1])
        if self.current_block_name not in self.output_data:
            self.output_data[self.current_block_name] = {}
        self.output_data[self.current_block_name][self.gauss_points_name] = {}

    def _process_end_of_gauss_points(self, line):
        self._process_end_of_block(line)
        self.gauss_points_name = None

    def _process_result_header(self, line):
        if "results" not in self.output_data:
            self.output_data["results"] = {}
        words = line.split()
        self.result_name = self._strip_off_quotes(words[1])
        if self.result_name not in self.output_data["results"]:
            self.output_data["results"][self.result_name] = []
        self.result_type = words[4]
        self.result_location = words[5]
        this_result = {
            "time": float(words[3]),
            "location": self.result_location,
            "values": [],
        }
        self.output_data["results"][self.result_name].append(this_result)
        if self.result_location == "OnGaussPoints":
            self.current_integration_point = 0
            self.gauss_points_name = self._strip_off_quotes(words[6])

    def _process_gauss_point_data(self, line):
        if line.startswith("Number Of Gauss Points:"):
            pos = line.index(":")
            num_gauss_points = int(line[pos + 1 :].strip())
            self.output_data[self.current_block_name][self.gauss_points_name][
                "size"
            ] = num_gauss_points

    def _process_value_data(self, line):
        words = line.split()
        if self.result_location == "OnNodes":
            self._process_nodal_result(words)
        elif self.result_location == "OnGaussPoints":
            self._process_gauss_point_result(words)

    def _process_nodal_result(self, words):
        value = {"node": int(words[0])}
        if self.result_type == "Scalar":
            value["value"] = float(words[1])
        elif self.result_type == "Vector" or self.result_type == "Matrix":
            value["value"] = [float(x) for x in words[1:]]
        self.output_data["results"][self.result_name][-1]["values"].append(value)

    def _process_gauss_point_result(self, words):
        self.current_integration_point %= self.output_data["GaussPoints"][
            self.gauss_points_name
        ]["size"]
        self.current_integration_point += 1
        if self.current_integration_point == 1:
            value = {"element": int(words[0]), "value": []}
            self.output_data["results"][self.result_name][-1]["values"].append(value)
            words.pop(0)

        value = self.output_data["results"][self.result_name][-1]["values"][-1]["value"]
        if self.result_type == "Scalar":
            value.append(float(words[0]))
        elif self.result_type == "Matrix" or self.result_type == "Vector":
            value.append([float(x) for x in words])
        else:
            raise RuntimeError(f'Unsupported result type "{self.result_type}"')

    def _process_begin_of_block(self, line):
        assert self.current_block_name is None  # nested blocks are not supported
        self.current_block_name = line.split()[0]

    def _process_end_of_block(self, line):
        words = line.split()
        assert words[0] == "End"
        assert self.current_block_name == words[1]
        self.current_block_name = None

    def _strip_off_quotes(self, quoted_string):
        assert quoted_string[0] == '"'
        assert quoted_string[-1] == '"'
        return quoted_string[1:-1]

    @staticmethod
    def nodal_values_at_time(result_item_name, time, output_data, node_ids=None):
        matching_item = None
        for item in output_data["results"][result_item_name]:
            if math.isclose(item["time"], time):
                matching_item = item
                break
        if matching_item is None:
            raise RuntimeError(
                f"'{result_item_name}' does not have results at time {time}"
            )

        if matching_item["location"] != "OnNodes":
            raise RuntimeError(f"'{result_item_name}' is not a nodal result")

        node_id_to_value_map = {
            item["node"]: item["value"] for item in matching_item["values"]
        }

        if node_ids is None:  # return all values
            node_ids = [item["node"] for item in matching_item["values"]]

        return [node_id_to_value_map[node_id] for node_id in node_ids]

    @staticmethod
    def element_integration_point_values_at_time(
        result_item_name,
        time,
        output_data,
        element_ids=None,
        integration_point_indices=None,
    ):
        if element_ids and element_ids != sorted(element_ids):
            raise RuntimeError("Element IDs must be sorted")

        matching_item = None
        for item in output_data["results"][result_item_name]:
            if math.isclose(item["time"], time):
                matching_item = item
                break
        if matching_item is None:
            raise RuntimeError(
                f"'{result_item_name}' does not have results at time {time}"
            )

        if matching_item["location"] != "OnGaussPoints":
            raise RuntimeError(
                f"'{result_item_name}' is not an integration point result"
            )

        if element_ids:
            element_results = [
                item["value"]
                for item in matching_item["values"]
                if item["element"] in element_ids
            ]
        else:
            element_results = [item["value"] for item in matching_item["values"]]

        if integration_point_indices:
            result = []
            for element_result in element_results:
                result.append(
                    [
                        item
                        for index, item in enumerate(element_result)
                        if index in integration_point_indices
                    ]
                )
            return result
        else:
            return element_results
