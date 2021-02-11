#!/usr/bin/python
#
#  Copyright 2002-2019 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
from exaqute.ExaquteTaskLocal import *


@ExaquteTask(returns=1)
def check_vector(*collection_in):
    base_string = ""
    for elem in from_args_to_vector(collection_in):
        base_string += str(len(elem))
    return base_string


def main():
    vec = [[1, 2, 3], [4, 5]]
    result = check_vector(*(from_vector_to_args(vec)))


if __name__ == "__main__":
    main()
