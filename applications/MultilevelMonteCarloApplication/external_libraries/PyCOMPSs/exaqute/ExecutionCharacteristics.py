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
class ExecutionCharacteristics(object):
    def __init__(self, generationNodes, generationCpus_per_node, runNodes,
                 runCpus_per_node, amountExecutions):
        self.generationNodes = generationNodes
        self.generationCpus_per_node = generationCpus_per_node
        self.runNodes = runNodes
        self.runCpus_per_node = runCpus_per_node
        self.amountExecutions = amountExecutions
