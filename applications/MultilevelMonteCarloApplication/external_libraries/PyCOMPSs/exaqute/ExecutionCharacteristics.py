class ExecutionCharacteristics(object):
    def __init__(self, generationNodes, generationCpus_per_node, runNodes,
                 runCpus_per_node, amountExecutions):
        self.generationNodes = generationNodes
        self.generationCpus_per_node = generationCpus_per_node
        self.runNodes = runNodes
        self.runCpus_per_node = runCpus_per_node
        self.amountExecutions = amountExecutions
