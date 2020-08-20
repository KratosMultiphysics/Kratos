class ExaquteExecutionConstraints(object):
    def __init__(self, nodes=1, cpus_per_task=1, amount_executions=1):
        self.nodes = nodes
        self.cpus_per_task = cpus_per_task
        self.amount_executions = amount_executions
