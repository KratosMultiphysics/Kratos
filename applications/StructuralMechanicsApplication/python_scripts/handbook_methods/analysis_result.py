from dataclasses import dataclass, field

@dataclass
class AnalysisResult:
    method_name: str
    category: str
    value: float
    output_variable: object
    metadata: dict = field(default_factory=dict)
