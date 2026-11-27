#!/usr/bin/env python3
"""Numerical comparison of two generated mortar condition ``.cpp`` files.

The generated ``CalculateLocalLHS`` / ``StaticCalculateLocalRHS`` bodies are plain arithmetic on
the condition data. This tool extracts, for every explicit specialisation, every node and every
active-set branch, the block of C++ statements, translates it to Python and evaluates it on
random (seeded) inputs, comparing the resulting local matrix/vector of both files. It allows
checking that a regenerated file is numerically equivalent to a previous one without compiling
(e.g. after a change of the sympy version, of the common-subexpression elimination or of the
printing), independently of the factor naming.

Usage:
    python3 compare_generated_conditions.py OLD.cpp NEW.cpp [--samples N] [--tolerance 1e-9]
"""
import argparse
import math
import random
import re
import sys

_BRANCH_MARKERS = ("// INACTIVE", "// ACTIVE", "// OBJECTIVE-SLIP", "// NONOBJECTIVE-SLIP", "// OBJECTIVE-STICK", "// NONOBJECTIVE-STICK", "// ACTIVE-SLIP", "// ACTIVE-STICK")
_SPECIALISATION = re.compile(r"^template<>\s*\nvoid (\w+)<(\d+),\s*(\d+),\s*(true|false),\s*(\d+)>::(\w+)\(", re.MULTILINE)
_NODE = re.compile(r"^\s*// NODE (\d+)\s*$")


class Block:
    """One leaf block of statements (factors and accumulations)."""

    def __init__(self):
        self.statements = []


def ExtractBlocks(text):
    """``{(class, dim, nnodes, nv, nnodes_master, method, node, branch): Block}``.

    Branches with nested ``if (is_objetive)`` are identified by the innermost marker.
    """
    blocks = {}
    for match in _SPECIALISATION.finditer(text):
        class_name, dim, nnodes, nv, nnodes_master, method = match.groups()
        key_base = (class_name, int(dim), int(nnodes), nv == "true", int(nnodes_master), method)
        body_start = text.index("{", match.end())
        body_end = text.index("\n}\n", body_start)
        body = text[body_start:body_end]
        node = None
        branch = None
        for line in body.split("\n"):
            node_match = _NODE.match(line)
            if node_match:
                node = int(node_match.group(1))
                branch = None
                continue
            marker = next((m for m in _BRANCH_MARKERS if line.rstrip().endswith(m)), None)
            if marker is not None:
                branch = marker[3:]
                continue
            stripped = line.strip()
            if node is None or branch is None:
                continue
            if stripped.startswith("const double") or stripped.startswith("rLocal"):
                blocks.setdefault(key_base + (node, branch), Block()).statements.append(stripped)
        if "rLocal" in body and not any(k[:6] == key_base for k in blocks):
            blocks[key_base + (-1, "FORWARDER")] = Block()
    return blocks


def _Translate(statement):
    """C++ statement to Python."""
    statement = statement.split("//")[0].strip()
    statement = statement.rstrip(";")
    statement = statement.replace("std::pow(", "pow(").replace("std::sqrt(", "sqrt(").replace("std::fabs(", "abs(").replace("std::abs(", "abs(")
    statement = re.sub(r"\]\((\d+),(\d+)\)", r"][\1,\2]", statement)  # Delta...[k](i,j)
    statement = re.sub(r"\b([A-Za-z]\w*)\((\d+),(\d+)\)", r"\1[\2,\3]", statement)  # Name(i,j)
    if statement.startswith("const double"):
        statement = statement[len("const double"):].strip()
    return statement


class _Matrix:
    """Minimal dense matrix with ``m[i,j]`` access."""

    def __init__(self, rows, cols, rng, scale=1.0):
        self.rows, self.cols = rows, cols
        self.data = [[rng.uniform(0.5, 1.5) * scale for _ in range(cols)] for _ in range(rows)]

    def __getitem__(self, index):
        i, j = index
        return self.data[i][j]


class _Accumulator:
    """Local matrix/vector being accumulated (only the touched entries are stored)."""

    def __init__(self):
        self.values = {}

    def __getitem__(self, index):
        return self.values.get(index, 0.0)

    def __setitem__(self, index, value):
        self.values[index] = value


def BuildInputs(dim, nnodes, nnodes_master, seed):
    rng = random.Random(seed)
    inputs = {
        "pow": pow, "sqrt": math.sqrt, "abs": abs,
        "u1": _Matrix(nnodes, dim, rng, 0.1), "u2": _Matrix(nnodes_master, dim, rng, 0.1),
        "u1old": _Matrix(nnodes, dim, rng, 0.1), "u2old": _Matrix(nnodes_master, dim, rng, 0.1),
        "X1": _Matrix(nnodes, dim, rng), "X2": _Matrix(nnodes_master, dim, rng),
        "LM": _Matrix(nnodes, dim, rng), "NormalSlave": _Matrix(nnodes, dim, rng), "TangentSlave": _Matrix(nnodes, dim, rng),
        "DOperator": _Matrix(nnodes, nnodes, rng), "MOperator": _Matrix(nnodes, nnodes_master, rng),
        "DOperatorold": _Matrix(nnodes, nnodes, rng), "MOperatorold": _Matrix(nnodes, nnodes_master, rng),
        "DeltaDOperator": [_Matrix(nnodes, nnodes, rng) for _ in range((nnodes + nnodes_master) * dim)],
        "DeltaMOperator": [_Matrix(nnodes, nnodes_master, rng) for _ in range((nnodes + nnodes_master) * dim)],
        "DeltaNormalSlave": [_Matrix(nnodes, dim, rng) for _ in range(nnodes * dim)],
        "PenaltyParameter": [rng.uniform(0.5, 1.5) for _ in range(nnodes)],
        "DynamicFactor": [rng.uniform(0.5, 1.5) for _ in range(nnodes)],
        "LMNormal": [rng.uniform(0.5, 1.5) for _ in range(nnodes)],
        "mu": [rng.uniform(0.1, 0.9) for _ in range(nnodes)],
        "ScaleFactor": rng.uniform(0.5, 1.5), "TangentFactor": rng.uniform(0.5, 1.5),
    }
    return inputs


def Evaluate(block, inputs):
    namespace = dict(inputs)
    namespace["rLocalLHS"] = _Accumulator()
    namespace["rLocalRHS"] = _Accumulator()
    code = "\n".join(_Translate(statement) for statement in block.statements)
    exec(code, namespace)
    return namespace["rLocalLHS"].values, namespace["rLocalRHS"].values


def Compare(old_text, new_text, samples=3, tolerance=1.0e-9, verbose=False, require_same_blocks=True):
    """Compare the blocks of two generated files; returns ``True`` when they are numerically equivalent.

    With ``require_same_blocks=False`` only the blocks present in both files are compared (e.g. a
    partial regeneration against the full committed file).
    """
    old_blocks = ExtractBlocks(old_text)
    new_blocks = ExtractBlocks(new_text)
    keys_old, keys_new = set(old_blocks), set(new_blocks)
    problems = 0
    if require_same_blocks and keys_old != keys_new:
        problems += 1
        print("Block sets differ:")
        for key in sorted(keys_old - keys_new):
            print("  only in old:", key)
        for key in sorted(keys_new - keys_old):
            print("  only in new:", key)
    checked = 0
    worst = 0.0
    for key in sorted(keys_old & keys_new):
        if key[-1] == "FORWARDER":
            continue
        _, dim, nnodes, _, nnodes_master, _, _, _ = key
        for sample in range(samples):
            inputs = BuildInputs(dim, nnodes, nnodes_master, 1000 * sample + 7)
            lhs_old, rhs_old = Evaluate(old_blocks[key], inputs)
            lhs_new, rhs_new = Evaluate(new_blocks[key], inputs)
            for label, a, b in (("LHS", lhs_old, lhs_new), ("RHS", rhs_old, rhs_new)):
                indices = set(a) | set(b)
                for index in indices:
                    va, vb = a.get(index, 0.0), b.get(index, 0.0)
                    error = abs(va - vb) / max(1.0, abs(va), abs(vb))
                    worst = max(worst, error)
                    if error > tolerance:
                        problems += 1
                        print("MISMATCH {} {} {}: old {:.12g} new {:.12g} (rel. error {:.3e})".format(key, label, index, va, vb, error))
        checked += 1
        if verbose:
            print("ok", key)
    print("Compared {} blocks x {} samples; worst relative error {:.3e}; {} problems".format(checked, samples, worst, problems))
    return problems == 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("old")
    parser.add_argument("new")
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--tolerance", type=float, default=1.0e-9)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--common-only", action="store_true", help="Compare only the blocks present in both files")
    args = parser.parse_args(argv)
    with open(args.old) as f:
        old_text = f.read()
    with open(args.new) as f:
        new_text = f.read()
    return 0 if Compare(old_text, new_text, args.samples, args.tolerance, args.verbose, not args.common_only) else 1


if __name__ == "__main__":
    sys.exit(main())
