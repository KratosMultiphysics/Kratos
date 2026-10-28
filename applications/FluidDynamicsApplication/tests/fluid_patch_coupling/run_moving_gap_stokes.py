#     |  /           |
#     ' /   __| _` | __|  _ \   __|
#     . \  |   (   | |   (   |\__ `
#    _|\_\_|  \__,_|\__|\___/ ____/
#                    Multi-Physics
#
#   Fluid Dynamics Application
#
#   License:         BSD License
#                    Kratos default license: kratos/license.txt
#
#   Main author:     Nicolò Antonelli
#

"""Run the translating Q1 Stokes MMS and save per-step diagnostics as JSON."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from test_moving_fluid_patch_gap import solve_moving_gap


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--divisions", type=int, default=12)
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--steps", type=int, default=6)
    parser.add_argument("--speed", type=float, default=0.01)
    parser.add_argument("--penalty", type=float, default=0.0)
    parser.add_argument("--angle", type=float, default=0.0, help="Initial fixed orientation; no angular motion.")
    parser.add_argument("--nonlinear", action="store_true", help="Use the cubic spatial velocity MMS.")
    parser.add_argument("--history-mode", choices=("exact", "transfer", "compare"),
                        default="transfer", help="Transfer permits bounded extrapolation; compare runs both modes.")
    parser.add_argument("--output", type=Path, default=Path("moving_gap_stokes.json"))
    args = parser.parse_args()
    if args.divisions < 8 or args.steps < 1 or args.dt <= 0:
        parser.error("Use divisions >= 8, steps >= 1, dt > 0.")
    kwargs = dict(divisions=args.divisions, dt=args.dt, steps=args.steps, speed=args.speed,
                  coefficient=args.penalty, nonlinear=args.nonlinear, angle=args.angle)
    if args.history_mode == "compare":
        rows = {"exact": solve_moving_gap(**kwargs, history_mode="exact"),
                "transfer": solve_moving_gap(
                    **kwargs, history_mode="transfer")}
    else:
        rows = solve_moving_gap(**kwargs, history_mode=args.history_mode)
    args.output.write_text(json.dumps(rows, indent=2)+"\n")
    print(f"Diagnostics saved in {args.output.resolve()}")


if __name__ == "__main__":
    main()
