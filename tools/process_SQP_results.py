"""Process simulation results."""

import sys
import os
import pathlib
import mmap

import numpy
import pandas

constraint_solvers = (
    "QP",
    "SQP",
)
timesteps = (1e-2, 1e-3, 1e-4, 1e-5)
constraint_offsets = (1e-2, 1e-3, 1e-4, 1e-5)
constraint_types = ("graphics", "volume", "Verschoor")
qp_solvers = (
    # "OSQP",
    "Gurobi",
)
active_set_convergence_policies = (
    # "noActiveSetConvergence",
    "useActiveSetConvergence",
)


def save_results_csv(results):
    """Save results to seperate CSV files."""
    with open("results-SQP.csv", "w", newline="") as f:
        scene_names = sorted(list(results.keys()))
        for constraint_solver in constraint_solvers:
            f.write(f"{constraint_solver}\n")
            for constraint_type in constraint_types:
                f.write(f"{constraint_type}\n")
                for qp_solver in qp_solvers:
                    f.write(f"{qp_solver}\n")
                    for active_set_convergence_policy in active_set_convergence_policies:
                        f.write(f"{active_set_convergence_policy}\n")
                        for scene_name in scene_names:
                            f.write(f"{scene_name}\n")
                            results[scene_name][constraint_solver][constraint_type][qp_solver][active_set_convergence_policy].to_csv(
                                f, header=False, index=False)
                            f.write("\n")


def check_error_file(log_path):
    err_log = log_path.with_suffix("").with_suffix(".err.txt")
    if not err_log.exists():
        err_log = log_path.with_suffix(".err")
        if not err_log.exists():
            return "Incomplete"
    with open(err_log) as err_file:
        s = mmap.mmap(err_file.fileno(), 0, access=mmap.ACCESS_READ)
        if s.find(b"out-of-memory") != -1 or s.find(b"what():  vector::reserve") != -1:
            return "Out-of-Memory"
        elif s.find(b"TIME LIMIT") != -1:
            return "Timeout"
        elif s.find(b"GRBException") != -1:
            return "GRBException"
    return "Incomplete"


def get_sim_status(log_path):
    with open(log_path) as log_file:
        s = mmap.mmap(log_file.fileno(), 0, access=mmap.ACCESS_READ)
        if s.find(b"intersecting state") != -1:
            return "Intersecting"
        if s.find(b"blow-up") != -1:
            return "Blow-Up"
        if s.find(b"simulation finished") != -1:
            return "Pass"
    error_status = check_error_file(log_path)
    if error_status != "Timeout":
        print("{:s}: {}".format(error_status, log_path))
    if error_status == "Out-of-Memory":
        if s.find(b"Gurobi_status=12") != -1:
            return "Gurobi_status=12"
        if s.find(b"Gurobi_status=12") != -1:
            return "Gurobi_status=4"
    return error_status


def add_scene_to_results(results, scene_name, default_result):
    if scene_name in results:
        return
    results[scene_name] = {}
    for constraint_solver in constraint_solvers:
        results[scene_name][constraint_solver] = {}
        for constraint_type in constraint_types:
            results[scene_name][constraint_solver][constraint_type] = {}
            for qp_solver in qp_solvers:
                results[scene_name][constraint_solver][constraint_type][qp_solver] = {}
                for active_set_convergence_policy in active_set_convergence_policies:
                    results[scene_name][constraint_solver][constraint_type][qp_solver][active_set_convergence_policy] = default_result.copy()


def main():
    """Process simulation results."""
    vals = numpy.full((len(timesteps), len(constraint_offsets)), "")
    default_result = pandas.DataFrame(
        vals, index=timesteps, columns=constraint_offsets)

    results = {}
    results_path = pathlib.Path(sys.argv[1])
    for log_dir in results_path.glob("**/logs"):
        for log in log_dir.iterdir():
            if (log.with_suffix("").suffix != ".out"
                    or os.stat(log).st_size == 0
                    or str(log.name)[:2] == "._"):
                continue
            # meshCO_pointTriangleCO__SQP_1e-2_1e-2_graphics_Gurobi_useActiveSetConvergence.out.txt
            scene_name, params = log.name.split("__")
            params = params.split(".")[0].split("_")
            if params[0] == "IP":
                continue
            (constraint_solver, timestep, constraint_offset, constraint_type,
             qp_solver, active_set_convergence_policy) = params
            timestep = float(timestep)
            constraint_offset = float(constraint_offset)
            if scene_name not in results:
                add_scene_to_results(results, scene_name, default_result)
            results[scene_name][constraint_solver][constraint_type][qp_solver][
                active_set_convergence_policy][constraint_offset][timestep] = get_sim_status(log)
    save_results_csv(results)


if __name__ == "__main__":
    main()
