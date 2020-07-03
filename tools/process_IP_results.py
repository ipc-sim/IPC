"""Process simulation results."""

import sys
import os
import pathlib
import mmap

import numpy
import pandas

timesteps = (1e-2, 1e-3, 1e-4, 1e-5)


def save_results_csv(results):
    """Save results to seperate CSV files."""
    with open("results-IP.csv", "w", newline="") as f:
        scene_names = sorted(list(results.keys()))
        for scene_name in scene_names:
            f.write(f"{scene_name}\n")
            results[scene_name].to_csv(f, header=False, index=False)
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
    return check_error_file(log_path)


def add_scene_to_results(results, scene_name, default_result):
    if scene_name in results:
        return
    results[scene_name] = default_result.copy()


def main():
    """Process simulation results."""
    vals = numpy.full((len(timesteps), 1), "")
    default_result = pandas.DataFrame(vals, index=timesteps)

    results = {}
    results_path = pathlib.Path(sys.argv[1])
    for log_dir in results_path.glob("**/logs"):
        for log in log_dir.iterdir():
            if (log.with_suffix("").suffix != ".out"
                    or os.stat(log).st_size == 0
                    or str(log.name)[:2] == "._"):
                continue
            # meshCO_pointTriangleCO__IP_1e-2.out.txt
            scene_name, params = log.name.split("__")
            params = params.split(".")[0].split("_")
            if params[0] != "IP":
                continue
            constraint_solver, timestep = params
            timestep = float(timestep)
            if scene_name not in results:
                add_scene_to_results(results, scene_name, default_result)
            results[scene_name][timestep] = get_sim_status(log)
    save_results_csv(results)


if __name__ == "__main__":
    main()
