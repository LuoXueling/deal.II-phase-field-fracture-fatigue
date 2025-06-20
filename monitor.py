#!/usr/bin/env python3
import subprocess
import sys
import time
import re
from datetime import datetime
from typing import List
import math
import os
import numpy as np
import pandas as pd
from argparse import ArgumentParser

SMOOTH_WINDOW = 10

parser = ArgumentParser()
parser.add_argument("--n_base", type=int, required=True)
parser.add_argument("--command", type=str, required=True)

args = parser.parse_args()
n_base = args.n_base
command = args.command

mpi_command = command.split(" ")
where_is_para = [".prm" in i for i in mpi_command].index(True)
para_name = mpi_command[where_is_para].split("/")[-1].split(".prm")[0]
output_file = (
    f"monitor-log-{para_name}-{datetime.now().strftime('%Y-%m-%d-%H-%M-%S')}.txt"
)


def log(*args, **kwargs):
    l = f"[{datetime.now().strftime('%m-%d-%Y,%H:%M:%S')}]"
    print(l, *args, **kwargs)
    if output_file is not None:
        _kwargs = kwargs.copy()
        with open(output_file, "a") as file:
            _kwargs.update(dict(file=file))
            print(l, *args, **_kwargs)


def window_slope_variance(n_cycles, data: List[float], pos: int, window: int) -> float:
    slopes = []
    start = max(0, pos - window + 1)
    for i in range(start, pos + 1):
        if i > 0:
            slopes.append((data[i] - data[i - 1]) / (n_cycles[i] - n_cycles[i - 1]))
    if not slopes:
        return 0.0
    mean = sum(slopes) / len(slopes)
    variance = sum((d - mean) ** 2 for d in slopes) / len(slopes)
    return variance


def is_last_point_rough(n_cycles, data: List[float]) -> bool:
    base_data = data[-10:-2]
    base_x = n_cycles[-10:-2]
    slopes = []
    for i in range(1, len(base_data)):
        slopes.append((base_data[i] - base_data[i - 1]) / (base_x[i] - base_x[i - 1]))

    mean_val = sum(slopes) / len(slopes)
    variance = sum((x - mean_val) ** 2 for x in slopes) / len(slopes)

    z_score_1 = abs(
        (data[-2] - data[-3]) / (n_cycles[-2] - n_cycles[-3]) - mean_val
    ) / (math.sqrt(variance) + 1e-50)
    z_score_2 = abs(
        (data[-1] - data[-2]) / (n_cycles[-1] - n_cycles[-2]) - mean_val
    ) / (math.sqrt(variance) + 1e-50)
    if z_score_1 > 20 and z_score_2 > 20:
        log(f"Rough. Z_score 1 = {z_score_1}, Z_score 2 = {z_score_2}")
        return True
    else:
        return False


def has_last_point_spike(n_cycles, data: List[float]) -> bool:
    base_data = data[:-2]
    base_x = n_cycles[:-2]
    slopes = []
    for i in range(1, len(base_data)):
        base_data[i]
        base_data[i - 1]
        base_x[i]
        base_x[i - 1]
        slopes.append((base_data[i] - base_data[i - 1]) / (base_x[i] - base_x[i - 1]))

    mean_val = sum(slopes) / len(slopes)
    variance = sum((x - mean_val) ** 2 for x in slopes) / len(slopes)

    z_score_1 = abs(
        (data[-2] - data[-3]) / (n_cycles[-2] - n_cycles[-3]) - mean_val
    ) / (math.sqrt(variance) + 1e-50)
    z_score_2 = abs(
        (data[-1] - data[-2]) / (n_cycles[-1] - n_cycles[-2]) - mean_val
    ) / (math.sqrt(variance) + 1e-50)

    if z_score_1 > 20 and z_score_2 > 20 and z_score_2 > z_score_1:
        log(f"Spike. Z_score 1 = {z_score_1}, Z_score 2 = {z_score_2}")
        return True
    else:
        return False


def is_stationary_failed(n_cycles, data: List[float]) -> bool:
    check_window = SMOOTH_WINDOW

    window_data = data[-check_window:]
    window_x = n_cycles[-check_window:]
    if not all(
        [window_data[i] > window_data[i - 1] for i in range(1, len(window_data))]
    ):
        # Not monotonically increasing. Happens when correcting residual after multiple cycle jumps
        return False
    sum_x = sum_y = sum_xy = sum_xx = 0.0
    for x, y in zip(window_x, window_data):
        sum_x += x
        sum_y += y
        sum_xy += x * y
        sum_xx += x * x

    denominator = check_window * sum_xx - sum_x**2

    slope = abs((check_window * sum_xy - sum_x * sum_y) / (denominator + 1e-50))
    if slope < 1e-5 * data[-1] / n_cycles[-1]:
        log(f"Stationary. Slope = {slope}. Reference slope = {data[-1]/n_cycles[-1]}.")
        return True
    else:
        return False


def is_last_point_anomaly(n_cycles: List[float], data: List[float]) -> bool:
    start_at = (
        [x > 1e-9 for x in data].index(True) if any([x > 1e-9 for x in data]) else None
    )
    if start_at is None:
        return False
    _n_cycles = n_cycles[start_at:]
    _data = data[start_at:]
    if len(_data) < SMOOTH_WINDOW:
        return False
    return (
        is_last_point_rough(_n_cycles, _data)
        or has_last_point_spike(_n_cycles, _data)
        or is_stationary_failed(_n_cycles, _data)
    )


def tune_parameter(s: str, d: dict):
    s_tmp = s
    for name, value in d.items():
        s_tmp = re.sub(f"{name} = " + r"(.*?)" + "\n", f"{name} = {value}\n", s_tmp)
    return s_tmp


possible_n_params = {
    "KristensenCLA": {
        "KristensenCLA": {"Name": "KristensenCLA", "at": [2, 3], "fa": [1]}
    },
    "CojocaruCycleJump": {
        "Cojocaru": {"Name": "Cojocaru", "at": [3], "fa": [2]},
        "CojocaruCLA": {"Name": "CojocaruCLA", "at": [3], "fa": [2]},
    },
    "LiCycleJump": {
        "Li": {"Name": "Li", "at": [3], "fa": [2]},
        "LiCLA": {"Name": "LiCLA", "at": [3], "fa": [2]},
    },
    "JacconCycleJump": {"Jaccon": {"Name": "Jaccon", "at": [2], "fa": [0]}},
    "JonasCycleJump": {
        "Jonas": {
            "Name": "Jonas",
            "at": [6],
            "fa": [0],
        },
        "JonasCLA": {
            "Name": "JonasCLA",
            "at": [6],
            "fa": [1],
        },
    },
    "YangCycleJump": {
        "Yang": {
            "Name": "Yang",
            "at": [5],
            "fa": [0],
        }
    },
    "Testing": {
        "Testing": {
            "Name": "Testing",
            "at": [3],
            "fa": [2],
        }
    },
}

N_propagate = n_base

changable_parameters = {
    "KristensenCLA": {
        "Adaptive timestep parameters": {
            3: {
                "name": "n_jump",
                "range": [1, N_propagate],  # None should be N when max(phi) reaches 1
                "type": int,
                "ratio": 0.5,
            }
        }
    },
    "Cojocaru": {
        "Adaptive timestep parameters": {
            3: {
                "name": "max_jumps",
                "range": [1, 5 * N_propagate],
                "type": int,
                "ratio": 0.5,
            }
        },
        "Fatigue accumulation parameters": {
            2: {
                "name": "q_jump",
                "range": [
                    1,
                    5 * N_propagate,
                ],  # None should be several times of N when max(phi) reaches 1
                "type": int,
                "ratio": 0.5,
            }
        },
    },
    "CojocaruCLA": {
        "Adaptive timestep parameters": {
            3: {
                "name": "max_jumps",
                "range": [1, 5 * N_propagate],
                "type": int,
                "ratio": 0.5,
            }
        },
        "Fatigue accumulation parameters": {
            2: {
                "name": "q_jump",
                "range": [
                    1,
                    5 * N_propagate,
                ],  # None should be several times of N when max(phi) reaches 1
                "type": int,
                "ratio": 0.5,
            }
        },
    },
    "Li": {
        "Adaptive timestep parameters": {
            3: {
                "name": "max_jump",
                "range": [N_propagate // 2],  # Could be a super large number
                "type": int,
                "ratio": 1,
            },
        },
        "Fatigue accumulation parameters": {
            2: {
                "name": "chi_cr",
                "range": [0.01, 1],
                "type": float,
                "ratio": 0.5,
            }
        },
    },
    "LiCLA": {
        "Adaptive timestep parameters": {
            3: {
                "name": "max_jump",
                "range": [N_propagate // 2],  # Could be a super large number
                "type": int,
                "ratio": 1,
            },
        },
        "Fatigue accumulation parameters": {
            2: {
                "name": "chi_cr",
                "range": [0.01, 1],
                "type": float,
                "ratio": 0.5,
            }
        },
    },
    "Jonas": {
        "Adaptive timestep parameters": {
            # Use the same value for lambda2 and lambda3, following the paper
            5: {
                "name": "lambda2",
                "range": [0.1, 10],
                "type": float,
                "ratio": 0.5,
            },
            6: {
                "name": "lambda3",
                "range": [0.1, 10],
                "type": float,
                "ratio": 0.5,
            },
        },
    },
    "JonasCLA": {
        "Adaptive timestep parameters": {
            # Use the same value for lambda2 and lambda3, following the paper
            5: {
                "name": "lambda2",
                "range": [0.1, 10],
                "type": float,
                "ratio": 0.5,
            },
            6: {
                "name": "lambda3",
                "range": [0.1, 10],
                "type": float,
                "ratio": 0.5,
            },
        },
    },
    "Yang": {
        "Adaptive timestep parameters": {
            2: {
                "name": "epsilon",
                "range": [1e-4, 1],
                "type": float,
                "ratio": 0.5,
            },
            5: {
                "name": "max_jump",
                "range": [5 * N_propagate],  # Could be a super large number
                "type": int,
                "ratio": 1,
            },
        }
    },
    "Jaccon": {
        "Adaptive timestep parameters": {
            2: {
                "name": "n_jump",
                "range": [1, N_propagate],
                "type": int,
                "ratio": 0.5,
            }
        }
    },
    "Testing": {
        "Stage": {
            1: {
                "name": "stage_value",
                "range": [1, 100],
                "type": int,
                "ratio": 0.01,
            }
        },
        "Adaptive timestep parameters": {
            2: {
                "name": "no_name",
                "range": [1e-5, 1],
                "type": float,
                "ratio": 0.5,
            },
            3: {
                "name": "no_name",
                "range": [1, 9999],
                "type": int,
                "ratio": 0.5,
            },
        },
    },
}


def which_scheme(s: str):
    at = re.findall(f"Adaptive timestep = " + r"(.*?)" + "\n", s)[0]
    fa = re.findall(f"Fatigue accumulation = " + r"(.*?)" + "\n", s)[0]
    at_param: str = re.findall(f"Adaptive timestep parameters = " + r"(.*?)" + "\n", s)[
        0
    ].strip()
    try:
        fa_param: str = re.findall(
            f"Fatigue accumulation parameters = " + r"(.*?)" + "\n", s
        )[0].strip()
    except:
        fa_param = None
    try:
        n_at = possible_n_params[at][fa]["at"]
        n_fa = possible_n_params[at][fa]["fa"]
        name = possible_n_params[at][fa]["Name"]
        if (
            len(at_param.split(" ")) not in n_at
            or (len(fa_param.split(" ")) if fa_param is not None else 0) not in n_fa
        ):
            raise Exception(
                f"Incorrect configuration for {name} ({at_param} and {fa_param}. "
                + f"Should be with length {n_at} and {n_fa}"
            )
        return name
    except KeyError:
        raise Exception(f"Unknown acceleration scheme ({at} and {fa})")


def get_param(s: str, loc: int, t):
    return t(s.split(" ")[loc])


def set_param(s: str, loc: int, new_val, t):
    seq = s.split(" ")
    new_val = float(new_val) if t == float else int(np.floor(new_val))
    seq[loc] = str(new_val)
    return " ".join(seq), new_val


if __name__ == "__main__":
    last_life = -1
    next_para = None
    job_start_time = time.time()

    ## Loop start from here
    while True:
        # Get current parameters
        current_para_path = mpi_command[where_is_para].split(".." + os.path.sep)[-1]
        file_in = open(current_para_path, "r")
        s = file_in.read()
        file_in.close()

        # What scheme are we using
        scheme = which_scheme(s)
        log(f"Using acceleration scheme: {scheme}")

        # Max number of cycles in the parameter file
        max_no_steps = int(
            re.findall(f"Max No of timesteps = " + r"(.*?)" + "\n", s)[0]
        )
        size1 = float(re.findall(f"Timestep size = " + r"(.*?)" + "\n", s)[0])
        size2 = float(
            re.findall(f"Timestep size to switch to = " + r"(.*?)" + "\n", s)[0]
        )
        n_size1 = float(
            re.findall(f"Switch timestep after steps = " + r"(.*?)" + "\n", s)[0]
        )
        max_cycles = int(np.floor(n_size1 * size1 + (max_no_steps - n_size1) * size2))
        log(f"The maximum number of cycles set in the paremeter file: {max_cycles}")

        # Initial parameters
        if next_para is None:
            next_para = {}
            for name in changable_parameters[scheme].keys():
                initial_params = re.findall(f"{name} = " + r"(.*?)" + "\n", s)[
                    0
                ].strip()
                for loc, info in changable_parameters[scheme][name].items():
                    initial_params, _ = set_param(
                        initial_params,
                        loc - 1,
                        info["range"][-1],
                        info["type"],
                    )
                next_para[name] = initial_params

        # Write new parameters
        t = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        project_name = re.findall(r"Project name = (.*?)\n", s)[0]
        _next_para = next_para.copy()
        _next_para.update({"Output sub-directory": f"{project_name}-{t}"})
        if "Output sub-directory" not in s:
            s = s.replace(
                "subsection Project\n  ",
                "subsection Project\n  set Output sub-directory = \n  ",
            )
        s_out = tune_parameter(s, _next_para)
        output_path = f"output/{project_name}-{t}"
        current_para_path = (
            current_para_path.split(".prm")[0].split("_T_")[0] + f"_T_{t}.prm"
        )
        file_out = open(current_para_path, "w")
        file_out.write(s_out)
        file_out.close()
        log(f"Parameter file: {current_para_path}")

        log(f"Monitoring folder: {output_path}")

        log(f"Parameters: {next_para}")

        # Get the new command
        mpi_command[where_is_para] = ".." + os.path.sep + current_para_path
        log(f"Current command: {' '.join(mpi_command)}")

        # Execute the process
        try:
            proc = subprocess.Popen(
                mpi_command,
                cwd=os.path.join(os.path.dirname(os.path.realpath(__file__)), "build"),
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
        except Exception as e:
            log(f"Error when executing {' '.join(mpi_command)}")
            raise e
        start_time = time.time()
        last_no_records = 0

        while True:
            return_code = proc.poll()
            # Monitor results
            if os.path.exists(os.path.join(output_path, "log-results.txt")):
                while True:
                    try:
                        res = pd.read_csv(
                            os.path.join(output_path, "log-results.txt"), sep=r"\s+"
                        )
                        break
                    except:
                        pass
                if len(res) > last_no_records:
                    log(
                        f"log-results.txt updated. Current # of records: {len(res)}. Cycle: {list(res['Step-Out'])[-1]}. Crack length: {list(res['Crack-length'])[-1]}"
                    )
                    last_no_records = len(res)
            else:
                res = None
            # The process terminates itself
            if return_code is not None:
                log(f"Process completed with return code: {return_code}")
                if return_code == 1:
                    log(
                        "Return code is 1. Something bad happened. Terminating the job."
                    )
                    sys.exit()
                else:
                    break
            # Anomaly detected and the process is terminated
            else:
                if res is not None:
                    res.drop_duplicates(subset="Step-Out", keep="last", inplace=True)
                    if is_last_point_anomaly(
                        list(res["Step-Out"]), list(res["Crack-length"])
                    ):
                        log("Crack length anomaly detected. Killing the process.")
                        proc.kill()
                        break
            # The job is gonna be terminated by HPC
            if (time.time() - job_start_time) / 60 / 60 / 24 > 6.95:
                log("Reaching 7 days. Terminating the job.")
                proc.kill()
                sys.exit(0)
            time.sleep(1)
        end_time = time.time()
        time_spent = end_time - start_time
        log(f"Time spent: {time_spent} s ({time_spent/60} min, {time_spent/3600} h)")
        if time_spent < 10:
            log(
                "Time spent < 10 s. Check if the setting is correct. Terminate the job."
            )
            break

        # Get current fatigue life and evaluate criteria
        life = list(
            sorted(
                [
                    int(i.split("solution_")[-1].split(".pvtu")[0])
                    for i in os.listdir(output_path)
                    if "pvtu" in i
                ]
            )
        )[-1]
        if last_life == -1:
            log(f"First attempt finished. Current life: {life}.")
        elif last_life >= max_cycles or life >= max_cycles:
            log(
                f"Last fatigue life {last_life} or current fatigue life {life} is greater than the maximum number of cycles set. Skipping termination criteria."
            )
        else:
            change_ratio = abs(life - last_life) / last_life
            log(
                f"Current life: {life}, last life: {last_life}, change ratio: {change_ratio}."
            )
            if change_ratio < 0.05:
                log(
                    f"Change ratio smaller than 5% ({change_ratio}). Terminating the job."
                )
                break
        last_life = life
        ############### decide the next parameters and check if it's legal ##############

        for name in changable_parameters[scheme].keys():
            initial_params = next_para[name]
            for loc, info in changable_parameters[scheme][name].items():
                old_val = get_param(initial_params, loc - 1, info["type"])
                initial_params, new_val = set_param(
                    initial_params,
                    loc - 1,
                    old_val * info["ratio"],
                    info["type"],
                )
                if new_val < info["range"][0]:
                    log(
                        f"Parameter {info['name']} reached the minimum allowed ({old_val}->{new_val} < {info['range'][0]}). Terminate the job."
                    )
                    sys.exit(0)
            next_para[name] = initial_params
