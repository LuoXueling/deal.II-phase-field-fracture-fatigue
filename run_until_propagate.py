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

parser = ArgumentParser()
parser.add_argument("--command", type=str, required=True)

args = parser.parse_args()
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


def tune_parameter(s: str, d: dict):
    s_tmp = s
    for name, value in d.items():
        s_tmp = re.sub(f"{name} = " + r"(.*?)" + "\n", f"{name} = {value}\n", s_tmp)
    return s_tmp


if __name__ == "__main__":
    last_life = -1

    # Get current parameters
    current_para_path = mpi_command[where_is_para].split(".." + os.path.sep)[-1]
    file_in = open(current_para_path, "r")
    s = file_in.read()
    file_in.close()

    # Write new parameters
    t = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    project_name = re.findall(r"Project name = (.*?)\n", s)[0]
    if "Output sub-directory" not in s:
        s = s.replace(
            "subsection Project\n  ",
            "subsection Project\n  set Output sub-directory = \n  ",
        )
    _next_para = {"Output sub-directory": f"{project_name}-{t}"}
    s_out = tune_parameter(s, _next_para)
    output_path = f"output/{project_name}-{t}"
    current_para_path = (
        current_para_path.split(".prm")[0].split("_T_")[0] + f"_T_{t}.prm"
    )
    file_out = open(current_para_path, "w")
    file_out.write(s_out)
    file_out.close()
    log(f"Parameter file: {current_para_path}")

    # Get the new command
    mpi_command[where_is_para] = ".." + os.path.sep + current_para_path
    log(f"Current command: {' '.join(mpi_command)}")

    log(f"Monitoring folder: {output_path}")

    # Execute the process
    try:
        proc = subprocess.Popen(
            mpi_command,
            cwd=os.path.join(os.path.dirname(os.path.realpath(__file__)), "build"),
            stdout=subprocess.DEVNULL,
            # stderr=subprocess.STDOUT,
        )
    except Exception as e:
        log(f"Error when executing {' '.join(mpi_command)}")
        raise e
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
                    f"log-results.txt updated. Current # of records: {len(res)}. Cycle: {list(res['Step-Out'])[-1]}. Maximum phi: {list(res['Max-phi'])[-1]}"
                )
                last_no_records = len(res)
        else:
            res = None
        # The process terminates itself
        if return_code is not None:
            log(f"Process completed with return code: {return_code}")
            if return_code == 1:
                log("Return code is 1. Something bad happened. Terminating the job.")
                sys.exit()
            else:
                break
        else:
            if res is not None and list(res["Max-phi"])[-1] > 0.95:
                log("phi reaches 0.95. Killing the process.")
                proc.kill()
                break
