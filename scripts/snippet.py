#!/usr/bin/env python3

import os

def run_silent(cmd:str, log:str) -> None:
    # Wrap the command in a bash subshell, redirect stdout and stderr
    full_cmd = f"bash -c \"{cmd}\" >> {log} 2>&1"
    os.system(full_cmd)
    return None