#!/usr/bin/env python3


import os
from datetime import datetime


def run_silent(cmd:str, log:str, log_ext:str='.log') -> None:
    """
    Wrap the command in a bash subshell, redirect stdout and stderr.
    """
    if log=='':
        print('Running command (no output redirection):', cmd)
        os.system(cmd)
    else:
        full_cmd = f"bash -c \"{cmd}\" >> {log}{log_ext} 2>&1"
        print('Running command:', cmd)
        os.system(full_cmd)
    return None

def run_silent(cmd:str, log:str) -> None:
    # Wrap the command in a bash subshell, redirect stdout and stderr
    full_cmd = f"bash -c \"{cmd}\" >> {log} 2>&1"
    os.system(full_cmd)
    return None