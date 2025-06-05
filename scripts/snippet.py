#!/usr/bin/env python3


import os
from datetime import datetime


def log_annotate(log:str, message:str, log_ext:str='.log') -> None:
    """
    Log a message with a timestamp to the specified log file.
    """
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log + log_ext, 'a') as f:
        f.write(f"{timestamp} - {message}\n")
    print(f"{timestamp} - {message}")
    return None

def remove_file(file_path:str) -> None:
    """
    Remove a file if it exists.
    """
    if os.path.exists(file_path):
        os.remove(file_path)
        print(f"Removed file: {file_path}")
    else:
        print(f"File not found, cannot remove: {file_path}")
    return None

def run_silent(cmd:str, log:str, log_ext:str='.log') -> None:
    """
    Wrap the command in a bash subshell, redirect stdout and stderr.
    """
    log_annotate(
        log+'_timestamp',
        f'Running command: {cmd}',
        log_ext=log_ext
        )
    if log=='':
        print('Running command (no output redirection):', cmd)
        os.system(cmd)
    else:
        full_cmd = f"bash -c \"{cmd}\" >> {log}{log_ext} 2>&1"
        print('Running command:', cmd)
        os.system(full_cmd)
    log_annotate(
        log+'_timestamp',
        f'Finished command',
        log_ext=log_ext
        )
    return None
