#!/usr/bin/env python3


import logging
import os
from datetime import datetime


def are_paths_same(path1:str, path2:str) -> bool:
    """
    Returns True if path1 and path2 refer to the same file/folder.
    """
    norm1 = os.path.abspath(os.path.normpath(path1))
    norm2 = os.path.abspath(os.path.normpath(path2))
    return norm1 == norm2

def change_output_ownership(output_dir:str) -> None:
    """
    Changes the ownership of all files and directories in output_dir.
    Requires HOST_UID and HOST_GID environment variables to be set.
    """
    # Get the UID and GID from environment variables
    try:
        uid = int(os.environ["HOST_UID"])
        gid = int(os.environ["HOST_GID"])
    except KeyError as e:
        w = 'Permissions could not be changed.'
        w += f' Missing environment variable: {e}'
        raise RuntimeError(w)
    # Walk through the output directory
    for root, dirs, files in os.walk(output_dir):
        for name in dirs + files:
            path = os.path.join(root, name)
            try:
                os.chown(path, uid, gid)
            except PermissionError as e:
                logging.warning(f"Permission denied on {path}: {e}")
            except FileNotFoundError:
                continue
    return None

def get_value(valname:str, sourcefile:str='config.env') -> str:
    """
    Gets the value valname from sourcefile.

    Args:
        valname (str): Name of the value to be found.
        sourcefile (str): Path to the file where the value is located.
                          Defaults to 'variables'.
    """
    # Initialize value that is returned
    val:str = ''
    # Open sourcefile
    with open(sourcefile, 'r') as f_src:
        # Go through lines in sourcefile
        for line in f_src.readlines():
            # If valname is in line
            if line.startswith(valname):
                # Split line by '='
                l_line = line.split('=')
                # Get the value and remove spaces and endline
                val = l_line[1].rstrip('\n').strip(' ')
    # Check if value was found
    if val=='':
        print(f'Value {valname} not found in "{sourcefile}".')
    return val

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
    if log=='':
        print('Running command (no output redirection):', cmd)
        os.system(cmd)
    else:
        log_annotate(
            log+'_timestamp',
            f'Running command: {cmd}',
            log_ext=log_ext
            )
        full_cmd = f"bash -c \"{cmd}\" >> {log}{log_ext} 2>&1"
        print('Running command:', cmd)
        os.system(full_cmd)
        log_annotate(
            log+'_timestamp',
            f'Finished command',
            log_ext=log_ext
            )
    return None
