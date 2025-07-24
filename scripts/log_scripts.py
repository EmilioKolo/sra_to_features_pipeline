#!/usr/bin/env python3

"""
Functions using for logging purposes.
"""

import logging
import time


# Define logging print level
# Options: DEBUG, INFO, WARNING, ERROR or CRITICAL
logging.basicConfig(level=logging.INFO)

def log_print(
        s:str,
        level:str|int='debug',
        log_file:str='logs.txt'
    ) -> None:
    """
    Function that logs and prints a message.
    """
    # Check by level
    if (level is int and level>=5) or str(level).lower()=='critical':
        logging.critical(s)
        lev = 'CRITICAL:'
    elif (level is int and level>=4) or str(level).lower()=='error':
        logging.error(s)
        lev = 'ERROR:'
    elif (level is int and level>=3) or str(level).lower()=='warn':
        logging.warning(s)
        lev = 'WARNING:'
    elif (level is int and level>=2) or str(level).lower()=='info':
        logging.info(s)
        lev = 'INFO:'
    elif (level is int and level>=1) or str(level).lower()=='debug':
        logging.debug(s)
        lev = 'DEBUG:'
    else:
        logging.warning(f'Log level {level} not recognised.')
        lev = 'UNRECOGNISED:'
    # Print the message
    print(lev, s)
    # Add message to log_file
    with open(log_file, 'a+') as f:
        f.write(lev+' '+s+'\n')
    return None


def log_code(
        l:str,
        log_file:str='logs_script.txt'
    ) -> None:
    """
    Function that logs and prints a line of code.
    """
    # Define time of execution
    t = time.strftime('%H:%M%p %Z on %b %d, %Y')
    # Add message to log_file
    with open(log_file, 'a+') as f:
        f.write('['+t+'] '+l+'\n')
    return None
