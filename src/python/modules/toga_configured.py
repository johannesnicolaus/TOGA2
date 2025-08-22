#!/usr/bin/env python3

"""
Wrapper over TogaMain for file-configured runs
"""

from modules.constants import TOGA2_SLOT2ARG
from modules.shared import CommandLineManager
from toga_main import __version__
from typing import Dict, List, Optional, Tuple, Union

import click
import os

__author__ = 'Yury V. Malovichko'
__year__ = '2024'
__all__ = (None)

LOCATION: str = os.path.dirname(os.path.abspath(__file__))
TOGA2: str = os.path.join(LOCATION, 'toga2.py')
COL_NUM: int = 2
REF: str = 'ref_2bit'
QUERY: str = 'query_2bit'
CHAIN: str = 'chain_file'
ANNOT: str = 'ref_annotation'
MANDATORY_ARGS: Tuple[str] = (REF, QUERY, CHAIN, ANNOT)
TRUE: str = 'True'
FALSE: str = 'False'
NONE: str = 'None'
PROJ_ARG_FILE: str = 'project_args.json'
VERSION: str = 'version'
SUPPORTED_ARGS: List[str] = list(TOGA2_SLOT2ARG.values())

class Toga2ConfiguredLauncher(CommandLineManager):
    __slots__ = ()

    def __init__(
        self, 
        config_file: click.File, 
        override: Optional[Union[str, None]]
    ) -> List[str]:
        self.v: bool = True
        args: Dict[str, str] = {}
        options: Dict[str, str] = {}
        for i, line in enumerate(config_file):
            data: List[str] = line.rstrip().split('\t')
            if not len(data):
                continue
            if len(data) != 2:
                self._die(
                    'Improper config file configuration at line %i; expected 2 columns, got %i' % (
                        i, len(data)
                    )
                )
            arg, value = data
            if arg == VERSION:
                self._echo('Specified TOGA2 version is %s' % arg)
                if value != __version__:
                    self._echo(
                        (
                            'Current TOGA2 version is %s but the configuration file was written by/for '
                            'version %s. Certain features listed in the provided file might differ or be '
                            'inaccessible in the current version'
                        ) % (__version__, value),
                        'warning'
                    )
                if arg not in TOGA2_SLOT2ARG.values():
                    self._echo(
                        'WARNING: Argument/option %s is not recognized; skipping' % arg,
                        'warning'
                    )
                    continue
                if arg in MANDATORY_ARGS:
                    args[arg] = value
                    continue
                if value == FALSE:
                    continue
                if value == NONE:
                    continue
                if value == TRUE:
                    # options += f' --{arg}'
                    options.append(f'--{arg}')
                    options[f'--{arg}'] = True
                options[f'--{arg}'] = value
        if override is not None:
            override = override.lstrip('"\'').rstrip('"\'')
            overriden_args: List[str] = override.split(' ')
            over_arg_num: int = len(overriden_args)
            curr: int = 0
            while curr < len(overriden_args):
                curr_arg: str = overriden_args[curr]
                curr_arg_name: str = curr_arg.lstrip('-')
                if curr_arg_name not in SUPPORTED_ARGS:
                    self._echo(
                        'Option %s is not supported by TOGA2 %s' % (curr_arg_name, __version__)
                    )
                    recognized_arg: bool = False
                else:
                    recognized_arg: bool = True
                if curr == over_arg_num - 1:
                    if recognized_arg:
                        options[curr_arg] = True
                    break
                next_arg: str = overriden_args[curr + 1]
                if next_arg.startswith('-'):
                    if recognized_arg:
                        options[curr_arg] = True
                    curr += 1
                    continue
                if recognized_arg:
                    options[curr_arg] = next_arg
                curr += 2
                continue
        cmd_args: List[str] = [args[REF], args[QUERY], args[CHAIN], args[ANNOT]]
        for option, value in options:
            if value is False or value is None:
                continue
            if value is True:
                cmd_args.append(option)
            else:
                cmd_args.extend((option, value))
        return cmd_args