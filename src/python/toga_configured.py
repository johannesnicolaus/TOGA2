#!/usr/bin/env python3

"""
Wrapper over TogaMain for file-configured runs
"""

from .shared import CommandLineManager
from typing import Optional, Union

import click

class Toga2ConfiguredLauncher(CommandLineManager):
    __slots__ = ()

    def __init__(self, config_file: click.File, override: Optional[Union[str, None]]) -> None:
        pass