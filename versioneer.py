#!/usr/bin/env python3

"""Version update routine"""


import click

@click.command()
def versioneer() -> None:
    pass
## update __version__.py

## fetch the recent update description from changelog.md

## replace the previous change description in README.md with a recent 

if __name__ == '__main__':
    versioneer()