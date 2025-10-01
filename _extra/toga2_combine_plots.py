#!/usr/bin/env python3

from typing import (
    Any, Dict, List, Optional, TypeVar, Union
)

import click
import json
import os
import subprocess
import sys

HL_PATH_STUB: str = '/projects/hillerlab/genome/gbdb-HL/{}/TOGA2/vs_{}'
ANNOT_UTR: str = 'query_annotation.with_utrs.bed'
ANNOT_CDS: str = 'query_annotation.bed'
UCSC: str = 'ucsc_browser_files'
BIGBED: str = 'HLTOGAannotVs{}v2.bb'
RETRO: str = 'retro'

LOCATION: str = os.path.dirname(os.path.dirname(__file__))
PLOTTER_SCRIPT: str = os.path.join(
    LOCATION, 'src', 'rust', 'target', 'release', 'combine_plots'
)

CONTEXT_SETTINGS: Dict[str, Any] = {
    'help_option_names': [
        '-h', '-help', '--help'
    ],
    'ignore_unknown_options': True,
    'allow_extra_args': True,
    'max_content_width': 150
}

JsonDict = TypeVar('JsonDict', bound=Dict[str, Dict[str, Dict[str, Union[str, List[str]]]]])

@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument(
    'ref',
    type=str,
    metavar='REF'
)
@click.argument(
    'queries',
    type=str,
    metavar='QUERIES'
)
@click.argument(
    'transcript',
    type=str,
    metavar='TRANSCRIPT'
)
@click.option(
    '--output',
    '-o',
    type=click.Path(exists=False),
    metavar='PATH',
    default=os.path.join(os.getcwd(), 'toga2_plot'),
    show_default=True,
    help=(
        """A path to the output plot. Note that:\n
\t* .png extension will be added automatically to the plot;\n
\t* All the parent directories must exist before the script is run\n
"""
    )
)

def main(
    ref: str, queries: str, transcript: str, output: Optional[str]
) -> None:
    """
    Plots all the projections corresponding to the progenitor transcript 
    for a given reference and a list of queries.\n
    Arguments are:\n
    \t*REF is a reference assembly name in Hiller Lab notation, e.g. hg38;\n
    \t*QUERIES is a comma-separated list of queries, e.g. mm10,HLrhiFer5,HLmyoViv2;\n
    \t*TRANSCRIPT is a name of reference transcript to plot. 
    For each query, the code will fetch all projections corresponding to this 
    transcripts from the query_annotation(.with_utrs).bed.\n
    \n
    Important caveats:\n
    1) Under the hood, he script runs the public TOGA2 script `combine_plots.rs` 
    which accepts input as a JSON file. The temporary JSON file produced by this script 
    is automatically removed on successful execution. If it crashed, do not forget 
    to remove the temporary files manually.\n
    2) The output file is an SVG format figure. This is intended behavior, and as of now 
    other output formats are not supported.\n
    3) Loss status filtering is currently not supported.
    """
    queries: List[str] = [x for x in queries.split(',') if x]
    output_file: str = output.split(os.sep)[-1]
    output_dir: str = os.path.dirname(output)
    json_input: JsonDict = {output_file: {}}
    for query in queries:
        query_path: str = HL_PATH_STUB.format(ref, query)
        if not os.path.exists(query_path):
            click.echo('ERROR: Input directory %s does not exist' % query_path)
            sys.exit(1)
        utr_path: str = os.path.join(query_path, ANNOT_UTR)
        if os.path.exists(utr_path):
            bed_path: str = utr_path
        else:
            cds_path: str = os.path.join(query_path, ANNOT_CDS)
            if os.path.exists(cds_path):
                click.echo(
                    (
                        'WARNING: UTR-containing annotation file %s '
                        'was not found for directory %s; '
                        'using %s instead'
                    ) % (ANNOT_UTR, query_path, ANNOT_CDS)
                )
                bed_path: str = cds_path
            else:
                click.echo(
                    'ERROR: Neither annotation file (%s or %s ) were found at %s' % (
                        ANNOT_UTR, ANNOT_CDS, query_path
                    )
                )
                sys.exit(1)
        prefix: str = ref.capitalize() if ref in ('hg38', 'mm10') else ref
        bigbed: str = BIGBED.format(prefix)
        bigbed_path: str = os.path.join(query_path, UCSC, bigbed)
        if not os.path.exists(bigbed_path):
            click.echo(
                'ERROR: BigBed file %s is missing from %s' % (bigbed, query_path)
            )
            sys.exit(1)
        projections: List[str] = []
        with open(bed_path, 'r') as h:
            for line in h:
                data: List[str] = line.strip().split('\t')
                if not data or not data[0]:
                    continue
                proj: str = data[3]
                proj_comps: List[str] = proj.split('#')
                if proj_comps[-1] == RETRO:
                    continue
                tr: str = '#'.join(proj_comps[:-1])
                if tr == transcript:
                    projections.append(proj)
        if projections:
            json_input[output_file][query] = {
                "bed": bed_path,
                "bigbed": bigbed_path,
                "projections": projections
            }
        else:
            click.echo(
                (
                    'WARNING: No projections corresponding to '
                    'transcript %s found for query %s'
                ) % (transcript, query)
            )
    json_name: str = f'{os.urandom(5).hex()}.json'
    with open(json_name, 'w') as h:
        json.dump(json_input, h, indent=4)
    cmd: str = f'{PLOTTER_SCRIPT} {json_name} {output_dir}'
    subprocess.call(cmd, shell=True)
    os.remove(json_name)

if __name__ == '__main__':
    main()