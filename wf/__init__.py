"""
Trim MSA using ClipKIT
"""

from enum import Enum
from pathlib import Path
import subprocess
from typing import Optional

from latch import small_task, workflow
from latch.types import LatchFile

class TrimmingMode(Enum):
    gappy = "gappy"
    smart_gap = "smart-gap"
    kpi = "kpi"
    kpi_gappy = "kpi-gappy"
    kpi_smart_gap = "kpi-smart-gap"
    kpic = "kpic"
    kpic_gappy = "kpic-gappy"
    kpic_smart_gap = "kpic-smart-gap"

@small_task
def trim_aln_task0(
    aln_fasta: LatchFile,
    output_file_name: Optional[str],
    gap_threshold: float = 0.9,
    trimming_mode: TrimmingMode = TrimmingMode.smart_gap
    ) -> LatchFile:

    # Set output file name
    if not output_file_name:
        trimmed_aln_fasta = Path("trimmed_aln.fna").resolve()
    else:
        trimmed_aln_fasta = output_file_name

    # Set gap threshold is user did not provide a value
    if not gap_threshold:
        gap_threshold=0.9

    _clipkit_cmd = [
        "clipkit",
        aln_fasta.local_path,
        "-o",
        str(trimmed_aln_fasta),
        "-m",
        trimming_mode.value,
        "-g",
        str(gap_threshold)
    ]

    subprocess.run(_clipkit_cmd)

    return LatchFile(str(trimmed_aln_fasta), f"latch:///{trimmed_aln_fasta}")


@workflow
def clipkit(
    aln_fasta: LatchFile,
    output_file_name: Optional[str] = None,
    gap_threshold: float = 0.9,
    trimming_mode: TrimmingMode = TrimmingMode.smart_gap
    ) -> LatchFile:

    """
    ClipKIT
    ----
    # ClipKIT, the multiple sequence alignment trimming toolkit
    ## About
    ClipKIT is a fast and flexible alignment trimming
    tool that keeps phylogenetically informative sites
    and removes those that display characteristics poor
    phylogenetic signal.

    <br /><br />

    If you found clipkit useful, please cite *ClipKIT:
    a multiple sequence alignment trimming software for
    accurate phylogenomic inference*. Steenwyk et al. 2020,
    PLoS Biology. doi:
    [10.1371/journal.pbio.3001007](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007).

    <br /><br />

    ## Modes
    Herein, we describe the various trimming modes implemented in ClipKIT. If you are unsure which is appropriate for you, we recommend using the default smart-gap trimming mode.

    <br /><br />

    ClipKIT can be run with eight different modes, which are specified with the -m/–mode argument. Default: ‘smart-gap’
    <br />
    - smart-gap: dynamic determination of gaps threshold
    - gappy: trim all sites that are above a threshold of gappyness (default: 0.9)
    - kpic: keep only parismony informative and constant sites
    - kpic-smart-gap: a combination of kpic- and smart-gap-based trimming
    - kpic-gappy: a combination of kpic- and gappy-based trimming
    - kpi: keep only parsimony informative sites
    - kpi-smart-gap: a combination of kpi- and smart-gap-based trimming
    - kpi-gappy: a combination of kpi- and gappy-based trimming

    __metadata__:
        display_name: Trim multiple sequence Alignments with ClipKIT
        author: Jacob L. Steenwyk
            name: Jacob L. Steenwyk
            email: jlsteenwyk@gmail.com
            github: https://github.com/JLSteenwyk
        repository: https://github.com/JLSteenwyk/ClipKIT
        license:
            id: MIT

    Args:

        aln_fasta:
            Input multiple sequence alignment in FASTA format
            __metadata__:
                display_name: "Input multiple sequence alignment"
                appearance:
					comment: "Input multiple sequence alignment in FASTA format"

        output_file_name:
            Output trimmed multiple sequence alignment in FASTA format
            __metadata__:
                display_name: "Output trimmed multiple sequence alignment"
                appearance:
					comment: "Output trimmed multiple sequence alignment in FASTA format"

        gap_threshold:
            Specifies gaps threshold (default: 0.9). Ignored if smart-gap is used.
			__metadata__:
				display_name: "Gappyness threshold"
				appearance:
					comment: "Specifies gaps threshold (default: 0.9). Ignored if smart-gap is used."

        trimming_mode:
            Mode used for trimming. See "About" for more information.
            __metadata__:
                display_name: "Trimming mode"
                appearance:
                    comment: "- smart-gap: dynamic determination of gaps threshold
                        - gappy: trim all sites that are above a threshold of gappyness (default: 0.9)
                        - kpic: keep only parismony informative and constant sites
                        - kpic-smart-gap: a combination of kpic- and smart-gap-based trimming
                        - kpic-gappy: a combination of kpic- and gappy-based trimming
                        - kpi: keep only parsimony informative sites
                        - kpi-smart-gap: a combination of kpi- and smart-gap-based trimming
                        - kpi-gappy: a combination of kpi- and gappy-based trimming"
    """


    return trim_aln_task0(
        aln_fasta=aln_fasta,
        output_file_name=output_file_name,
        gap_threshold=gap_threshold,
        trimming_mode=trimming_mode
        )
