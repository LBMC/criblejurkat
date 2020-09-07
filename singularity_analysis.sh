#!/bin/bash
R_CMD="-e criblejurkat::analysis(\"${1}\")"
singularity exec docker://lbmc/criblejurkat:1.0.0 R $R_CMD
