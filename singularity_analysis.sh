#!/bin/bash
R_CMD="-e criblejurkat::analysis(\"$(pwd)\")"
singularity exec docker://lbmc/criblejurkat:1.0.0 R $R_CMD
