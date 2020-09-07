#!/bin/bash
R_CMD="-e criblejurkat::analysis(\"/root/${1}\")"
docker run -it --volume $(pwd):/root/ lbmc/criblejurkat:1.0.0 R $R_CMD
BASH_CMD="-e system(\"chown -R $(id -u ${USER}) /root/\")"
docker run -it --volume $(pwd):/root/ lbmc/criblejurkat:1.0.0 R $BASH_CMD
