#!/bin/sh
docker pull lbmc/criblejurkat:1.0.0
docker build ./ -t 'lbmc/criblejurkat:1.0.0'
docker push lbmc/criblejurkat:1.0.0
