#!/bin/bash
set -xev

[[ ! -f snap_024 ]] && wget 'https://docs.google.com/uc?export=download&id=1b4u1HLaXyBIXLfnqOYhZsA_iAWme3KWA' -O snap_024
[[ ! -f sub_060.0 ]] && wget 'https://docs.google.com/uc?export=download&id=1BDjiLpDpj5P-jwQrdSWWALuwhSpE5aVy' -O sub_060.0
python test_ci.py
