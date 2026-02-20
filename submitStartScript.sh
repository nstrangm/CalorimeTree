#! /bin/bash

# load alice env
eval `$(which alienv) -w /software/flo/alice/sw --no-refresh printenv O2Physics/latest-masterbuild-o2`

# python3 start.py --config RunConfigFlorian.yaml

python3 start.py --config RunConfigIsopp2023skimmed.yaml --nodelist=clmachine
