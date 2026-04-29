#! /bin/bash

# load alice env
# eval `$(which alienv) -w /software/flo/alice/sw --no-refresh printenv O2Physics/latest-masterbuild-o2`
# use instead  source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.34.04/x86_64-ubuntu22.04-gcc114-opt/bin/thisroot.sh  
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.34.04/x86_64-ubuntu22.04-gcc114-opt/bin/thisroot.sh  
# python3 start.py --config RunConfigFlorian.yaml

# python3 start.py --config RunConfigIsopp2023skimmed.yaml
python3 start.py --config RunConfigIsopp2024skimmed.yaml
