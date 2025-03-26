#!/bin/bash
eval `$(which alienv) -w /software/flo/alice/sw --no-refresh printenv O2Physics/latest-masterbuild-o2`
export JALIEN_TOKEN_CERT=/alf/data/flo/token/tokencert_${UID}.pem
export JALIEN_TOKEN_KEY=/alf/data/flo/token/tokenkey_${UID}.pem
export ALIENPY_DEBUG_FILE=/tmp/${UID}_alien_py.log

o2-analysis-centrality-table -b --configuration json://configuration.json |
o2-analysis-event-selection -b --configuration json://configuration.json |
o2-analysis-multiplicity-table -b --configuration json://configuration.json |
o2-analysis-ft0-corrected-table -b --configuration json://configuration.json |
o2-analysis-timestamp -b --configuration json://configuration.json |
o2-analysis-trackselection -b --configuration json://configuration.json |
o2-analysis-track-propagation -b --configuration json://configuration.json |
o2-analysis-tracks-extra-v002-converter -b --configuration json://configuration.json |
o2-analysis-je-jet-deriveddata-producer -b --configuration json://configuration.json |
o2-analysis-je-jet-deriveddata-producer-dummy -b --configuration json://configuration.json |
o2-analysis-je-jet-finder-data-charged -b --configuration json://configuration.json |
o2-analysis-je-estimator-rho -b --configuration json://configuration.json |
o2-analysis-je-gamma-jet-tree-producer -b --configuration json://configuration.json |
o2-analysis-je-emcal-correction-task -b --configuration json://configuration.json --aod-file @input_data.txt --aod-writer-json OutputDirector.json
# o2-analysis-timestamp -b --configuration json://configuration.json \
# | o2-analysis-trackselection -b --configuration json://configuration.json \
# | o2-analysis-track-propagation -b --configuration json://configuration.json \
# | o2-analysis-event-selection -b --configuration json://configuration.json \
# | o2-analysis-je-emc-eventselection-qa -b --configuration json://configuration.json \
# | o2-analysis-je-emcal-correction-task -b --configuration json://configuration.json --aod-file @input_data.txt