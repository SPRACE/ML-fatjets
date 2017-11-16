#!/bin/bash

# Remove the old executable
if [ -e generateJetImage ];
then 
  echo rm -f generateJetImage
  rm -f generateJetImage
fi

# Check if PYTHIA8DIR is set
if [ -z ${PYTHIA8DIR+x} ];
then
  echo "PYTHIA8DIR is unset, please set it to the location of the PYTHIA8 directory";
  return -1
else
  echo "PYTHIA8DIR is set to '$PYTHIA8DIR'";
fi

# Check for ROOT
if hash root-config 2>/dev/null;
then
  echo -n "ROOT version: "
  root-config --version
else
  echo "ROOT not found, please install ROOT (https://root.cern.ch)"
  return -1
fi

# Check for FASTJET
if hash fastjet-config 2>/dev/null;
then
  echo -n "FASTJET version: "
  fastjet-config --version
else
  echo "FASTJET not found, please install FASTJET (http://fastjet.fr)"
  return -1
fi

echo "Compiling..."
g++ generateJetImage.cc ${PYTHIA8DIR}/lib/libpythia8.a -I${PYTHIA8DIR}/include \
`root-config --cflags` `root-config --libs` \
`fastjet-config --cxxflags` `fastjet-config --libs` \
-I. -I/Users/trtomei/work/HEP/Software/Delphes/Delphes-3.4.1 \
-L/Users/trtomei/work/HEP/Software/Delphes/Delphes-3.4.1 -lDelphesNoFastJet \
-lz -o generateJetImage
echo "Done"
