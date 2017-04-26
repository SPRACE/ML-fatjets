#!/bin/bash
g++ generateJetImage.cc ${PYTHIA8DIR}/lib/libpythia8.a -I${PYTHIA8DIR}/include \
`root-config --cflags` `root-config --libs` `${FASTJETDIR}/bin/fastjet-config --cxxflags` \
`${FASTJETDIR}/bin/fastjet-config --libs` -DGZIPSUPPORT -lz -o generateJetImage