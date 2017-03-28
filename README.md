# README

## Overview

This is Thiago's initial code for studies in ML. Probably nothing in this repo works well, but it is good to document it while we work...

## Requisites
   * ROOT
   * Pythia8

## Files

   * compile.sh: script to compile the executables that generate the signal and background samples.
   * redoJets.C: fixes the output of these files from (isSignal, jetPt, vector(Px_i, Py_i, Pz_i)/jetPt) to (isSignal, jetPt, jetEta, jetPhi, vector(relPt_i, dEta_i, dPhi_i)), where relPt_i = pt_i/jetPt, dEta_i = eta_i - jetEta, dPhi_i = phi_i - jetPhi mod 2pi. "i" refers to the i-th particle. Run inside root.
   * makeTree.C: format shift the table from .txt to .root . Run inside root, and save the tree somewhere.
   * drawJet.C: make plot of a given jet. Run inside root for a given jet in the tree with .L drawJet.C; drawJet(nJet)

## TODO

   * Try for a smart way to generate the LHE file for the signal (background can come from Pythia itself)
   * Try to make this compile in some other places other than in Thiago's machine
   * Get rid of the hardcoded names
   * Get rid of the stupid multiple loop drawJet.C