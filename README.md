# README

## Overview

Code for studies on jet imaging and identification with Machine Learning.

## Requisites
   * ROOT6
   * Pythia8
   * FastJet 3.2.1

## Files

   * compileGenerateJetImage.sh: checks for the presence of ROOT, Pythia and Fastjet, and compiles the generateJetImage.cc code.
   * generateJetImage.cc: generates jet images according to a series of switches. Run with ./generateJetImage --help to see the options.
   
## TODO

(everything after here is old and may be deprecated)

   * Try for a smart way to generate the LHE file for the signal (background can come from Pythia itself)
   * Try to make this compile in some other places other than in Thiago's machine
   * Get rid of the hardcoded names
   * Get rid of the stupid multiple loop drawJet.C
