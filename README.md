# README

## Overview

Code for studies on jet imaging and identification with Machine Learning.

## Requisites
   * ROOT6
   * Pythia8
   * FastJet 3.2.1

## Files

   * compileGenerateJetImage.sh: checks for the presence of ROOT, Pythia and Fastjet, and compiles the generateJetImage.cc code.
   * generateJetImage.cc: generates jet images according to the switches:

USAGE: generateJetImage [options]

Options:
 --outFileName   Output file name
 --numEvents, -n Number of tried events
 --averagePU     Average pileup
 --sqrts         Center of mass energy in GeV
 --minJetMass    Min jet mass in GeV
 --maxJetMass    Max jet mass in GeV
 --minJetPt      Min jet mass in GeV
 --maxJetPt      Max jet mass in GeV
 --process       Process: 0=bkg, 1=sig

## TODO

(everything after here is old and may be deprecated)

   * Try for a smart way to generate the LHE file for the signal (background can come from Pythia itself)
   * Try to make this compile in some other places other than in Thiago's machine
   * Get rid of the hardcoded names
   * Get rid of the stupid multiple loop drawJet.C
