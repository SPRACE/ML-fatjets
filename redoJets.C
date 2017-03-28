struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& struct1, const TLorentzVector& struct2)
  {
    return (struct1.Perp() > struct2.Perp());
  }
};

void redoJets() {
  std::ifstream infile;
  infile.open("mixed.txt");
  
  std::ofstream outfile;
  outfile.open("mixed_new.txt");

  std::string sep = "\t";
  std::string line;
  int isSignal;
  double jetPt;
  std::vector<TLorentzVector> theParticles;
  theParticles.reserve(200);
  int nEvents = 0;
  while (std::getline(infile, line)) {
      if(nEvents%500==0) std::cout << "nEvents = " << nEvents << std::endl;
      nEvents++;
      std::istringstream iss(line);
      iss >> isSignal >> jetPt;

      TLorentzVector part;
      TLorentzVector jet;
      double px,py,pz;
      theParticles.clear();
      for(int i=0; i!=200; ++i) {
	iss >> px >> py >> pz;
	part.SetXYZM(px*jetPt,py*jetPt,pz*jetPt,0);
	theParticles.push_back(part);
	jet += part;
      }
      
      std::sort(theParticles.begin(),theParticles.end(),greater_than_pt());  
      
      outfile << isSignal << sep;
      outfile << jet.Perp() << sep;
      outfile << jet.Eta() << sep;
      outfile << jet.Phi() << sep;
      
      double relPt,deltaEta,deltaPhi;
      for(int i=0; i!=200; ++i) {
	TLorentzVector& thisPart = theParticles.at(i);
	relPt = thisPart.Perp()/jet.Perp();
	deltaEta = thisPart.Eta() - jet.Eta();
	deltaPhi = ROOT::Math::VectorUtil::DeltaPhi(thisPart,jet);
	outfile << relPt << sep;
	outfile << deltaEta << sep;
	outfile << deltaPhi << sep;
      }
      outfile << std::endl;

  }

}
