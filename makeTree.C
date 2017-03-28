TTree* makeTree() {

  TTree* tree;
  std::string isSignal = "isSignal";
  std::string jetPt = "jetPt";
  std::string jetEta = "jetEta";
  std::string jetPhi = "jetPhi";
  std::string delim = ":";
  std::string varType = "/D";
  std::string branchDescriptor =
    isSignal + varType + delim +
    jetPt + delim +
    jetEta + delim +
    jetPhi + delim;

  char buffer[256];
  //std::cout << branchDescriptor << std::endl;
  for(int i=0; i!=200; ++i) {
    if(i!=199)
      sprintf (buffer,"relPt%i%sdEta%i%sdPhi%i%s",i,delim.c_str(),i,delim.c_str(),i,delim.c_str());
    else
      sprintf (buffer,"relPt%i%sdEta%i%sdPhi%i",i,delim.c_str(),i,delim.c_str(),i);
    string tempString(buffer);
    branchDescriptor = branchDescriptor + tempString;
  }
  std::cout << branchDescriptor << std::endl;
  tree = new TTree("jets","Jets");
  tree->ReadFile("mixed_new.txt",branchDescriptor.c_str());
  return tree;
}
