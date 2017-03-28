void drawJet(int jet) {

  TFile *f = new TFile("mixed_new.root");
  TTree *t1 = (TTree*)f->Get("jets");

  double pt, eta, phi;
  t1->SetBranchAddress("relPt0",&pt);
  t1->SetBranchAddress("dEta0",&eta);
  t1->SetBranchAddress("dPhi0",&phi);

  TGraph* jetGraph = new TGraph(200);
  char buffer[256];
  for(int i=0; i!=200; ++i) {
    sprintf(buffer,"relPt%i",i);
    t1->SetBranchAddress(buffer,&pt);
    sprintf(buffer,"dEta%i",i);
    t1->SetBranchAddress(buffer,&eta);
    sprintf(buffer,"dPhi%i",i);
    t1->SetBranchAddress(buffer,&phi);
    t1->GetEntry(jet);
    if(pt < 1E-5) {
      cout << "Only " << i << " particles loaded" << std::endl; 
      break;
    }
    jetGraph->SetPoint(i,eta,phi);
  }

  TCanvas* cv = new TCanvas("cv","cv",600,600);
  TH2D* histo = new TH2D("histo","histo",10,-1,1,10,-1,1);
  histo->GetXaxis()->SetTitle("#Delta #eta"); 
  histo->GetYaxis()->SetTitle("#Delta #phi"); 
  histo->Draw();
  jetGraph->Draw("P");
  TEllipse* ell = new TEllipse(0,0,0.8,0.8);
  ell->SetFillStyle(0);
  ell->SetLineColor(kRed);
  ell->SetLineStyle(kDashed);
  ell->SetLineWidth(2);
  ell->Draw();
  cv->SetGridx();
  cv->SetGridy();
}
