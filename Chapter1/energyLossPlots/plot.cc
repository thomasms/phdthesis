const int NumberOfBins = 5000;
const double lowEnergy = 1e5;
const double highEnergy = 1e12;

/// CONSTANTS ///
const double Na 	= 6.022e23; //Avagadros number
const double PI 	= 3.141;
const double RadiusE 	= 2.81794032528e-13; //classical electron radius in cm
//const double MassE	= 9.10938e-34; //mass of electron in grams
const double MassE	= 0.510998918e6; //mass of electron in eV
const double C		= 3.0e10; //speed of light in cm/s

struct material{
//properties of material
  int Z;
  double A;
  double IonEne;
};

struct particle{
  double z;
  double TMax;
  double mass;
};

double get_dEdx(const double energy,const particle &part,const material &mat, const double delta=0.);
double getBeta(const double energy,const double mass);
double getLorentzFactor(const double beta);
double getTMax(const double energy,const double mass);

void plot(){

  //define Argon
  material argon;
  argon.Z = 18;
  argon.A = 39.948;
  argon.IonEne = 11.6;

  //define particle
  particle muon;
  muon.z = 1;
  muon.mass = 105.7e6; // in eV

  particle proton;
  proton.z =1;
  proton.mass = 0.937e9;

  particle alpha;
  alpha.z = 2;
  alpha.mass = 3.727379e9;

  particle pion;
  pion.z = 1;
  pion.mass = 139.5701835e6;

  particle electron;
  electron.z = 1;
  electron.mass = 510.998e3;

  double energy[NumberOfBins];
  double dEdxMuon[NumberOfBins];
  double dEdxProton[NumberOfBins];
  double dEdxAlpha[NumberOfBins];
  double dEdxPion[NumberOfBins];
  double dEdxElectron[NumberOfBins];

  //loop over energy range
  for(int bin=0;bin<NumberOfBins;bin++){

	//plot in log x axis
    	double binWidth = (TMath::Log10(highEnergy) - TMath::Log10(lowEnergy))/double(NumberOfBins);

	double energyBinLog = TMath::Log10(lowEnergy) + binWidth*bin;
	double energyBin = TMath::Power(10,energyBinLog);

	//muon
	double mass = muon.mass;
	double beta = getBeta(energyBin,mass);
	double gamma = getLorentzFactor(beta);
	double Tmax = 2*MassE*beta*beta*gamma*gamma/(1. + (2*gamma*MassE/mass) + TMath::Power((MassE/mass),2) );
  	//get the Tmax for a muon
	muon.TMax = Tmax;

	//Proton
	mass = proton.mass;
	beta = getBeta(energyBin,mass);
	gamma = getLorentzFactor(beta);
	Tmax = 2*MassE*beta*beta*gamma*gamma/(1. + (2*gamma*MassE/mass) + TMath::Power((MassE/mass),2) );
  	//get the Tmax for a muon
	proton.TMax = Tmax;

	//Alpha
	mass = alpha.mass;
	beta = getBeta(energyBin,mass);
	gamma = getLorentzFactor(beta);
	Tmax = 2*MassE*beta*beta*gamma*gamma/(1. + (2*gamma*MassE/mass) + TMath::Power((MassE/mass),2) );
  	//get the Tmax for a muon
	alpha.TMax = Tmax;

	//Pion
	mass = pion.mass;
	beta = getBeta(energyBin,mass);
	gamma = getLorentzFactor(beta);
	Tmax = 2*MassE*beta*beta*gamma*gamma/(1. + (2*gamma*MassE/mass) + TMath::Power((MassE/mass),2) );
  	//get the Tmax for a muon
	pion.TMax = Tmax;

	//Electron
	mass = electron.mass;
	beta = getBeta(energyBin,mass);
	gamma = getLorentzFactor(beta);
	Tmax = 2*MassE*beta*beta*gamma*gamma/(1. + (2*gamma*MassE/mass) + TMath::Power((MassE/mass),2) );
  	//get the Tmax for a muon
	electron.TMax = Tmax;

	//get dE/dx
	double dEdxMuonBin = get_dEdx(energyBin,muon,argon);
	double dEdxProtonBin = get_dEdx(energyBin,proton,argon);
	double dEdxAlphaBin = get_dEdx(energyBin,alpha,argon);
	double dEdxPionBin = get_dEdx(energyBin,pion,argon);
	double dEdxElectronBin = get_dEdx(energyBin,electron,argon);

	//fill the arrays
	energy[bin] = energyBin/1e9; 		//convert to GeV

	dEdxMuon[bin] = dEdxMuonBin/1e6; 	//convert to MeV 
	dEdxProton[bin] = dEdxProtonBin/1e6;
	dEdxAlpha[bin] = dEdxAlphaBin/1e6;
	dEdxPion[bin] = dEdxPionBin/1e6;
	dEdxElectron[bin] = dEdxElectronBin/1e6;

//	std::cout<< "\nEnergy: [eV]" << energyBin << ", mass [eV]: " << mass << ", beta: " << beta << ", gamma: " <<gamma << ", dE/dxMuon =" <<dEdxMuonBin << ", Tmax= " << muon.TMax;
  }

 
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogx();
  gPad->SetLogy();
  //make a graph of these
  TGraph * grMuon = new TGraph(NumberOfBins,energy,dEdxMuon);
  TGraph * grProton = new TGraph(NumberOfBins,energy,dEdxProton);
  TGraph * grAlpha = new TGraph(NumberOfBins,energy,dEdxAlpha);
  TGraph * grPion = new TGraph(NumberOfBins,energy,dEdxPion);
  TGraph * grElectron = new TGraph(NumberOfBins,energy,dEdxElectron);
  grMuon->SetLineColor(kAzure+7);
  grMuon->SetFillColor(kWhite);
  grMuon->SetMarkerColor(kAzure+7);
  grMuon->SetTitle("Muon");
  grMuon->Draw("AC");
  grMuon->GetXaxis()->SetTitle("Energy [GeV]");
  grMuon->GetYaxis()->SetTitle("-dE/dx [MeVcm^{2}g^{-1}]");
  grProton->SetLineColor(kRed);
  grProton->SetMarkerColor(kRed);
  grProton->SetTitle("Proton");
  grProton->Draw("C");
  grAlpha->SetLineColor(kBlack);
  grAlpha->SetMarkerColor(kBlack);
  grAlpha->SetTitle("Alpha");
  grAlpha->Draw("C");
  grPion->SetLineColor(kTeal+3);
  grPion->SetMarkerColor(kTeal+3);
  grPion->SetTitle("Pion");
  grPion->Draw("C");
  grElectron->SetLineColor(kBlue);
  grElectron->SetMarkerColor(kBlue);
  grElectron->SetTitle("Electron");
  grElectron->Draw("C");

}

double get_dEdx(const double energy,const particle &part,const material &mat, const double delta){
  
  double result = 0.;

  double A = mat.A;
  int Z = mat.Z;
  double IonEne = mat.IonEne;

  double z = part.z;
  double Tmax = part.TMax;
  double mass = part.mass;

  double beta = getBeta(energy,mass);
  double gamma = getLorentzFactor(beta);

//  std::cout << "\n\tTEST: Tmax: " << Tmax << ", mass:" << mass;

  if(2*MassE*beta*beta*gamma*gamma*Tmax/(IonEne*IonEne) <= 0)return result;
  if(beta<=0)return result;

  result = 4*PI*Na*RadiusE*RadiusE*MassE*z*z*(Z/A)*(1/(beta*beta))*(0.5*TMath::Log(2*MassE*beta*beta*gamma*gamma*Tmax/(IonEne*IonEne)) - (beta*beta) - (delta/2.));

  if(result>1e11 || result <1e-3)result =0.;

  return result;

}

double getBeta(const double energy,const double mass){

  double result = 0.;

  if(energy<mass)return result;

  result = TMath::Sqrt(1 - (mass*mass/(energy*energy)) );

  return result;
}

double getLorentzFactor(const double beta){

  double result = 0.;

  double inverseResult = TMath::Sqrt(1-(beta*beta));
  if(inverseResult <= 0)return 1.;

  result = 1/(TMath::Sqrt(1-(beta*beta)));
  return result;

}
/*
double getTMax(const double energy,const double mass){
 
  double result = 0.;

  double beta = getBeta(energy,mass);
  double gamma = getLorentzFactor(beta);

  result = 2*MassE*beta*beta*gamma*gamma/(1 + (2*gamma*MassE/mass) + TMath::Power((MassE/mass),2) );
  return result;

}
*/
