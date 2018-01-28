//declare PI
const double PI = 3.141;

//delcare the size of the array
const Int_t n = 1000;

const double maxEnergy = 10.;

//baseline in km
const double L = 2300;
const double baseline = L;

  //oscillation parameters
const double s12 = 0.559;
const double s23 = 0.714;
const double s13 = 0.155;
const double c12 = 0.829;
const double c23 = 0.700;
const double c13 = 0.988;
    
//Mass squared parameters eV^2 (set for NH)
const double m21 = 7.59e-5;
const double m31 = 2.54e-3;
const double m32 = 2.46e-3;

//IH
const double m12 = m21;
const double m13 = -2.34e-3;
const double m23 = -2.42e-3;
  
//unknown parameter delta cp

const double delta0 = 0;
const double delta1 = PI/2.0;
const double delta2 = PI;
const double delta3 = 3*PI/2.0;
 
const double delta_cp = PI;
const double delta_cp_low = PI/2.0;
const double delta_cp_high = 3*PI/2.0;
  
//matter parameter
const double density = 2.8; //in g/cm^3
const double alpha = 7.56e-5*density;
  
//factor for sin and cos conversion
const double factor = 1.27;

double NuMuToNuEMatterProb(double energy,const double delta){

  double numuToNuEProb = 0.;

  numuToNuEProb = 0.002 + (4*c13*c13*s13*s13*s23*s23*(1+((2*(alpha*energy)/(m31))*(1-(2*s13*s13))))*sin(factor*m31*L/energy)*sin(factor*m31*L/energy))
        + (8*c13*c13*s12*s13*s23*(c12*c23*cos(delta) - s12*s13*s23)*cos(factor*m32*L/energy)*sin(factor*m31*L/energy)*sin(factor*m21*L/energy))
        - (8*c13*c13*c12*c23*s12*s13*s23*sin(delta)*sin(factor*m32*L/energy)*sin(factor*m31*L/energy)*sin(factor*m21*L/energy))
        + (4*s12*s12*c13*c13*(c12*c12*c23*c23 + s12*s12*s23*s23*s13*s13 - 2*c12*c23*s12*s23*s13*cos(delta))*sin(factor*m21*L/energy)*sin(factor*m21*L/energy))
        - (8*c13*c13*s13*s13*s23*s23*(alpha*energy)*(L/energy)*cos(factor*m32*L/energy)*sin(factor*m31*L/energy)*(1-(2*s13*s13))/4.0); 

  return numuToNuEProb;
}

double NuMuToNuEMatterProbIH(double energy,const double delta){

  double numuToNuEProb = 0.;

  numuToNuEProb = 0.007 + (4*c13*c13*s13*s13*s23*s23*(1+((2*(alpha*energy)/(m13))*(1-(2*s13*s13))))*sin(factor*m13*L/energy)*sin(factor*m13*L/energy))
        + (8*c13*c13*s12*s13*s23*(c12*c23*cos(delta) - s12*s13*s23)*cos(factor*m23*L/energy)*sin(factor*m13*L/energy)*sin(factor*m12*L/energy))
        - (8*c13*c13*c12*c23*s12*s13*s23*sin(delta)*sin(factor*m23*L/energy)*sin(factor*m13*L/energy)*sin(factor*m12*L/energy))
        + (4*s12*s12*c13*c13*(c12*c12*c23*c23 + s12*s12*s23*s23*s13*s13 - 2*c12*c23*s12*s23*s13*cos(delta))*sin(factor*m12*L/energy)*sin(factor*m12*L/energy))
        - (8*c13*c13*s13*s13*s23*s23*(alpha*energy)*(L/energy)*cos(factor*m23*L/energy)*sin(factor*m13*L/energy)*(1-(2*s13*s13))/4.0); 

  return numuToNuEProb;
}

double FindFirstMaximum(TGraph * graph, double Xlower,double Xupper){

   double maxX = 0.;
   double maxY = 0.;

   int size = graph->GetN();
   double* x = graph->GetX();

   //loop over x values to find values above xlower
   for(int i=0;i<size;i++){
	if(x[i] <Xlower)continue;
	if(x[i]>Xupper)continue;
	
	double evaluatedValue = graph->Eval(x[i]);

	if(maxY<=evaluatedValue){
	  maxY = evaluatedValue;
	  maxX = x[i];
	}
   }

   std::cout << "\nMaxX: " << maxX << ", maxY: " << maxY;

   return maxX;
}

void osc_prob()
{

  //x = L/E in km/GeV
  Double_t x[n];
  Double_t E[n];

  Double_t P_numu_numu0[n];
  Double_t P_numu_numu1[n];
  Double_t P_numu_numu2[n];
  Double_t P_numu_numu3[n];

  Double_t P_numu_nue[n];
  Double_t P_numu_nue_vac_NH[n];
  Double_t P_numu_nue_vac_IH[n];
  Double_t P_numu_nue_NH[n];
  Double_t P_numu_nue_low_NH[n];
  Double_t P_numu_nue_high_NH[n];
  Double_t P_numu_nue_IH[n];
  Double_t P_numu_nue_low_IH[n];
  Double_t P_numu_nue_high_IH[n];

  for(Int_t i =0;i<n;i++)
    {
      double binWidth = maxEnergy/( double(n));

      E[i] = binWidth*((2*i) + 1)/2.;
      x[i] = L/E[i];
      //x[i] = i;

      //
      //-----------------------numu to numu---------------------------
      //

      P_numu_numu0[i] = 1 - (4*s23*s23*c13*c13*sin(factor*m23*x[i])*sin(factor*m23*x[i])*(c12*c12*c23*c23 + s12*s12*s13*s13*s23*s23 - 2*c12*c23*s12*s13*s23*cos(delta0)))
	- (4*s23*s23*c13*c13*sin(factor*m13*x[i])*sin(factor*m13*x[i])*(s12*s12*c23*c23 + c12*c12*s13*s13*s23*s23 + 2*c12*c23*s12*s13*s23*cos(delta0)))
	- (4*c12*c12*c23*c23 + s12*s12*s13*s13*s23*s23 - 2*c12*c23*s12*s13*s23*cos(delta0))*(s12*s12*c23*c23 + c12*c12*s13*s13*s23*s23 + 2*c12*c23*s12*s13*s23*cos(delta0))*sin(factor*m12*x[i])*sin(factor*m12*x[i]);

      P_numu_numu1[i] = 1 - (4*s23*s23*c13*c13*sin(factor*m23*x[i])*sin(factor*m23*x[i])*(c12*c12*c23*c23 + s12*s12*s13*s13*s23*s23 - 2*c12*c23*s12*s13*s23*cos(delta1)))
	- (4*s23*s23*c13*c13*sin(factor*m13*x[i])*sin(factor*m13*x[i])*(s12*s12*c23*c23 + c12*c12*s13*s13*s23*s23 + 2*c12*c23*s12*s13*s23*cos(delta1)))
	- (4*c12*c12*c23*c23 + s12*s12*s13*s13*s23*s23 - 2*c12*c23*s12*s13*s23*cos(delta1))*(s12*s12*c23*c23 + c12*c12*s13*s13*s23*s23 + 2*c12*c23*s12*s13*s23*cos(delta1))*sin(factor*m12*x[i])*sin(factor*m12*x[i]);

      P_numu_numu2[i] = 1 - (4*s23*s23*c13*c13*sin(factor*m23*x[i])*sin(factor*m23*x[i])*(c12*c12*c23*c23 + s12*s12*s13*s13*s23*s23 - 2*c12*c23*s12*s13*s23*cos(delta2)))
	- (4*s23*s23*c13*c13*sin(factor*m13*x[i])*sin(factor*m13*x[i])*(s12*s12*c23*c23 + c12*c12*s13*s13*s23*s23 + 2*c12*c23*s12*s13*s23*cos(delta2)))
	- (4*c12*c12*c23*c23 + s12*s12*s13*s13*s23*s23 - 2*c12*c23*s12*s13*s23*cos(delta2))*(s12*s12*c23*c23 + c12*c12*s13*s13*s23*s23 + 2*c12*c23*s12*s13*s23*cos(delta2))*sin(factor*m12*x[i])*sin(factor*m12*x[i]);

      P_numu_numu3[i] = 1 - (4*s23*s23*c13*c13*sin(factor*m23*x[i])*sin(factor*m23*x[i])*(c12*c12*c23*c23 + s12*s12*s13*s13*s23*s23 - 2*c12*c23*s12*s13*s23*cos(delta3)))
	- (4*s23*s23*c13*c13*sin(factor*m13*x[i])*sin(factor*m13*x[i])*(s12*s12*c23*c23 + c12*c12*s13*s13*s23*s23 + 2*c12*c23*s12*s13*s23*cos(delta3)))
	- (4*c12*c12*c23*c23 + s12*s12*s13*s13*s23*s23 - 2*c12*c23*s12*s13*s23*cos(delta3))*(s12*s12*c23*c23 + c12*c12*s13*s13*s23*s23 + 2*c12*c23*s12*s13*s23*cos(delta3))*sin(factor*m12*x[i])*sin(factor*m12*x[i]);

      //
      //-----------------------numu to nume---------------------------
      //

      //vacuum oscillation
      P_numu_nue_vac_NH[i] = (4*c13*c13*s13*s13*s23*s23*sin(factor*m31*x[i])*sin(factor*m31*x[i]))
	+ (8*c13*c13*s12*s13*s23*(c12*c23*cos(delta_cp) - s12*s13*s23)*cos(factor*m32*x[i])*sin(factor*m31*x[i])*sin(factor*m21*x[i]))
	- (8*c13*c13*c12*c23*s12*s13*s23*sin(delta_cp)*sin(factor*m32*x[i])*sin(factor*m31*x[i])*sin(factor*m21*x[i]))
       	+ (4*s12*s12*c13*c13*(c12*c12*c23*c23 + s12*s12*s23*s23*s13*s13 - 2*c12*c23*s12*s23*s13*cos(delta_cp))*sin(factor*m21*x[i])*sin(factor*m21*x[i])); 

      P_numu_nue_vac_IH[i] = (4*c13*c13*s13*s13*s23*s23*sin(factor*m13*x[i])*sin(factor*m13*x[i]))
	+ (8*c13*c13*s12*s13*s23*(c12*c23*cos(delta_cp) - s12*s13*s23)*cos(factor*m32*x[i])*sin(factor*m13*x[i])*sin(factor*m21*x[i]))
	- (8*c13*c13*c12*c23*s12*s13*s23*sin(delta_cp)*sin(factor*m23*x[i])*sin(factor*m13*x[i])*sin(factor*m21*x[i]))
       	+ (4*s12*s12*c13*c13*(c12*c12*c23*c23 + s12*s12*s23*s23*s13*s13 - 2*c12*c23*s12*s23*s13*cos(delta_cp))*sin(factor*m21*x[i])*sin(factor*m21*x[i])); 

      //normal hierarchy
      P_numu_nue_NH[i] = NuMuToNuEMatterProb(E[i],delta_cp);
      //inverted hierarchy
      P_numu_nue_IH[i] = NuMuToNuEMatterProbIH(E[i],delta_cp);

      P_numu_nue_low_NH[i] = NuMuToNuEMatterProb(E[i],delta_cp_low);
      P_numu_nue_low_IH[i] = NuMuToNuEMatterProbIH(E[i],delta_cp_low);

      P_numu_nue_high_NH[i] = NuMuToNuEMatterProb(E[i],delta_cp_high);
      P_numu_nue_low_IH[i] = NuMuToNuEMatterProbIH(E[i],delta_cp_high);

    }

  TCanvas *c1 = new TCanvas("c1","c1");
  TCanvas *c2 = new TCanvas("c2","c2");
  TCanvas *c3 = new TCanvas("c3","c3");

  gr0 = new TGraph(n,E,P_numu_numu0);
  gr1 = new TGraph(n,E,P_numu_numu1);
  gr2 = new TGraph(n,E,P_numu_numu2);
  gr3 = new TGraph(n,E,P_numu_numu3);

  gr = new TGraph(n,x,P_numu_nue_NH);

  TGraphErrors * gr_NH = new TGraphErrors(n,E,P_numu_nue_NH,0,0);
  gr_low_NH = new TGraph(n,E,P_numu_nue_low_NH);
  gr_high_NH = new TGraph(n,E,P_numu_nue_high_NH);

  gr_vac_NH = new TGraph(n,E,P_numu_nue_vac_NH);
  gr_vac_IH = new TGraph(n,E,P_numu_nue_vac_IH);

  TGraphErrors *gr_IH = new TGraphErrors(n,E,P_numu_nue_IH,0,0);
  gr_low_IH = new TGraph(n,E,P_numu_nue_low_IH);
  gr_high_IH = new TGraph(n,E,P_numu_nue_high_IH);

  double maximumEnergy = FindFirstMaximum(gr_NH,2,maxEnergy);

  gr_NH->GetXaxis()->SetRangeUser(0,maxEnergy);
  gr_NH->GetYaxis()->SetRangeUser(0,0.2);
  gr_NH->SetLineWidth(3);
  gr_NH->SetFillColor(0);
  gr_low_NH->SetFillColor(0);
  gr_high_NH->SetFillColor(0);
  gr_low_NH->SetLineColor(1);
  gr_high_NH->SetLineColor(1);
  gr_IH->SetLineColor(4);
  gr_IH->SetLineWidth(2);
  gr_low_IH->SetLineColor(4);
  gr_high_IH->SetLineColor(4);
  gr_IH->SetFillColor(0);
  gr_low_IH->SetFillColor(0);
  gr_high_IH->SetFillColor(0);

  c1->cd();

  //titles
  gr_NH->SetTitle("NH #delta = 180");
  gr_low_NH->SetTitle("NH #delta = 90");
  gr_high_NH->SetTitle("NH #delta = 270");
  gr_IH->SetTitle("IH #delta = 180");
  gr_low_IH->SetTitle("IH #delta = 90");
  gr_high_IH->SetTitle("IH #delta = 270");

  gr_NH->GetXaxis()->SetTitle("Energy [GeV]");
  gr_NH->GetYaxis()->SetTitle("P(#nu_{#mu} to #nu_{e})");

  //Draw
  gr_NH->Draw("AC");
  gr_low_NH->Draw("C");
  gr_high_NH->Draw("C");
  gr_IH->Draw("C");
  gr_low_IH->Draw("C");
  gr_high_IH->Draw("C"); 

  c2->cd();

  //titles
  gr_vac_NH->SetTitle("NH #delta = 180");
  gr_vac_IH->SetTitle("IH #delta = 180");
  gr_vac_NH->GetXaxis()->SetTitle("Energy [GeV]");
  gr_vac_NH->GetYaxis()->SetTitle("P(#nu_{#mu} to #nu_{e})");


  gr_vac_NH->SetLineColor(2);
  gr_vac_IH->SetLineColor(3);

  gr_vac_NH->GetXaxis()->SetRangeUser(0,maxEnergy);
  gr_vac_NH->GetYaxis()->SetRangeUser(0,0.2);

  //Draw
  gr_vac_NH->Draw("AC");
  gr_vac_IH->Draw("C");

  c3->cd();
  gr0->GetXaxis()->SetTitle("Energy [GeV]");
  gr0->GetYaxis()->SetTitle("P(#nu_{#mu} to #nu_{#mu})");
  gr0->Draw("AC");
  //gr1->Draw("C");
  //gr2->Draw("C");
  //gr3->Draw("C");

}

