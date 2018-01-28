{
  int size =6;

  //read in files
  string filename_material1_tot = "neutronHe3TotalXSecData.txt";
  string filename_material1_ela = "neutronHe4TotalXSecData.txt";
  string filename_material1_ine = "neutronHe3ProtonXSecData.txt";

  string filename_material2_tot = "neutronB10TotalXSecData.txt";
  string filename_material2_ela = "neutronGd157TotalXSecData.txt";
  string filename_material2_ine = "neutronLi6NonelasticXSecData.txt";

  //filename array 
  string filenames[size] = {filename_material1_tot,filename_material2_tot,filename_material1_ela,filename_material2_ela,filename_material1_ine,filename_material2_ine};

  //create an array of the entries
  int entries[size];

  //loop over file array once to get number of entries for each
  for(int i=0;i<size;i++){

    ifstream openfile(filenames[i].c_str());

    int counter = 0;

    std::string line ="";  
    while(std::getline(openfile,line))
    {
        std::istringstream is(line);
        double energy;
        double xsec;

        std::string energyAsString;
        std::string xsecAsString;

        is >> energyAsString >> xsecAsString;
	if( energyAsString == "")continue;
//	cout << "\nEnergy: " << energyAsString << " [eV], XSec: " << xsecAsString << " [barns]"; 
	counter++;

    }
    entries[i] = counter;
    cout << "\n----- Entries: " << counter <<std::endl;
  }

  //make arrays for xsec 
  double energy1[entries[0]];
  double xsec1[entries[0]];

  double energy2[entries[1]];
  double xsec2[entries[1]];

  double energy3[entries[2]];
  double xsec3[entries[2]];

  double energy4[entries[3]];
  double xsec4[entries[3]];

  if(size==6){
    double energy5[entries[4]];
    double xsec5[entries[4]];

    double energy6[entries[5]];
    double xsec6[entries[5]];
  }

  //loop over file array again to get values
  for(int i=0;i<size;i++){

    ifstream openfile(filenames[i].c_str());

    int index = 0;
    std::string line ="";  
    while(std::getline(openfile,line))
    {
        std::istringstream is(line);
        double energy;
        double xsec;

        std::string energyAsString;
        std::string xsecAsString;

        is >> energyAsString >> xsecAsString;
	if( energyAsString == "")continue;

	//loop over energy string to turn to double
	for(int j=0; j<energyAsString.length();++j){
			if( energyAsString[j] == 'e')break;
			if( energyAsString[j] == 'E')break;
			if( (energyAsString[j] == '-') || (energyAsString[j] == '+') ){
				energyAsString[j-1] ='e';
			}
	}

	//loop over xsec string to turn to double
	for(int j=0; j<xsecAsString.length();++j){
			if( xsecAsString[j] == 'e')break;
			if( xsecAsString[j] == 'E')break;
			if( (xsecAsString[j] == '-') || (xsecAsString[j] == '+') ){
				xsecAsString[j-1] ='e';
			}
	}

	std::istringstream stringToDouble(energyAsString);	
	stringToDouble >> energy;
	stringToDouble.clear();
	stringToDouble.str(xsecAsString);
	stringToDouble >> xsec;
	stringToDouble.clear();	

	//cout << "\nEnergy = " << energy << "[eV], Xsec = " << xsec <<"[barns]";

	if(i==0){
		energy1[index] = energy;
		xsec1[index] = xsec;
	}
	if(i==1){
		energy2[index] = energy;
		xsec2[index] = xsec;
	}
	if(i==2){
		energy3[index] = energy;
		xsec3[index] = xsec;
	}
	if(i==3){
		energy4[index] = energy;
		xsec4[index] = xsec;
	}
	if(i==4){
		energy5[index] = energy;
		xsec5[index] = xsec;
	}
	if(i==5){
		energy6[index] = energy;
		xsec6[index] = xsec;
	}

	index++;
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();

  //make the graphs
  TGraph* gr_mat1_tot = new TGraph(entries[0],energy1,xsec1);
  TGraph* gr_mat2_tot = new TGraph(entries[1],energy2,xsec2);
  TGraph* gr_mat1_ela = new TGraph(entries[2],energy3,xsec3);
  TGraph* gr_mat2_ela = new TGraph(entries[3],energy4,xsec4);
  TGraph* gr_mat1_pro = new TGraph(entries[4],energy5,xsec5);
  TGraph* gr_mat2_ine = new TGraph(entries[5],energy6,xsec6);

  //draw on multigraph
  gPad->SetTickx();  
  gPad->SetTicky();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->Update();
  TMultiGraph *mg = new TMultiGraph();

  int linewidth =2;
  gr_mat1_tot->SetLineColor(kRed);
  gr_mat1_tot->SetFillColor(0);
  gr_mat1_tot->SetLineWidth(linewidth);
  gr_mat1_tot->SetTitle("^{3}He");
  gr_mat2_tot->SetFillColor(0);
  gr_mat2_tot->SetLineColor(kRed);
  gr_mat2_tot->SetTitle("^{10}B");
  gr_mat2_tot->SetLineWidth(linewidth);
  gr_mat1_ela->SetLineColor(kBlue);
  gr_mat1_ela->SetFillColor(0);
  gr_mat1_ela->SetLineStyle(1);
  gr_mat1_ela->SetTitle("^{4}He");
  gr_mat1_ela->SetLineWidth(linewidth);
  gr_mat2_ela->SetFillColor(0);
  gr_mat2_ela->SetLineColor(kBlue);
  gr_mat2_ela->SetLineStyle(1);
  gr_mat2_ela->SetLineWidth(linewidth-1);
  gr_mat2_ela->SetTitle("^{157}Gd");
  gr_mat2_ine->SetLineColor(kOrange);
  gr_mat2_ine->SetLineStyle(2);
  gr_mat2_ine->SetLineWidth(linewidth);
  mg->Add(gr_mat1_tot);
//  mg->Add(gr_mat2_tot);
  mg->Add(gr_mat1_ela);
//  mg->Add(gr_mat2_ela);
//  mg->Add(gr_mat1_pro);
//  mg->Add(gr_mat2_ine); 

  mg->Draw("A");
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitleOffset(1.1);
  mg->GetXaxis()->SetTitle("Energy [eV]");
  mg->GetYaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitle("Cross Section [barns]");
  mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetYaxis()->SetRangeUser(1e-1,1e6);
  mg->Draw("C");

}
