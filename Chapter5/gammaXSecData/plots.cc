{
  int size =1;

  //read in files
  string filename_1 = "gammaHeXSecData.txt";

  //filename array 
  string filenames[1] = {filename_1};

  //create an array of the entries
  int entries[1];

  //loop over file array once to get number of entries for each
  for(int i=0;i<size;i++){

    ifstream openfile(filenames[i].c_str());

    int counter = 0;

    std::string line ="";  
    while(std::getline(openfile,line))
    {
        std::istringstream is(line);
        double energy;
	double cohXSec,inCohXSec,pairProd1XSec,pairProd2XSec;
        double comptonXSec,photoXSec,pairProdXSec,totalXSec;

        std::string energyAsString,cohScatAsString,inCohAsString,pairProd1AsString,pairProd2AsString;
        std::string photoAsString,totalAsString;

        is >> energyAsString >> cohScatAsString >> inCohAsString >> photoAsString >> pairProd1AsString >> pairProd2AsString >> totalAsString;
	counter++;

    }
    entries[i] = counter;
    cout << "\n----- Entries: " << counter <<std::endl;
  }

  //make arrays for xsecs 
  double energy1[entries[0]];
  double compXSec1[entries[0]];
  double pairXSec1[entries[0]];
  double photoXSec1[entries[0]];
  double totalXSec1[entries[0]];

  //loop over file array again to get values
  for(int i=0;i<size;i++){

    ifstream openfile(filenames[i].c_str());

    int index = 0;
    std::string line ="";  
    while(std::getline(openfile,line))
    {
        std::istringstream is(line);
        double energy;
	double cohXSec,inCohXSec,pairProd1XSec,pairProd2XSec;
        double comptonXSec,photoXSec,pairProdXSec,totalXSec;

        std::string energyAsString,cohScatAsString,inCohAsString,pairProd1AsString,pairProd2AsString;
        std::string photoAsString,totalAsString;

        is >> energyAsString >> cohScatAsString >> inCohAsString >> photoAsString >> pairProd1AsString >> pairProd2AsString >> totalAsString;

	std::istringstream stringToDouble(energyAsString);	
	stringToDouble >> energy;
	stringToDouble.clear();
	stringToDouble.str(cohScatAsString);
	stringToDouble >> cohXSec;
	stringToDouble.clear();	
	stringToDouble.str(inCohAsString);
	stringToDouble >> inCohXSec;
	stringToDouble.clear();	
	stringToDouble.str(pairProd1AsString);
	stringToDouble >> pairProd1XSec;
	stringToDouble.clear();	
	stringToDouble.str(pairProd2AsString);
	stringToDouble >> pairProd2XSec;
	stringToDouble.clear();	
	stringToDouble.str(photoAsString);
	stringToDouble >> photoXSec;
	stringToDouble.clear();	
	stringToDouble.str(totalAsString);
	stringToDouble >> totalXSec;
	stringToDouble.clear();	

	if(i==0){
		energy1[index] = energy;
  		compXSec1[index] = inCohXSec;//cohXSec + inCohXSec;
 		pairXSec1[index] = pairProd1XSec + pairProd2XSec;
 		photoXSec1[index] = photoXSec;
		totalXSec1[index] = totalXSec;
	}

	index++;
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();

  //make the graphs
  TGraph* gr1 = new TGraph(entries[0],energy1,compXSec1);
  TGraph* gr2 = new TGraph(entries[0],energy1,pairXSec1);
  TGraph* gr3 = new TGraph(entries[0],energy1,photoXSec1);
  TGraph* gr4 = new TGraph(entries[0],energy1,totalXSec1);

  //draw on multigraph
  gPad->SetTickx();  
  gPad->SetTicky();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->Update();
  TMultiGraph *mg = new TMultiGraph();

  int linewidth =2;

  gr1->SetLineColor(kBlack);
  gr1->SetLineStyle(1);
  gr1->SetFillColor(0);
  gr1->SetLineWidth(linewidth);
  gr1->SetTitle("compton");

  gr2->SetFillColor(0);
  gr2->SetLineColor(kRed);
  gr2->SetLineStyle(1);
  gr2->SetTitle("pair production");
  gr2->SetLineWidth(linewidth);

  gr3->SetLineColor(kOrange-3);
  gr3->SetFillColor(0);
  gr3->SetLineStyle(1);
  gr3->SetTitle("photoelectric");
  gr3->SetLineWidth(1);//linewidth);

  gr4->SetFillColor(0);
  gr4->SetLineColor(kBlue);
  gr4->SetLineStyle(1);
  gr4->SetLineWidth(linewidth);
  gr4->SetTitle("total");

  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr3,"1pl");
  mg->Add(gr4,"1pl");

  mg->Draw("A");
  mg->GetXaxis()->SetTitle("Energy [MeV]");
  mg->GetYaxis()->SetTitle("Cross section [cm^{2}/g]");
  mg->GetYaxis()->SetTitleOffset(1.2);
//  mg->GetXaxis()->SetRangeUser(0,50);
  mg->GetYaxis()->SetRangeUser(0.01,1e5);
  mg->Draw("C");

}
