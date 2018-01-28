{
  int size =6;

  //read in files
  string filename_1 = "isochoric_2o0gml.txt";
  string filename_2 = "isochoric_0o4gml.txt";
  string filename_3 = "isochoric_0o8gml.txt";
  string filename_4 = "isochoric_1o0gml.txt";
  string filename_5 = "isochoric_1o2gml.txt";
  string filename_6 = "isochoric_1o6gml.txt";

  //filename array 
  string filenames[size] = {filename_1,filename_2,filename_3,filename_4,filename_5,filename_6};

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
        double pressure;
        double temperature;
	std::string type;
	std::string col3, col4, col5,col6,col7,col8,col9,col10,col11,col12,col13,col14;

        std::string tempAsString;
        std::string pressureAsString;

        is >> tempAsString >> pressureAsString >> col3 >> col4 >> col5 >> col6 >> col7
		>> col8 >> col9 >>col10 >> col11 >> col12 >> col13>> col14;
	if( col14 == "liquid" && i!=0 )continue;
	counter++;

    }
    entries[i] = counter;
    cout << "\n----- Entries: " << counter <<std::endl;
  }

  //make arrays for pressure 
  double temp1[entries[0]];
  double pressure1[entries[0]];

  double temp2[entries[1]];
  double pressure2[entries[1]];

  double temp3[entries[2]];
  double pressure3[entries[2]];

  double temp4[entries[3]];
  double pressure4[entries[3]];

  if(size==6){
    double temp5[entries[4]];
    double pressure5[entries[4]];

    double temp6[entries[5]];
    double pressure6[entries[5]];
  }

  //loop over file array again to get values
  for(int i=0;i<size;i++){

    ifstream openfile(filenames[i].c_str());

    int index = 0;
    std::string line ="";  
    while(std::getline(openfile,line))
    {
        std::istringstream is(line);
        double pressure;
        double temperature;
	std::string type;
	std::string col3, col4, col5,col6,col7,col8,col9,col10,col11,col12,col13,col14;

        std::string tempAsString;
        std::string pressureAsString;

        is >> tempAsString >> pressureAsString >> col3 >> col4 >> col5 >> col6 >> col7
		>> col8 >> col9 >>col10 >> col11 >> col12 >> col13>> col14;
	if( col14 == "liquid" && i!=0 )continue;

	std::istringstream stringToDouble(tempAsString);	
	stringToDouble >> temperature;
	stringToDouble.clear();
	stringToDouble.str(pressureAsString);
	stringToDouble >> pressure;
	stringToDouble.clear();	

	if(i==0){
		temp1[index] = temperature;
		pressure1[index] = pressure;
	}
	if(i==1){
		temp2[index] = temperature;
		pressure2[index] = pressure;
	}
	if(i==2){
		temp3[index] = temperature;
		pressure3[index] = pressure;
	}
	if(i==3){
		temp4[index] = temperature;
		pressure4[index] = pressure;
	}
	if(i==4){
		temp5[index] = temperature;
		pressure5[index] = pressure;
	}
	if(i==5){
		temp6[index] = temperature;
		pressure6[index] = pressure;
	}

	index++;
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();

  //make the graphs
  TGraph* gr1 = new TGraph(entries[0],temp1,pressure1);
  TGraph* gr2 = new TGraph(entries[1],temp2,pressure2);
  TGraph* gr3 = new TGraph(entries[2],temp3,pressure3);
  TGraph* gr4 = new TGraph(entries[3],temp4,pressure4);
  TGraph* gr5 = new TGraph(entries[4],temp5,pressure5);
  TGraph* gr6 = new TGraph(entries[5],temp6,pressure6);

  //draw on multigraph
  gPad->SetTickx();  
  gPad->SetTicky();
//  gPad->SetLogx();
//  gPad->SetLogy();
  gPad->Update();
  TMultiGraph *mg = new TMultiGraph();

  int linewidth =2;
  gr1->SetLineColor(kBlack);
  gr1->SetLineStyle(1);
  gr1->SetFillColor(0);
  gr1->SetLineWidth(linewidth);
  gr1->SetTitle("2.0 g/mL");

  gr2->SetFillColor(0);
  gr2->SetLineColor(kRed);
  gr2->SetLineStyle(1);
  gr2->SetTitle("0.4 g/mL");
  gr2->SetLineWidth(linewidth);

  gr3->SetLineColor(kBlue);
  gr3->SetFillColor(0);
  gr3->SetLineStyle(1);
  gr3->SetTitle("0.8 g/mL");
  gr3->SetLineWidth(linewidth);

  gr4->SetFillColor(0);
  gr4->SetLineColor(kBlack);
  gr4->SetLineStyle(2);
  gr4->SetLineWidth(linewidth);
  gr4->SetTitle("1.0 g/mL");

  gr5->SetTitle("1.2 g/mL");
  gr5->SetLineColor(kRed);
  gr5->SetLineStyle(2);
  gr5->SetLineWidth(linewidth);
  gr5->SetFillColor(0);

  gr6->SetTitle("1.6 g/mL");
  gr6->SetLineColor(kOrange-3);
  gr6->SetLineStyle(2);
  gr6->SetLineWidth(linewidth);
  gr6->SetFillColor(0);

  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr3);
  mg->Add(gr4);
  mg->Add(gr5);
  mg->Add(gr6); 

  mg->Draw("A");
  mg->GetXaxis()->SetTitle("Temperature [#circ C]");
  mg->GetYaxis()->SetTitle("Pressure [bar]");
  mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->SetRangeUser(0,50);
  mg->GetYaxis()->SetRangeUser(0,200);
  mg->Draw("C");

}
