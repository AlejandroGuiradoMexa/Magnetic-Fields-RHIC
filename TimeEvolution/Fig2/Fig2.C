
void StylePlots(){
  gStyle->SetOptFit(11111);
  gStyle->Reset("Plain");
  gStyle->SetOptStat(0);

  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10); //
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.18); //0.12
  gStyle->SetPadLeftMargin(0.18); // 0.12
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);

  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.06,"xyz"); //0.5
  //  gStyle->SetTitleOffset(1.5,"y");  //0.95
  //  gStyle->SetTitleOffset(1.3,"x"); //0-95

 gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);

  gStyle->SetLegendBorderSize(0); //
  
  gStyle->SetLegendFillColor(0);
    gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
  }

void Fig2(){
 StylePlots();

  const int tI = 0.0;
  const int tF = 10.0;
  const float tSteps = 0.1;
  const int tTotal = 101;
  const int protonNumber1 = 166.0;
  const int protonNumber2 = 158.0;

  gROOT->Reset();
  TChain mychainBiBi("T");
    mychainBiBi.Add("../Datos/BiBi4GeV.root");
  TChain mychainAuAu("T");
    mychainAuAu.Add("../Datos/AuAu4GeV.root");
  TChain mychainBiBi0("T");
    mychainBiBi0.Add("../Datos/BiBi1000_4GeV_t0.root");
  TChain mychainAuAu0("T");
    mychainAuAu0.Add("../Datos/AuAu1000_4GeV_t0.root");
    struct particula_t 
    {
      Float_t time,X,Y,Z,E,Px,Py,Pz,Pt,P,m,id,isoespin,charge,lastcoll,numbercoll,history,frezetime,frezeX,frezeY,frezeZ,frezeE,frezePx,frezePy,frezePz,b,nspec,R,PXR,eBx,eBy,eBz;
    } PARTICLE;
    particula_t  particle;
    
    mychainBiBi.SetBranchAddress("particle",&particle);
    mychainAuAu.SetBranchAddress("particle",&particle);
    mychainBiBi0.SetBranchAddress("particle",&particle);
    mychainAuAu0.SetBranchAddress("particle",&particle);

    Int_t nevent1 = mychainBiBi.GetEntries();
    Int_t nevent2 = mychainAuAu.GetEntries();
    Int_t nevent1_0 = mychainBiBi0.GetEntries();
    Int_t nevent2_0 = mychainAuAu0.GetEntries();

    const int NumberOfDivisions = 3.0;
    double_t aImpacParI[NumberOfDivisions] = {0,7,12};
    double_t aImpacParF[NumberOfDivisions] = {5,10,16};

    double_t aTime[tTotal] = {0.};    
    double_t particlesTotal1[NumberOfDivisions] = {0.},particlesTotal2[NumberOfDivisions] = {0.};
    double_t particlesTotal1_0[NumberOfDivisions] = {0.},particlesTotal2_0[NumberOfDivisions] = {0.};

    double_t aNorm1[NumberOfDivisions] = {0.};
    double_t aNorm1_0[NumberOfDivisions] = {0.};
    double_t aNorm2[NumberOfDivisions] = {0.};
    double_t aNorm2_0[NumberOfDivisions] = {0.};
    
    double eBy1AuAu[tTotal+1] = {}, eBy2AuAu[tTotal+1] = {}, eBy3AuAu[tTotal+1] = {};
    double eBy1BiBi[tTotal+1] = {}, eBy2BiBi[tTotal+1] = {}, eBy3BiBi[tTotal+1] = {};
    double_t  eBy1_0[NumberOfDivisions] = {0.}, eBy2_0[NumberOfDivisions] = {0.};
                

    //Calculate Number of particles in each range of impact parameter
    for(Int_t j=0; j<nevent1; j++){
      mychainBiBi.GetEvent(j);
        if(particle.time == tSteps && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            for( Int_t bID = 0; bID<NumberOfDivisions; bID++ ){
                if(particle.b > aImpacParI[bID] && particle.b < aImpacParF[bID] )
                    particlesTotal1[bID]++;
            }
        }
    }
    
    for(Int_t j=0; j<nevent2; j++){
      mychainAuAu.GetEvent(j);
        if(particle.time == tSteps && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            for( Int_t bID = 0; bID<NumberOfDivisions; bID++ ){
                if(particle.b > aImpacParI[bID] && particle.b < aImpacParF[bID] ){
                    particlesTotal2[bID]++;
                }
            }
        }
    }

    for(Int_t j=0; j<nevent1_0; j++){
      mychainBiBi0.GetEvent(j);
        if(particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            for( Int_t bID = 0; bID<NumberOfDivisions; bID++ ){
                if(particle.b > aImpacParI[bID] && particle.b < aImpacParF[bID] ){
                    particlesTotal1_0[bID]++;
                }
            }
        }
    }
    
    for(Int_t j=0; j<nevent2_0; j++){
      mychainAuAu0.GetEvent(j);
        if(particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            for( Int_t bID = 0; bID<NumberOfDivisions; bID++ ){
                if(particle.b > aImpacParI[bID] && particle.b < aImpacParF[bID] )
                    particlesTotal2_0[bID]++;
            }
        }
    }

    //Divide number of particles between the number of protons of each ion to
    //get the number of events
    for( Int_t bID = 0; bID<NumberOfDivisions; bID++ ){
                    particlesTotal1[bID] = particlesTotal1[bID]/protonNumber1;
                    particlesTotal2[bID] = particlesTotal2[bID]/protonNumber2;
                    }
                    
    for( Int_t bID = 0; bID<NumberOfDivisions; bID++ ){
                    particlesTotal1_0[bID] = particlesTotal1_0[bID]/protonNumber1;
            particlesTotal2_0[bID] = particlesTotal2_0[bID]/protonNumber2;
    }

    //Calculate the normalization factor of events and alpha EM
    for( Int_t NormID = 0; NormID<NumberOfDivisions; NormID++ ){
        aNorm1[NormID] = particlesTotal1[NormID]*137.;
        aNorm1_0[NormID] = particlesTotal1_0[NormID]*137.;
        aNorm2[NormID] = particlesTotal2[NormID]*137.;
        aNorm2_0[NormID] = particlesTotal2_0[NormID]*137.;
    }

    

    //Getting the magnetic field in y direction
    
    //Bi+Bi
    for (Int_t iEv=0; iEv<nevent1_0; iEv++) 
    {
        mychainBiBi0.GetEvent(iEv);
        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){
            for( Int_t bID = 0; bID<NumberOfDivisions; bID++ ){
                if(particle.b > aImpacParI[bID] && particle.b < aImpacParF[bID] ){
                    eBy1_0[bID] += particle.eBy;
                }
            }
        }
    }    

    
    for (Int_t iEv=0; iEv<nevent1; iEv++) 
    {
        mychainBiBi.GetEvent(iEv);

        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){
            for( Int_t iTime = 1; iTime<tTotal; ++iTime){
                Float_t iTC = iTime*0.1;
                if(particle.time == iTC){
                    if(particle.b>aImpacParI[0] && particle.b<aImpacParF[0] ){
                        eBy1BiBi[iTime] += particle.eBy;
                    }
                    if(particle.b>aImpacParI[1] && particle.b<aImpacParF[1] ){
                        eBy2BiBi[iTime] += particle.eBy;
                    }
                    if(particle.b>aImpacParI[2] && particle.b<aImpacParF[2] ){
                        eBy3BiBi[iTime] += particle.eBy;
                    }
                }
            }
        }
    }     

    
    //Au+Au
    for (Int_t iEv=0; iEv<nevent2_0; iEv++) 
    {
        mychainAuAu0.GetEvent(iEv);
        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){
            for( Int_t bID = 0; bID<NumberOfDivisions; bID++ ){
                if(particle.b > aImpacParI[bID] && particle.b < aImpacParF[bID] ){
                    eBy2_0[bID] += particle.eBy;
                }
            }
        }
    }    

    double Counter = 0.;

    for (Int_t iEv=0; iEv<nevent2; iEv++)
    {
        mychainAuAu.GetEvent(iEv);
        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){

            for( Int_t iTime = 1; iTime<tTotal; ++iTime){
                Float_t iTC = iTime*0.1;
//                  cout << particle.time <<"   " << iTC<< "   " << particle.eBy<<"   " << particle.b <<endl;

                if(particle.time == iTC){

//                     cout << particle.time << "   " << iTC << "   " << particle.eBy<< "   " << Counter<<endl;

                    if(particle.b>aImpacParI[0] && particle.b<aImpacParF[0] ){
                        eBy1AuAu[iTime] += particle.eBy;
                    }
                    if(particle.b>aImpacParI[1] && particle.b<aImpacParF[1] ){
                        eBy2AuAu[iTime] += particle.eBy;
                    }
                    if(particle.b>aImpacParI[2] && particle.b<aImpacParF[2] ){
                        eBy3AuAu[iTime] += particle.eBy;
                    }
                }
            }
        }
    }   
  
    //Normalize MF 

    eBy1BiBi[0] =eBy1_0[0]/aNorm1_0[0];
    eBy2BiBi[0] =eBy1_0[1]/aNorm1_0[1];
    eBy3BiBi[0] =eBy1_0[2]/aNorm1_0[2];
    
    eBy1AuAu[0] =eBy2_0[0]/aNorm2_0[0];
    eBy2AuAu[0] =eBy2_0[1]/aNorm2_0[1];
    eBy3AuAu[0] =eBy2_0[2]/aNorm2_0[2];
  
    for( Int_t iTime = 1; iTime<tTotal; ++iTime){
        eBy1BiBi[iTime]= eBy1BiBi[iTime]/aNorm1[0];
        eBy2BiBi[iTime]= eBy2BiBi[iTime]/aNorm1[1];
        eBy3BiBi[iTime]= eBy3BiBi[iTime]/aNorm1[2];
        
        eBy1AuAu[iTime]= eBy1AuAu[iTime]/aNorm2[0];
        eBy2AuAu[iTime]= eBy2AuAu[iTime]/aNorm2[1];
        eBy3AuAu[iTime]= eBy3AuAu[iTime]/aNorm2[2];            
        }
    
    for( Int_t iTime = 0; iTime<tTotal; ++iTime){
        double_t iTC = iTime*tSteps;
        aTime[iTime] = iTC;
    }    

    
    for( Int_t iTime = 0; iTime<tTotal; ++iTime){
       cout << "  "<<  eBy1BiBi[iTime] <<"  "<<  eBy2BiBi[iTime] <<"  "<<  eBy3BiBi[iTime] <<"  "<<  iTime <<endl;
    }

    for( Int_t iTime = 0; iTime<tTotal; ++iTime){
       cout << "  "<<  eBy1AuAu[iTime] <<"  "<<  eBy2AuAu[iTime] <<"  "<<  eBy3AuAu[iTime] <<"  "<<  iTime <<endl;
    }
  
   TCanvas *c1 =  new TCanvas("c1","",1024,640);
      c1->SetFillStyle(4000);    
    gStyle->SetOptStat(false);
    c1->SetRightMargin(0.0465116);
    c1->SetTopMargin(0.1);
    c1->SetFillColor(0);
    c1->Divide(2,1);
    c1->cd(1);
    TGraph *gr = new TGraph(tTotal,aTime,eBy1BiBi); 
   gr->SetLineColor(1);
   gr->SetLineWidth(2);
   gr->SetMarkerColor(1);
   gr->SetMarkerStyle(21);
   gr->SetMarkerSize(1.5);
   gr->SetTitle("");
   gr->GetXaxis()->SetTitle("Time[fm]");
   gr->GetYaxis()->SetTitle("-B_{y}/m_{#pi}^{2}");
   gr->GetXaxis()->SetLimits(0.,11.);
   gr->SetMinimum(0.);
   gr->SetMaximum(0.35);
   gr->Draw("AL");

   
   TGraph *gr1 = new TGraph(tTotal,aTime,eBy2BiBi);
   gr1->SetLineColor(2);
   gr1->SetLineWidth(2);
   gr1->SetMarkerColor(1);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.5);
   gr1->Draw("L");
   
   TGraph *gr2 = new TGraph(tTotal,aTime,eBy3BiBi);
   gr2->SetLineColor(3);
   gr2->SetLineWidth(2);
   gr2->SetMarkerColor(1);
   gr2->SetMarkerStyle(23);
   gr2->SetMarkerSize(1.5);
   gr2->Draw("L");
 
   
    TLegend *leg1 = new TLegend(0.60,0.60,0.80,0.85);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.045);
    leg1->SetLineColor(0);
    leg1->SetLineStyle(0);
    leg1->SetLineWidth(2);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1012);
    leg1->SetHeader("Impact Parameter[fm]","C");
    leg1->AddEntry(gr,"0-5","L");
    leg1->AddEntry(gr1,"7-10","L");
    leg1->AddEntry(gr2,"12-16"  ,"L");
    leg1->Draw();

    
    TLegend *leg2 = new TLegend(0.60,0.42,0.80,0.54);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.055);
    leg2->SetLineColor(0);
    leg2->SetLineStyle(0);
    leg2->SetLineWidth(2);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1012);
    leg2->SetHeader("BiBi, #sqrt{S_{NN}}=4GeV ","C");
    leg2->Draw(); 
    c1->cd(2);

   TGraph *grA = new TGraph(tTotal,aTime,eBy1AuAu);
   grA->SetLineColor(1);
   grA->SetLineWidth(2);
   grA->SetMarkerColor(1);
   grA->SetMarkerStyle(21);
   grA->SetMarkerSize(1.5);
   grA->SetTitle("");
   grA->GetXaxis()->SetTitle("Time[fm]");
   grA->GetYaxis()->SetTitle("-B_{y}/m_{#pi}^{2}");
   grA->GetXaxis()->SetLimits(0.,11.);
   grA->SetMinimum(0.);
   grA->SetMaximum(0.35);
   grA->Draw("AL");
   
   TGraph *gr1A = new TGraph(tTotal,aTime,eBy2AuAu);
   gr1A->SetLineColor(2);
   gr1A->SetLineWidth(2);
   gr1A->SetMarkerColor(1);
   gr1A->SetMarkerStyle(22);
   gr1A->SetMarkerSize(1.5);
   gr1A->Draw("L");
   
   TGraph *gr2A = new TGraph(tTotal,aTime,eBy3AuAu);
   gr2A->SetLineColor(3);
   gr2A->SetLineWidth(2);
   gr2A->SetMarkerColor(1);
   gr2A->SetMarkerStyle(23);
   gr2A->SetMarkerSize(1.5);
   gr2A->Draw("L");
   
    TLegend *leg1A = new TLegend(0.60,0.60,0.80,0.85);
    leg1A->SetTextFont(42);
    leg1A->SetTextSize(0.045);
    leg1A->SetLineColor(0);
    leg1A->SetLineStyle(0);
    leg1A->SetLineWidth(2);
    leg1A->SetFillColor(0);
    leg1A->SetFillStyle(1012);
    leg1A->SetHeader("Magnetic Field","C");
    leg1A->SetHeader("Impact Parameter[fm]","C");
    leg1A->AddEntry(gr,"0-5","L");
    leg1A->AddEntry(gr1,"7-10","L");
    leg1A->AddEntry(gr2,"12-16"  ,"L");
    leg1A->Draw();

    TLegend *leg2A = new TLegend(0.60,0.42,0.80,0.54);
    leg2A->SetTextFont(42);
    leg2A->SetTextSize(0.055);
    leg2A->SetLineColor(0);
    leg2A->SetLineStyle(0);
    leg2A->SetLineWidth(2);
    leg2A->SetFillColor(0);
    leg2A->SetFillStyle(1012);
    leg2A->SetHeader("AuAu,#sqrt{S_{NN}}=4GeV ","C");
    leg2A->Draw(); 


   c1->SaveAs("fig2-4GeV.pdf");
   c1->SaveAs("fig2-4GeV.png");


  }
