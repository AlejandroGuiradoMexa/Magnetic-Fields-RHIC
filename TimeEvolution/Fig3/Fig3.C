
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
  gStyle->SetPadTopMargin(0.10);
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
gStyle->SetFrameFillStyle(0);
  gStyle->SetLegendFillColor(0);
    gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
  }

void Fig3(){
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
    
    double eBy1AuAu[tTotal] = {}, eBy2AuAu[tTotal] = {}, eBy3AuAu[tTotal] = {};
    double eBy1BiBi[tTotal] = {}, eBy2BiBi[tTotal] = {}, eBy3BiBi[tTotal] = {};
    double_t  eBy1_0[NumberOfDivisions] = {0.}, eBy2_0[NumberOfDivisions] = {0.};
    
                

    /**************************************************************************************
     * 
     *                                   SPECTATORS 
     * 
     ***************************************************************************************
     */
    
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


    for (Int_t iEv=0; iEv<nevent2; iEv++)
    {
        mychainAuAu.GetEvent(iEv);
        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){
            for( Int_t iTime = 1; iTime<tTotal; ++iTime){
                Float_t iTC = iTime*0.1;
                if(particle.time == iTC){
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
  
      /**************************************************************************************
     * 
     *                                   PARTICIPANTS 
     * 
     ***************************************************************************************
     */
    double eBy1AuAuP[tTotal] = {}, eBy2AuAuP[tTotal] = {}, eBy3AuAuP[tTotal] = {};
    double_t particlesTotal1P[NumberOfDivisions] = {0.},particlesTotal2P[NumberOfDivisions] = {0.};

     for (Int_t iEv=0; iEv<nevent2; iEv++)
    {
        mychainAuAu.GetEvent(iEv);
        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        if(particle.numbercoll!=0 && particle.charge==1 && particle.id == 1 && R>0.3){
            for( Int_t iTime = 1; iTime<tTotal; ++iTime){
                Float_t iTC = iTime*0.1;
                if(particle.time == iTC){
                    if(isnan(particle.eBy) == true )continue;
                    else{ 
                        if(particle.b>aImpacParI[0] && particle.b<aImpacParF[0] ){
                            eBy1AuAuP[iTime] += particle.eBy;
                        }
                        if(particle.b>aImpacParI[1] && particle.b<aImpacParF[1] ){
                            eBy2AuAuP[iTime] += particle.eBy;
                        }
                        if(particle.b>aImpacParI[2] && particle.b<aImpacParF[2] ){
                            eBy3AuAuP[iTime] += particle.eBy;
                        }
                    }
                }
            }
        }
    }   
        
    eBy1AuAuP[0] = 0.;
    eBy2AuAuP[0] = 0.;
    eBy3AuAuP[0] = 0.;
      
    for( Int_t iTime = 1; iTime<tTotal; ++iTime){
        eBy1AuAuP[iTime]= eBy1AuAuP[iTime]/aNorm2[0];
        eBy2AuAuP[iTime]= eBy2AuAuP[iTime]/aNorm2[1];
        eBy3AuAuP[iTime]= eBy3AuAuP[iTime]/aNorm2[2];            
        }
    
    
    cout << "----------------------------------------------------------" <<endl;
    cout << "                  PARTICIPANTS                            " <<endl;
    cout << "----------------------------------------------------------" <<endl;
    for( Int_t iTime = 0; iTime<tTotal; ++iTime){
       cout << "  "<<  eBy1AuAuP[iTime] <<"  "<<  eBy2AuAuP[iTime] <<"  "<<  eBy3AuAuP[iTime] <<"  "<<  iTime <<endl;
    }
  
    gStyle->SetFrameFillStyle(0);

   TCanvas *c1 =  new TCanvas("c1","",1024,640);
//       c1->SetFillStyle(4000);    
//     gStyle->SetOptStat(false);
    c1->SetRightMargin(0.0465116);
    c1->SetTopMargin(0.1);
    c1->SetFillColor(0);
    c1->Divide(3,1);
    c1->cd(1);
gPad->SetFillStyle(0);
TGraph *gr = new TGraph(tTotal,aTime,eBy1AuAu); 
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
   gr->SetMaximum(0.3);
   gr->Draw("AL");

   
   TGraph *gr1 = new TGraph(tTotal,aTime,eBy1AuAuP);
   gr1->SetLineColor(2);
   gr1->SetLineWidth(2);
   gr1->SetMarkerColor(1);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.5);
   gr1->Draw("L");


    TLegend *leg1 = new TLegend(0.60,0.50,0.80,0.70);
    leg1->SetTextFont(62);
    leg1->SetTextSize(0.053);
    leg1->SetLineColor(0);
    leg1->SetLineStyle(0);
    leg1->SetLineWidth(2);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1012);
    leg1->AddEntry(gr1,"Participants","l");
    leg1->AddEntry(gr,"Spectators","l");
    leg1->Draw();
    TLegend *legHead = new TLegend(0.3,0.80,0.50,0.85);
    legHead->SetTextFont(62);
    legHead->SetTextSize(0.055);
    legHead->SetLineColor(0);
    legHead->SetLineStyle(0);
    legHead->SetLineWidth(2);
    legHead->SetFillColor(0);
    legHead->SetFillStyle(1012);
    legHead->SetHeader("a)  b=0-5[fm]","L");
    legHead->Draw(); 

    
     c1->cd(2);
gPad->SetFillStyle(0);
TGraph *gr2 = new TGraph(tTotal,aTime,eBy2AuAu); 
   gr2->SetLineColor(1);
   gr2->SetLineWidth(2);
   gr2->SetMarkerColor(1);
   gr2->SetMarkerStyle(21);
   gr2->SetMarkerSize(1.5);
   gr2->SetTitle("");
   gr2->GetXaxis()->SetTitle("Time[fm]");
   gr2->GetYaxis()->SetTitle("-B_{y}/m_{#pi}^{2}");
   gr2->GetXaxis()->SetLimits(0.,11.);
   gr2->SetMinimum(0.);
   gr2->SetMaximum(0.30);
   gr2->Draw("AL");

   
   TGraph *gr3 = new TGraph(tTotal,aTime,eBy2AuAuP);
   gr3->SetLineColor(2);
   gr3->SetLineWidth(2);
   gr3->SetMarkerColor(1);
   gr3->SetMarkerStyle(22);
   gr3->SetMarkerSize(1.5);
   gr3->Draw("L");
   
   
    TLegend *leg2 = new TLegend(0.60,0.50,0.80,0.7);
    leg2->SetTextFont(62);
    leg2->SetTextSize(0.053);
    leg2->SetLineColor(0);
    leg2->SetLineStyle(0);
    leg2->SetLineWidth(2);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1012);
    leg2->AddEntry(gr3,"Participants","l");
    leg2->AddEntry(gr2,"Spectators","l");
    leg2->Draw();
    TLegend *legHead2 = new TLegend(0.3,0.80,0.50,0.85);
    legHead2->SetTextFont(62);
    legHead2->SetTextSize(0.055);
    legHead2->SetLineColor(0);
    legHead2->SetLineStyle(0);
    legHead2->SetLineWidth(2);
    legHead2->SetFillColor(0);
    legHead2->SetFillStyle(1012);
    legHead2->SetHeader("b)  b=7-10[fm]","L");
    legHead2->Draw(); 

     c1->cd(3);
gPad->SetFillStyle(0);
TGraph *gr4 = new TGraph(tTotal,aTime,eBy3AuAu); 
   gr4->SetLineColor(1);
   gr4->SetLineWidth(2);
   gr4->SetMarkerColor(1);
   gr4->SetMarkerStyle(21);
   gr4->SetMarkerSize(1.5);
   gr4->SetTitle("");
   gr4->GetXaxis()->SetTitle("Time[fm]");
   gr4->GetYaxis()->SetTitle("-B_{y}/m_{#pi}^{2}");
   gr4->GetXaxis()->SetLimits(0.,11.);
   gr4->SetMinimum(0.);
   gr4->SetMaximum(0.3);
   gr4->Draw("AL");

   
   TGraph *gr5 = new TGraph(tTotal,aTime,eBy3AuAuP);
   gr5->SetLineColor(2);
   gr5->SetLineWidth(2);
   gr5->SetMarkerColor(1);
   gr5->SetMarkerStyle(22);
   gr5->SetMarkerSize(1.5);
   gr5->Draw("L");
   
   
    TLegend *leg3 = new TLegend(0.60,0.50,0.80,0.70);
    leg3->SetTextFont(62);
    leg3->SetTextSize(0.053);
    leg3->SetLineColor(0);
    leg3->SetLineStyle(0);
    leg3->SetLineWidth(2);
    leg3->SetFillColor(0);
    leg3->SetFillStyle(1012);
    leg3->AddEntry(gr5,"Participants","l");
    leg3->AddEntry(gr4,"Spectators","l");
    leg3->Draw();
    TLegend *legHead3 = new TLegend(0.3,0.80,0.50,0.85);
    legHead3->SetTextFont(62);
    legHead3->SetTextSize(0.055);
    legHead3->SetLineColor(0);
    legHead3->SetLineStyle(0);
    legHead3->SetLineWidth(2);
    legHead3->SetFillColor(0);
    legHead3->SetFillStyle(1012);
    legHead3->SetHeader("c)  b=12-16[fm]","L");
    legHead3->Draw(); 
         
    
     c1->cd(0);
     gStyle->SetOptTitle(0);
     TPaveLabel *title = new TPaveLabel(.2,.93,.8,.98,"Au+Au #sqrt{S_{NN}}=4GeV","brndc"); title->Draw();
     title->SetFillColor(0);
     title->SetTextSize(1.3);
     title->SetTextFont(62); 
     title->SetBorderSize(0);
     c1->Update();
   c1->SaveAs("fig3P-4GeV.pdf");
   c1->SaveAs("fig3P-4GeV.png");


  }
