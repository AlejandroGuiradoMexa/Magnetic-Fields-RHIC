/*
This macro read the test.f14 file from urqmd_3.4 after used the urqmdtoroot.C macro 
to create a ROOT file. The magnetic fields are calculated in the particle point aproximation
using the Lienard-Wiechert potentials. 

Author: Alejandro Guirado.
Date: 25/03/20.
*/

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

void MagneticField_TemporalEvolution(){
 StylePlots();

  const int tI = 0.0;
  const int tF = 10.0;
  const float tSteps = 0.5;
  const int tTotal = 21;
  const int numB = 5.;
  const int protonNumber1 = 166.0;
  const int protonNumber2 = 158.0;

  gROOT->Reset();
  TChain mychainBiBi("T");
    mychainBiBi.Add("../Datos/BiBi100_11GeV.root");
  TChain mychainAuAu("T");
    mychainAuAu.Add("../Datos/AuAu100_11GeV.root");
  TChain mychainBiBi0("T");
    mychainBiBi0.Add("../Datos/BiBi1000_11GeV_t0.root");
  TChain mychainAuAu0("T");
    mychainAuAu0.Add("../Datos/AuAu1000_11GeV_t0.root");
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

    double_t eBx1[tTotal] = {}, eBx2[tTotal] = {};
    double_t eBy1[tTotal] = {},eBy2[tTotal] = {};
    double_t eBz1[tTotal] = {},eBz2[tTotal] = {};
    double_t eB1[tTotal] = {}, eB2[tTotal] = {};
    double_t  eBx1_0 = 0., eBy1_0 = 0., eBz1_0 = 0., eB1_0 = 0.;
    double_t  eBx2_0 = 0., eBy2_0 = 0., eBz2_0 = 0., eB2_0 = 0.;

    double_t aTime[tTotal] = {0.};    
    double_t particlesTotal1 = 0.,particlesTotal2 = 0.;
    double_t particlesTotal1_0 = 0.,particlesTotal2_0 = 0.;




    for(Int_t j=0; j<nevent1; j++){
      mychainBiBi.GetEvent(j);
        if(particle.time == tSteps && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            particlesTotal1++;
        }
    }
    for(Int_t j=0; j<nevent2; j++){
      mychainAuAu.GetEvent(j);
        if(particle.time == tSteps && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            particlesTotal2++;
        }
    }


    Double_t norm1=(particlesTotal1*137)/protonNumber1;
    Double_t norm2=(particlesTotal2*137)/protonNumber2;
    Double_t norm1_0= 1000.*137.; 
    Double_t norm2_0= 1000.*137.; 
    

    for (Int_t iEv=0; iEv<nevent1_0; iEv++) 
    {
        mychainBiBi0.GetEvent(iEv);
        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){
                    eBx1_0 += particle.eBx;
                    eBy1_0 += particle.eBy;
                    eBz1_0 += particle.eBz;
                    eB1_0 += TMath::Sqrt(B); 
        }
    }    

    
    for (Int_t iEv=0; iEv<nevent2_0; iEv++) 
    {
        mychainAuAu0.GetEvent(iEv);
        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){
                    eBx2_0 += particle.eBx;
                    eBy2_0 += particle.eBy;
                    eBz2_0 += particle.eBz;
                    eB2_0 += TMath::Sqrt(B); 
        }
    }    


    for (Int_t iEv=0; iEv<nevent1; iEv++) 
    {
        mychainBiBi.GetEvent(iEv);

        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){
            for( Int_t iTime = 1; iTime<tTotal; ++iTime){
                double_t iTC = iTime*tSteps;
                if(particle.time == iTC){
                    eBx1[iTime] += particle.eBx;
                    eBy1[iTime] += particle.eBy;
                    eBz1[iTime] += particle.eBz;
                    eB1[iTime] += TMath::Sqrt(B);                
                }
            }
        }
    }     

    for (Int_t iEv=0; iEv<nevent2; iEv++) {
        mychainAuAu.GetEvent(iEv);
        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){
            for( Int_t iTime = 1; iTime<tTotal; ++iTime){
                double_t iTC = iTime*tSteps;
                if(particle.time == iTC){
                    eBx2[iTime] += particle.eBx;
                    eBy2[iTime] += particle.eBy;
                    eBz2[iTime] += particle.eBz;
                    eB2[iTime] += TMath::Sqrt(B);                
                }
            }
        }
    }   
  

    eBx1[0] = eBx1_0/norm1_0;
    eBy1[0] = eBy1_0/norm1_0;
    eBz1[0] = eBz1_0/norm1_0;
    eB1[0] = eB1_0/norm1_0;

    
    for( Int_t iTime = 1; iTime<tTotal; ++iTime){
        eBx1[iTime] =eBx1[iTime]/norm1;
        eBy1[iTime] =eBy1[iTime]/norm1;
        eBz1[iTime] =eBz1[iTime]/norm1;
        eB1[iTime] =eB1[iTime]/norm1;
    }

    eBx2[0] = eBx2_0/norm2_0;
    eBy2[0] = eBy2_0/norm2_0;
    eBz2[0] = eBz2_0/norm2_0;
    eB2[0] = eB2_0/norm2_0;

    for( Int_t iTime = 1; iTime<tTotal; ++iTime){
        eBx2[iTime] =eBx2[iTime]/norm2;
        eBy2[iTime] =eBy2[iTime]/norm2;
        eBz2[iTime] =eBz2[iTime]/norm2;
        eB2[iTime] =eB2[iTime]/norm2;
    }

    for( Int_t iTime = 0; iTime<tTotal; ++iTime){
        double_t iTC = iTime*tSteps;
        aTime[iTime] = iTC;
    }

      for( Int_t iTime = 0; iTime<tTotal; ++iTime){
       cout << "  "<<  eBx1[iTime] <<"  "<<  eBy1[iTime] <<"  "<<  eBz1[iTime] <<"  "<<  eB1[iTime] <<"  "<<  iTime <<endl;
    }
    for( Int_t iTime = 0; iTime<tTotal; ++iTime){
       cout << "  "<<  eBx2[iTime] <<"  "<<  eBy2[iTime] <<"  "<<  eBz2[iTime] <<"  "<<  eB2[iTime] <<"  "<<  iTime <<endl;
    }

    
   gStyle->SetOptStat(0);
   TCanvas *C = (TCanvas*) gROOT->FindObject("C");
   if (C) delete C;
   C = new TCanvas("C","canvas",1024,640);
   C->SetFillStyle(4000);    
    gStyle->SetOptStat(false);
    c1->SetRightMargin(0.0465116);
    c1->SetTopMargin(0.1);
    c1->SetFillColor(0);
    c1->Divide(2,1,0.,0.);
    c1->cd(1);

   TCanvas *c1 =  new TCanvas("c1","",1024,640);
      c1->SetFillStyle(4000);    
    gStyle->SetOptStat(false);
    c1->SetRightMargin(0.0465116);
    c1->SetTopMargin(0.1);
    c1->SetFillColor(0);
    c1->Divide(2,1);
    c1->cd(1);
    TGraph *gr = new TGraph(tTotal+1,aTime,eBx1); 
   gr->SetLineColor(1);
   gr->SetLineWidth(2);
   gr->SetMarkerColor(1);
   gr->SetMarkerStyle(21);
   gr->SetMarkerSize(1.5);
   gr->SetTitle("");
   gr->GetXaxis()->SetTitle("Time[fm]");
   gr->GetYaxis()->SetTitle("B/m_{#pi}^{2}");
   gr->GetXaxis()->SetLimits(0.,11.);
   gr->SetMinimum(0.);
   gr->SetMaximum(0.35);
   gr->Draw("ACP");

   TGraph *gr1 = new TGraph(tTotal+1,aTime,eBy1);
   gr1->SetLineColor(2);
   gr1->SetLineWidth(2);
   gr1->SetMarkerColor(1);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.5);
   gr1->Draw("CP");
   
   TGraph *gr2 = new TGraph(tTotal+1,aTime,eBz1);
   gr2->SetLineColor(3);
   gr2->SetLineWidth(2);
   gr2->SetMarkerColor(1);
   gr2->SetMarkerStyle(23);
   gr2->SetMarkerSize(1.5);
   gr2->Draw("CP");
 
   TGraph *gr3 = new TGraph(tTotal+1,aTime,eB1);
   gr3->SetLineColor(4);
   gr3->SetLineWidth(2);
   gr3->SetMarkerColor(1);
   gr3->SetMarkerStyle(33);
   gr3->SetMarkerSize(1.5);
   gr3->Draw("CP");
   
    TLegend *leg1 = new TLegend(0.60,0.60,0.80,0.85);
    leg1->SetTextFont(62);
    leg1->SetTextSize(0.045);
    leg1->SetLineColor(0);
    leg1->SetLineStyle(0);
    leg1->SetLineWidth(2);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1012);
    leg1->SetHeader("Magnetic Field","C");
    leg1->AddEntry(gr,"eB_{x}","lp");
    leg1->AddEntry(gr1,"-eB_{y}","lp");
    leg1->AddEntry(gr2,"eB_{z}"  ,"lp");
    leg1->AddEntry(gr3,"eB","lp");
    leg1->Draw();

    TLegend *leg2 = new TLegend(0.60,0.20,0.80,0.4);
    leg2->SetTextFont(62);
    leg2->SetTextSize(0.055);
    leg2->SetLineColor(0);
    leg2->SetLineStyle(0);
    leg2->SetLineWidth(2);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1012);
    leg2->SetHeader("BiBi, #sqrt{S_{NN}}=11GeV ","C");
    leg2->Draw(); 

    c1->cd(2);

   TGraph *grA = new TGraph(tTotal,aTime,eBx2);
   grA->SetLineColor(1);
   grA->SetLineWidth(2);
   grA->SetMarkerColor(1);
   grA->SetMarkerStyle(21);
   grA->SetMarkerSize(1.5);
   grA->SetTitle("");
   grA->GetXaxis()->SetTitle("Time[fm]");
   grA->GetYaxis()->SetTitle("B/m_{#pi}^{2}");
   grA->GetXaxis()->SetLimits(0.,11.);
   grA->SetMinimum(0.);
   grA->SetMaximum(0.35);
   grA->Draw("ACP");
   
   TGraph *gr1A = new TGraph(tTotal,aTime,eBy2);
   gr1A->SetLineColor(2);
   gr1A->SetLineWidth(2);
   gr1A->SetMarkerColor(1);
   gr1A->SetMarkerStyle(22);
   gr1A->SetMarkerSize(1.5);
   gr1A->Draw("CP");
   
   TGraph *gr2A = new TGraph(tTotal,aTime,eBz2);
   gr2A->SetLineColor(3);
   gr2A->SetLineWidth(2);
   gr2A->SetMarkerColor(1);
   gr2A->SetMarkerStyle(23);
   gr2A->SetMarkerSize(1.5);
   gr2A->Draw("CP");
 
   TGraph *gr3A = new TGraph(tTotal,aTime,eB2);
   gr3A->SetLineColor(4);
   gr3A->SetLineWidth(2);
   gr3A->SetMarkerColor(1);
   gr3A->SetMarkerStyle(33);
   gr3A->SetMarkerSize(1.5);
   gr3A->Draw("CP");
   
    TLegend *leg1A = new TLegend(0.60,0.60,0.80,0.85);
    leg1A->SetTextFont(62);
    leg1A->SetTextSize(0.045);
    leg1A->SetLineColor(0);
    leg1A->SetLineStyle(0);
    leg1A->SetLineWidth(2);
    leg1A->SetFillColor(0);
    leg1A->SetFillStyle(1012);
    leg1A->SetHeader("Magnetic Field","C");
    leg1A->AddEntry(grA,"eB_{x}","lp");
    leg1A->AddEntry(gr1A,"-eB_{y}","lp");
    leg1A->AddEntry(gr2A,"eB_{z}"  ,"lp");
    leg1A->AddEntry(gr3A,"eB","lp");
    leg1A->Draw();

    TLegend *leg2A = new TLegend(0.60,0.20,0.80,0.4);
    leg2A->SetTextFont(62);
    leg2A->SetTextSize(0.055);
    leg2A->SetLineColor(0);
    leg2A->SetLineStyle(0);
    leg2A->SetLineWidth(2);
    leg2A->SetFillColor(0);
    leg2A->SetFillStyle(1012);
    leg2A->SetHeader("AuAu,#sqrt{S_{NN}}=11GeV ","C");
    leg2A->Draw(); 


   c1->SaveAs("fig2.pdf");
   c1->SaveAs("fig2.png");
  delete c1;
  delete gr,gr1,gr2,gr3,leg1,leg2;
  delete grA,gr1A,gr2A,gr3A,leg1A,leg2A;


}