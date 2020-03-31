/*
This macro read the test.f14 file from urqmd_3.4 after used the urqmdtoroot.C macro 
to create a ROOT file. The magnetic fields are calculated in the particle point aproximation
using the Lienard-Wiechert potentials.  It calculates the magnetic fields in y-direction which
is normal to reaction plane formed by the impact parameter and beam direction for each centrality bin
determined by the total number of particles and the number of spectators.


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

void MagneticField_y_CentralityBin(){
StylePlots();

  const int tI = 0.0;
  const int tF = 10.0;
  const float tSteps = 0.5;
  const int tTotal = 20;
  const int numB = 5.;
  const int protonNumber = 158.;

  gROOT->Reset();
  //for(Int_t j=1; j<=161; j++){
  TChain mychain("T");
    mychain.Add("../Datos/AuAu100_11GeV.root");
    Int_t nlines = 0;
    struct particula_t 
    {
      Float_t time,X,Y,Z,E,Px,Py,Pz,Pt,P,m,id,isoespin,charge,lastcoll,numbercoll,history,frezetime,frezeX,frezeY,frezeZ,frezeE,frezePx,frezePy,frezePz,b,nspec,R,PXR,eBx,eBy,eBz;
    } PARTICLE;
    
    particula_t  particle;
    mychain.SetBranchAddress("particle",&particle);
    Int_t nevent = mychain.GetEntries();
    Double_t c20=0, c40=0, c60=0, c80=0, c100=0;
    double_t eBx[tTotal][numB] = {0.};
    double_t eBy[tTotal][numB] = {0.};
    double_t eBz[tTotal][numB] = {0.};
    double_t eB[tTotal][numB] = {0.};
    double_t aTime[tTotal] = {0.};    
    double_t impacPar[numB+1] = {0,5,8,11,14,16};    

    double_t particlesTotal  = 0.;
    for(Int_t j=0; j<nevent; j++){
      mychain.GetEvent(j);
        if(particle.time == tSteps && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            particlesTotal++;
            if(particle.b>=0 && particle.b<=5.){
            c20++;
            }
            if(particle.b>=5.7 && particle.b<=8.4){
            c40++;
            }
            if(particle.b>=8.4 && particle.b<=10.1){
            c60++;
            }
            if(particle.b>=10.1 && particle.b<=11.7){
            c80++;
            }
            if(particle.b>=11.7 && particle.b<=15.7){
            c100++;
            }
        }
    }
    
    
    Double_t norm20=(c20*137)/protonNumber;
    Double_t norm40=(c40*137)/protonNumber;
    Double_t norm60=(c60*137)/protonNumber;
    Double_t norm80=(c80*137)/protonNumber;
    Double_t norm100=(c100*137)/protonNumber;
    
    cout << norm20 << "  " << norm40 << "  " << norm60 << "  " << norm80 << "  " << norm100 << "  " <<endl;
    
    
    for (Int_t iEv=0; iEv<nevent; iEv++) 
    {
        mychain.GetEvent(iEv);
        Double_t x = 0;
        Double_t R = particle.R;
        double_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){
            for( int iB = 0.;iB<numB;++iB){
                if(particle.b>=impacPar[iB] && particle.b<=impacPar[iB+1]){
                    for( Int_t iTime = tI; iTime<=tTotal; ++iTime){
                        double_t iTC = iTime*tSteps;
                        if(particle.time == iTC){
                            eBx[iTime][iB] += particle.eBx;
                            eBy[iTime][iB] += particle.eBy;
                            eBz[iTime][iB] += particle.eBz;
                            eB[iTime][iB] += TMath::Sqrt(B);
                          //  cout << iTime << "  "<< particle.eBx<<   "  " << eBx[iTime][iB] <<  "  " << iB <<endl;
                        }
                    }  
                }
            }
        }
    }        
     
     for( int iB = 0.;iB<numB;++iB){
        for( Int_t iTime = tI; iTime<=tTotal; ++iTime){
            cout<< eBx[iTime][iB] <<  "  "<< eBy[iTime][iB] <<  "  "<< eBz[iTime][iB]<<"  " << eB[iTime][iB] <<  "  " <<  "  " << iB << " " << iTime <<endl;
        }
     }
     
    double_t aBxC1[tTotal] = {0},aBxC2[tTotal] = {0},aBxC3[tTotal] = {0},aBxC4[tTotal] = {0},aBxC5[tTotal] = {0}; 
    for( Int_t iTime = tI; iTime<=tTotal; ++iTime){
        double_t iTC = iTime*tSteps;
        aTime[iTime] = iTC;
        aBxC1[iTime] = eBy[iTime][0]/(norm20);
        aBxC2[iTime] = eBy[iTime][1]/(norm40);
        aBxC3[iTime] = eBy[iTime][2]/(norm60);
        aBxC4[iTime] = eBy[iTime][3]/(norm80);
        aBxC5[iTime] = eBy[iTime][4]/(norm100);

        
        //cout << iTC << "   "<<  aBxC1[iTime] 
//         << "   "<<  aBxC2[iTime] 
//         << "   "<<  aBxC3[iTime] 
//         << "   "<<  aBxC4[iTime] <<endl;
  //  cout <<  aBxC1[iTime] << "    " << eBy[iTime][0]/(norm20) << "  "<<eBy[iTime][0]<<endl;
        
    }

    
   TCanvas* c1 = new TCanvas("c1","UrQMD test example",800,800);
    gStyle->SetOptStat(false);
    c1->SetRightMargin(0.0465116);
    c1->SetTopMargin(0.1);
    c1->SetFillColor(0);
    
    
   TGraph *gr = new TGraph(tTotal,aTime,aBxC1);
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
   gr->SetMaximum(0.5);
   gr->Draw("ACP");
   
   TGraph *gr1 = new TGraph(tTotal,aTime,aBxC2);
   gr1->SetLineColor(2);
   gr1->SetLineWidth(2);
   gr1->SetMarkerColor(1);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.5);
   gr1->Draw("CP");
   
   TGraph *gr2 = new TGraph(tTotal,aTime,aBxC3);
   gr2->SetLineColor(3);
   gr2->SetLineWidth(2);
   gr2->SetMarkerColor(1);
   gr2->SetMarkerStyle(23);
   gr2->SetMarkerSize(1.5);
   gr2->Draw("CP");
 
   TGraph *gr3 = new TGraph(tTotal,aTime,aBxC4);
   gr3->SetLineColor(4);
   gr3->SetLineWidth(2);
   gr3->SetMarkerColor(1);
   gr3->SetMarkerStyle(33);
   gr3->SetMarkerSize(1.5);
   gr3->Draw("CP");

   TGraph *gr4 = new TGraph(tTotal,aTime,aBxC5);
   gr4->SetLineColor(5);
   gr4->SetLineWidth(2);
   gr4->SetMarkerColor(1);
   gr4->SetMarkerStyle(29);
   gr4->SetMarkerSize(1.5);
   gr4->Draw("CP");
   
    TLegend *leg1 = new TLegend(0.60,0.60,0.80,0.85);
    leg1->SetTextFont(62);
    leg1->SetTextSize(0.04);
    leg1->SetLineColor(0);
    leg1->SetLineStyle(0);
    leg1->SetLineWidth(2);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1012);
    leg1->SetHeader("Centrality","C");
    leg1->AddEntry(gr,"0-20%","lp");
    leg1->AddEntry(gr1,"20-40%","lp");
    leg1->AddEntry(gr2,"40-60%"  ,"lp");
    leg1->AddEntry(gr3,"60-80%","lp");
    leg1->AddEntry(gr4,"80-100%","lp");
    leg1->Draw();

    TLegend *leg2 = new TLegend(0.60,0.20,0.80,0.40);
    leg2->SetTextFont(62);
    leg2->SetTextSize(0.045);
    leg2->SetLineColor(0);
    leg2->SetLineStyle(0);
    leg2->SetLineWidth(2);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1012);
    leg2->SetHeader("#sqrt{S_{NN}}=9GeV ","C");
    leg2->Draw(); 
   
   
   c1->SaveAs("9GeV_Centrality_By.pdf");
   c1->SaveAs("9GeV_Centrality_By.png");

}

