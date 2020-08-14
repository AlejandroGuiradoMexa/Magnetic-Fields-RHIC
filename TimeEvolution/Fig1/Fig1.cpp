/*void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);

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

void Fig1(){
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

    double eBx1[tTotal] = {}, eBx2[tTotal] = {};
    double eBy1[tTotal] = {},eBy2[tTotal] = {};
    double eBz1[tTotal] = {},eBz2[tTotal] = {};
    double eB1[tTotal] = {}, eB2[tTotal] = {};
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
            for( int iTime = 1; iTime<tTotal; ++iTime){
                Float_t iTC = iTime*0.1;
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
//              if(particle.time  == 0.4){
//              }
        Float_t R = particle.R;
        Float_t B =(particle.eBx)*(particle.eBx)+(particle.eBy)*(particle.eBy)+(particle.eBz)*(particle.eBz);
        if(particle.numbercoll==0 && particle.charge==1 && R>0.3){
//            if(particle.time ==  8.2){
//                   cout<< particle.time  <<"   " <<particle.eBy<<   "   "<< particle.X<<endl;
//              }
            for( int iTime = 1; iTime<tTotal; ++iTime){
                float_t iTC = iTime*0.1;
            //cout << particle.time << "   " << iTC<<endl;

                if(iTC == particle.time){
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
       cout << "  "<<  eBx1[iTime] <<"  "<<  eBy1[iTime] <<"  "<<  eBz1[iTime] <<"  "<<  eB1[iTime] <<"  "<<  iTime*tSteps <<endl;
    }
    for( Int_t iTime = 0; iTime<tTotal; ++iTime){
       cout << "  "<<  eBx2[iTime] <<"  "<<  eBy2[iTime] <<"  "<<  eBz2[iTime] <<"  "<<  eB2[iTime] <<"  "<<  iTime*tSteps <<endl;
    }

    /*
   gStyle->SetOptStat(0);
   TCanvas *C = (TCanvas*) gROOT->FindObject("L");
   if (C) delete C;
   C = new TCanvas("L","canvas",1024,640);
   C->SetFillStyle(4000);    
    gStyle->SetOptStat(false);
    c1->SetRightMargin(0.0465116);
    c1->SetTopMargin(0.1);
    c1->SetFillColor(0);
    c1->Divide(2,1,0.,0.);
    c1->cd(1);
 */

 /*
   // Number of PADS
   const Int_t Nx = 2;
   const Int_t Ny = 1;
   // Margins
   Float_t lMargin = 0.12;
   Float_t rMargin = 0.05;
   Float_t bMargin = 0.15;
   Float_t tMargin = 0.05;
   // Canvas setup
   CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

   // Dummy histogram.

   TPad *pad[Nx][Ny];
  TGraph *gr = new TGraph(tTotal,aTime,eBy1); 

  TH1F *h1 = gr->GetHistogram();
   h1->GetXaxis()->SetTitle("Time[fm]");
   h1->GetYaxis()->SetTitle("B/m_{#pi}^{2}");
   for (Int_t i=0;i<Nx;i++) {
      for (Int_t j=0;j<Ny;j++) {
         C->cd(0);
         // Get the pads previously created.
         char pname[16];
         sprintf(pname,"pad_%i_%i",i,j);
         pad[i][j] = (TPad*) gROOT->FindObject(pname);
         pad[i][j]->Draw();
         pad[i][j]->SetFillStyle(4000);
         pad[i][j]->SetFrameFillStyle(4000);
         pad[i][j]->cd();
         // Size factors
         Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
         Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();
         char hname[16];
         sprintf(hname,"h_%i_%i",i,j);
         TH1F *hFrame = (TH1F*) h1->Clone(hname);
         hFrame->Reset();
         hFrame->Draw();
         // y axis range
         hFrame->GetYaxis()->SetRangeUser(0.0001,1.2*h1->GetMaximum());
         // Format for y axis
         hFrame->GetYaxis()->SetLabelFont(43);
         hFrame->GetYaxis()->SetLabelSize(16);
         hFrame->GetYaxis()->SetLabelOffset(0.02);
         hFrame->GetYaxis()->SetTitleFont(43);
         hFrame->GetYaxis()->SetTitleSize(16);
         hFrame->GetYaxis()->SetTitleOffset(5);
         hFrame->GetYaxis()->CenterTitle();
         hFrame->GetYaxis()->SetNdivisions(505);
         // TICKS Y Axis
         hFrame->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
         // Format for x axis
         hFrame->GetXaxis()->SetLabelFont(43);
         hFrame->GetXaxis()->SetLabelSize(16);
         hFrame->GetXaxis()->SetLabelOffset(0.02);
         hFrame->GetXaxis()->SetTitleFont(43);
         hFrame->GetXaxis()->SetTitleSize(16);
         hFrame->GetXaxis()->SetTitleOffset(5);
         hFrame->GetXaxis()->CenterTitle();
         hFrame->GetXaxis()->SetNdivisions(505);
         // TICKS X Axis
         hFrame->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
         gr->Draw("ACP");
      }
   }
   C->cd();
   */
  //     /*
   TCanvas *c1 =  new TCanvas("c1","",1024,640);
      c1->SetFillStyle(4000);    
    gStyle->SetOptStat(false);
    c1->SetRightMargin(0.0465116);
    c1->SetTopMargin(0.1);
    c1->SetFillColor(0);
    c1->Divide(2,1);
    c1->cd(1);
    TGraph *gr = new TGraph(tTotal,aTime,eBx1); 
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
   gr->Draw("AL");

   TGraph *gr1 = new TGraph(tTotal,aTime,eBy1);
   gr1->SetLineColor(2);
   gr1->SetLineWidth(2);
   gr1->SetMarkerColor(1);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.5);
   gr1->Draw("L");
   
   TGraph *gr2 = new TGraph(tTotal,aTime,eBz1);
   gr2->SetLineColor(3);
   gr2->SetLineWidth(2);
   gr2->SetMarkerColor(1);
   gr2->SetMarkerStyle(23);
   gr2->SetMarkerSize(1.5);
   gr2->Draw("L");
 
   TGraph *gr3 = new TGraph(tTotal,aTime,eB1);
   gr3->SetLineColor(4);
   gr3->SetLineWidth(2);
   gr3->SetMarkerColor(1);
   gr3->SetMarkerStyle(33);
   gr3->SetMarkerSize(1.5);
   gr3->Draw("L");
   
    TLegend *leg1 = new TLegend(0.60,0.60,0.80,0.85);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.045);
    leg1->SetLineColor(0);
    leg1->SetLineStyle(0);
    leg1->SetLineWidth(2);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1012);
    leg1->SetHeader("Magnetic Field","L");
    leg1->AddEntry(gr,"eB_{x}","L");
    leg1->AddEntry(gr1,"-eB_{y}","L");
    leg1->AddEntry(gr2,"eB_{z}"  ,"L");
    leg1->AddEntry(gr3,"eB","L");
    leg1->Draw();

    TLegend *leg2 = new TLegend(0.45,0.42,0.65,0.54);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.055);
    leg2->SetLineColor(0);
    leg2->SetLineStyle(0);
    leg2->SetLineWidth(2);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1012);
    leg2->SetHeader("BiBi, #sqrt{S_{NN}}=4GeV ","L");
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
   grA->Draw("AL");
   
   TGraph *gr1A = new TGraph(tTotal,aTime,eBy2);
   gr1A->SetLineColor(2);
   gr1A->SetLineWidth(2);
   gr1A->SetMarkerColor(1);
   gr1A->SetMarkerStyle(22);
   gr1A->SetMarkerSize(1.5);
   gr1A->Draw("L");
   
   TGraph *gr2A = new TGraph(tTotal,aTime,eBz2);
   gr2A->SetLineColor(3);
   gr2A->SetLineWidth(2);
   gr2A->SetMarkerColor(1);
   gr2A->SetMarkerStyle(23);
   gr2A->SetMarkerSize(1.5);
   gr2A->Draw("L");
 
   TGraph *gr3A = new TGraph(tTotal,aTime,eB2);
   gr3A->SetLineColor(4);
   gr3A->SetLineWidth(2);
   gr3A->SetMarkerColor(1);
   gr3A->SetMarkerStyle(33);
   gr3A->SetMarkerSize(1.5);
   gr3A->Draw("L");
   
    TLegend *leg1A = new TLegend(0.60,0.60,0.80,0.85);
    leg1A->SetTextFont(42);
    leg1A->SetTextSize(0.045);
    leg1A->SetLineColor(0);
    leg1A->SetLineStyle(0);
    leg1A->SetLineWidth(2);
    leg1A->SetFillColor(0);
    leg1A->SetFillStyle(1012);
    leg1A->SetHeader("Magnetic Field","L");
    leg1A->AddEntry(grA,"eB_{x}","L");
    leg1A->AddEntry(gr1A,"-eB_{y}","L");
    leg1A->AddEntry(gr2A,"eB_{z}"  ,"L");
    leg1A->AddEntry(gr3A,"eB","L");
    leg1A->Draw();

    TLegend *leg2A = new TLegend(0.45,0.42,0.65,0.54);
    leg2A->SetTextFont(42);
    leg2A->SetTextSize(0.055);
    leg2A->SetLineColor(0);
    leg2A->SetLineStyle(0);
    leg2A->SetLineWidth(2);
    leg2A->SetFillColor(0);
    leg2A->SetFillStyle(1012);
    leg2A->SetHeader("AuAu,#sqrt{S_{NN}}=4GeV ","L");
    leg2A->Draw(); 


   c1->SaveAs("Fig1-4GeV.pdf");
   c1->SaveAs("Fig1-4GeV.png");
  delete c1;
  delete gr,gr1,gr2,gr3,leg1,leg2;
  delete grA,gr1A,gr2A,gr3A,leg1A,leg2A;

  //*/

}
/*
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin){
   if (!C) return;
   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;
   for (Int_t i=0;i<Nx;i++) {
      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }
      for (Int_t j=0;j<Ny;j++) {
         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }
         C->cd(0);
         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);
         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);
         pad->Draw();
      }
   }

}
*/
