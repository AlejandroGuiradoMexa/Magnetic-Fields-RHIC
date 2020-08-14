
    double MagneticFieldX(
        int x0, int y0, double X, double Y, double Z,
        double Px,double Py, double Pz, double E){

	    Float_t  Pt= -999.0 ,P= -999.0;

        Pt = TMath::Sqrt(Px*Px+ Py*Py);
	    P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
	    
        double R = 0., PXR= 0.;
        double eBx = 0., eBy = 0., eBz = 0.;
        X = x0-X;
	    R= TMath::Sqrt(X*X+Y*Y+Z*Z);
	    PXR= TMath::Sqrt((Py*Z-Pz*Y)*(Py*Z-Pz*Y)+(X*Pz-Px*Z)*(X*Pz-Px*Z)+(Px*Y-Py*X)*(Px*Y-Py*X));
	    if(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3)!=0){
        eBx = (100/49)*(E*E-P*P)*(Py*Z-Pz*Y)/(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3));
        eBy = (100/49)*(E*E-P*P)*(Pz*X-Px*Z)/(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3));
        eBz = (100/49)*(E*E-P*P)*(Px*Y-Py*X)/(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3));
    }
    return eBy;
}

       double MagneticFieldY(
        int x0, int y0, double X, double Y, double Z,
        double Px,double Py, double Pz, double E){

	    Float_t  Pt= -999.0 ,P= -999.0;

        Pt = TMath::Sqrt(Px*Px+ Py*Py);
	    P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
	    
        double R = 0., PXR= 0.;
        double eBx = 0., eBy = 0., eBz = 0.;
        Y = y0-Y;
	    R= TMath::Sqrt(X*X+Y*Y+Z*Z);
	    PXR= TMath::Sqrt((Py*Z-Pz*Y)*(Py*Z-Pz*Y)+(X*Pz-Px*Z)*(X*Pz-Px*Z)+(Px*Y-Py*X)*(Px*Y-Py*X));
	    if(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3)!=0){
        eBx = (100/49)*(E*E-P*P)*(Py*Z-Pz*Y)/(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3));
        eBy = (100/49)*(E*E-P*P)*(Pz*X-Px*Z)/(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3));
        eBz = (100/49)*(E*E-P*P)*(Px*Y-Py*X)/(pow(TMath::Sqrt((E*R)*(E*R)-(PXR)*(PXR)),3));
    }
    return eBy;
}

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

void Fig4(){
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
    
    mychainAuAu.SetBranchAddress("particle",&particle);
    mychainAuAu.SetBranchAddress("particle",&particle);
    mychainAuAu0.SetBranchAddress("particle",&particle);
    mychainAuAu0.SetBranchAddress("particle",&particle);

    Int_t nevent1 = mychainAuAu.GetEntries();
    Int_t nevent2 = mychainAuAu.GetEntries();
    Int_t nevent1_0 = mychainAuAu0.GetEntries();
    Int_t nevent2_0 = mychainAuAu0.GetEntries();
    
    const int NumberOfPositions = 6.;
    double B_X0[NumberOfPositions] = {0.}, B_Y0[NumberOfPositions] = {0.};
    double_t B_X[tTotal][NumberOfPositions] = {};
    double_t B_Y[tTotal][NumberOfPositions] = {};

    double_t eBy1[tTotal] = {},eBy2[tTotal] = {};

    //Graphs
    TGraph *gr_eBy_X[6] = {NULL};
    TGraph *gr_eBy_Y[6] = {NULL};

    
    for(Int_t j=0; j<nevent1_0; j++){
      mychainAuAu0.GetEvent(j);
      double_t refX = 0.,refY = 0.;
      double X = particle.X, Y = particle.Y, Z = particle.Z;
      double Px = particle.Px, Py = particle.Py, Pz = particle.Pz;
      double E = particle.E; 
      double R = particle.R;
      if(particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
        for (int iPosition = 0.; iPosition<NumberOfPositions; iPosition++){
            refX = MagneticFieldX( iPosition, 0., X,Y,Z,Px,Py,Pz,E );
            refY = MagneticFieldY( 0.,iPosition, X,Y,Z,Px,Py,Pz,E );
            B_X0[iPosition] += refX;
            B_Y0[iPosition] += refY;
        }
      }
    }
    
    
    for(Int_t j=0; j<nevent1; j++){
      mychainAuAu.GetEvent(j);
      double_t refX = 0.,refY = 0.;
      double X = particle.X, Y = particle.Y, Z = particle.Z;
      double Px = particle.Px, Py = particle.Py, Pz = particle.Pz;
      double E = particle.E; 
      double R = particle.R;

      if(particle.numbercoll==0 && particle.charge==1 && particle.id == 1. && R>0.3){
          for( Int_t iTime = 1; iTime<tTotal; ++iTime){
                Float_t iTC = iTime*0.1;
                if(particle.time == iTC){
                    for (int iPosition = 0.; iPosition<NumberOfPositions; iPosition++){
                        refX = MagneticFieldX( iPosition, 0., X,Y,Z,Px,Py,Pz,E );
                        refY = MagneticFieldY( 0.,iPosition, X,Y,Z,Px,Py,Pz,E );
                            B_X[iTime][iPosition] += refX;
                            B_Y[iTime][iPosition] += refY;
                    }   
                }
            }
        }
    }

    //Constants for MF
    const int nevAuAu0 = 1000.;
    const int nevAuAu = 1000.;
    const double normAuAu0 = nevAuAu0*137.;
    const double normAuAu = nevAuAu*137.;   
    
//     const int nevAuAu0 = 1000.;
//     const int nevAuAu = 100.;
//     const double normAuAu0 = nevAuAu0*137.;
//     const double normAuAu = nevAuAu*137.;   

    
    
    for (int iPosition = 0.; iPosition<NumberOfPositions; iPosition++){
        for( Int_t iTime = 1; iTime<tTotal; ++iTime){
            B_X[iTime][iPosition] = B_X[iTime][iPosition]/normAuAu;
            B_Y[iTime][iPosition] = B_Y[iTime][iPosition]/normAuAu;
            cout << B_Y[iTime][iPosition] << "  " << iTime << "   " << iPosition <<endl;
        }
    }
    for (int iPosition = 0.; iPosition<NumberOfPositions; iPosition++){
        B_X0[iPosition] = B_X0[iPosition]/normAuAu0;    
        B_Y0[iPosition] = B_Y0[iPosition]/normAuAu0;
        cout << B_X0[iPosition]<< 
          "   " << B_Y0[iPosition]<<endl;     
        B_X[0][iPosition] = B_X0[iPosition];
        B_Y[0][iPosition] = B_Y0[iPosition];
    }    
    for (int iPosition = 0.; iPosition<NumberOfPositions; iPosition++){
        for( Int_t iTime = 0; iTime<tTotal; ++iTime){
                cout << B_X[iTime][iPosition] << "  " <<  "   " << B_Y[iTime][iPosition] <<"   " <<iTime << "   " << iPosition <<endl;
        }
    }
    for (int iPosition = 0.; iPosition<NumberOfPositions; iPosition++){
        gr_eBy_X[iPosition] = new TGraph();
        gr_eBy_Y[iPosition] = new TGraph();
        for( Int_t iTime = 0; iTime<tTotal; ++iTime){
            Float_t iTC = iTime*0.1;
            gr_eBy_X[iPosition]->SetPoint(iTime,iTC,-B_X[iTime][iPosition]);
            gr_eBy_Y[iPosition]->SetPoint(iTime,iTC,B_Y[iTime][iPosition]);
            cout << iTime << "   " << B_Y[iTime][0]<<endl;
        }
    }
        
        
    //Defining Graph settings
    
    TCanvas *c1 =  new TCanvas("c1","",1024,640);
    c1->SetFillStyle(4000);    
    gStyle->SetOptStat(false);
    c1->SetRightMargin(0.0465116);
    c1->SetTopMargin(0.1);
    c1->SetFillColor(0);
    c1->Divide(2,1);
    c1->cd(1);
    
   gr_eBy_X[0]->SetLineColor(1);
   gr_eBy_X[0]->SetLineWidth(2);

   gr_eBy_X[0]->SetTitle("");
   gr_eBy_X[0]->GetXaxis()->SetTitle("Time[fm]");
   gr_eBy_X[0]->GetYaxis()->SetTitle("-B_{y}/m_{#pi}^{2}");
   gr_eBy_X[0]->GetXaxis()->SetLimits(0.,11.);
   gr_eBy_X[0]->SetMinimum(0.);
   gr_eBy_X[0]->SetMaximum(0.35);
   gr_eBy_X[0]->Draw("ACP");
   
    for (int iPosition = 1.; iPosition<NumberOfPositions; iPosition++){
        gr_eBy_X[iPosition]->SetLineColor(iPosition+1);
        gr_eBy_X[iPosition]->SetLineWidth(2);
        gr_eBy_X[iPosition]->SetLineStyle(iPosition+1);
        gr_eBy_X[iPosition]->Draw("CP");

    }
    TLegend *leg1 = new TLegend(0.50,0.50,0.90,0.85);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.055);
    leg1->SetLineColor(0);
    leg1->SetLineStyle(0);
    leg1->SetLineWidth(2);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1012);
    leg1->AddEntry(gr_eBy_X[0],"x = 0fm","lp");
    leg1->AddEntry(gr_eBy_X[1],"x = 1fm","lp");
    leg1->AddEntry(gr_eBy_X[2],"x = 2fm","lp");
    leg1->AddEntry(gr_eBy_X[3],"x = 3fm","lp");
    leg1->AddEntry(gr_eBy_X[4],"x = 4fm","lp");
    leg1->AddEntry(gr_eBy_X[5],"x = 5fm","lp");
    leg1->Draw();

    TLegend *leg4 = new TLegend(0.50,0.35,0.90,0.48);
    leg4->SetTextFont(42);
    leg4->SetTextSize(0.055);
    leg4->SetLineColor(0);
    leg4->SetLineStyle(0);
    leg4->SetLineWidth(2);
    leg4->SetFillColor(0);
    leg4->SetFillStyle(1012);
    leg4->SetHeader("AuAu, #sqrt{S_{NN}}=4GeV ","C");
    leg4->Draw();
    
    
    c1->cd(2);
    
   gr_eBy_Y[0]->SetLineColor(1);
   gr_eBy_Y[0]->SetLineWidth(2);
   gr_eBy_Y[0]->SetTitle("");
   gr_eBy_Y[0]->GetXaxis()->SetTitle("Time[fm]");
   gr_eBy_Y[0]->GetYaxis()->SetTitle("-B_{y}/m_{#pi}^{2}");
   gr_eBy_Y[0]->GetXaxis()->SetLimits(0.,11.);
   gr_eBy_Y[0]->SetMinimum(0.);
   gr_eBy_Y[0]->SetMaximum(0.35);
   gr_eBy_Y[0]->Draw("ACP");
   
    for (int iPosition = 1.; iPosition<NumberOfPositions; iPosition++){
        gr_eBy_Y[iPosition]->SetLineColor(iPosition+1);
        gr_eBy_Y[iPosition]->SetLineWidth(2);
        gr_eBy_Y[iPosition]->SetLineStyle(iPosition+1);
        gr_eBy_Y[iPosition]->Draw("CP");
    }
    TLegend *leg2 = new TLegend(0.60,0.50,0.90,0.85);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.055);
    leg2->SetLineColor(0);
    leg2->SetLineStyle(0);
    leg2->SetLineWidth(2);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1012);
    leg2->AddEntry(gr_eBy_Y[0],"y = 0fm","lp");
    leg2->AddEntry(gr_eBy_Y[1],"y = 1fm","lp");
    leg2->AddEntry(gr_eBy_Y[2],"y = 2fm","lp");
    leg2->AddEntry(gr_eBy_Y[3],"y = 3fm","lp");
    leg2->AddEntry(gr_eBy_Y[4],"y = 4fm","lp");
    leg2->AddEntry(gr_eBy_Y[5],"y = 5fm","lp");
    leg2->Draw();
   

    TLegend *leg3 = new TLegend(0.50,0.35,0.90,0.48);
    leg3->SetTextFont(42);
    leg3->SetTextSize(0.055);
    leg3->SetLineColor(0);
    leg3->SetLineStyle(0);
    leg3->SetLineWidth(2);
    leg3->SetFillColor(0);
    leg3->SetFillStyle(1012);
    leg3->SetHeader("AuAu, #sqrt{S_{NN}}=4GeV ","C");
    leg3->Draw();
    
    
   c1->SaveAs("fig4_AuAu-4GeV.pdf");
   c1->SaveAs("fig4_AuAu-4GeV.png");
}
