/*
This macro read the test.f14 file from urqmd_3.4 after used the urqmdtoroot.C macro 
to create a ROOT file. The magnetic fields are calculated in the particle point aproximation
using the Lienard-Wiechert potentials. It calculates the magnetic fields in y-direction which
is normal to reaction plane formed by the impact parameter and beam direction. 

Author: Alejandro Guirado.
Date: 25/03/20.
*/
void MagneticField_TemporalEvol_OutReactionPlaneGIF(){
 //StylePlots();

  const int tI = 0.0;
  const int tF = 10.0;
  const float tSteps = 0.5;
  const int tTotal = 20;
  const int numB = 5.;
  const int protonNumber = 158.;

  gROOT->Reset();
  TChain mychain("T");
    mychain.Add("/home/capyc/NICA/TESIS_MAESTRIA/Cuarto_Semestre/Tesis/MagneticFields/Espectadores/EvolTemporal/Datos/10000_t0-10_E4test.root");
    Int_t nlines = 0;
    struct particula_t 
    {
      Float_t time,X,Y,Z,E,Px,Py,Pz,Pt,P,m,id,isoespin,charge,lastcoll,numbercoll,history,frezetime,frezeX,frezeY,frezeZ,frezeE,frezePx,frezePy,frezePz,b,nspec,R,PXR,eBx,eBy,eBz;
    } PARTICLE;
    
    particula_t  particle;
    mychain.SetBranchAddress("particle",&particle);
    Int_t nevent = mychain.GetEntries();

    int particlesTotal  = 0.;

    double_t nPoints1 = 0,nPoints2 = 0, nPoints3 = 0,nPoints5 = 0,nPoints7 = 0,nPoints9 = 0.;
    double_t xPos1 = 0, xPos2 = 0, xPos3 = 0, xPos5 = 0, xPos7 = 0, xPos9 = 0.;
    double_t zPos1 = 0, zPos2 = 0, zPos3 = 0, zPos5 = 0, zPos7 = 0, zPos9 = 0.;
    double_t eBy1 = 0, eBy2 = 0,  eBy3 = 0, eBy5 = 0, eBy7 = 0, eBy9 = 0.;
    
    TGraph2D *graph1 = new TGraph2D();
    TGraph2D *graph2 = new TGraph2D();
    TGraph2D *graph3 = new TGraph2D();
    TGraph2D *graph5 = new TGraph2D();
    TGraph2D *graph7 = new TGraph2D();
    TGraph2D *graph9 = new TGraph2D();

     for(Int_t j=0; j<nevent; j++){
       mychain.GetEvent(j);
       
        if(particle.b> 0 && particle.b <2 && particle.time == tSteps && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            nPoints1++;
            xPos1 = particle.X;
            zPos1 = particle.Z;
            eBy1 = -particle.eBy;
            graph1->SetPoint(nPoints1,xPos1,zPos1,eBy1);
        }
        
        if(particle.b> 0 && particle.b <2 && particle.time == 2.0 && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            nPoints2++;
            xPos2 = particle.X;
            zPos2 = particle.Z;
            eBy2 = -particle.eBy;
            graph2->SetPoint(nPoints2,xPos2,zPos2,eBy2);
        }

        if(particle.b> 0 && particle.b <2 && particle.time == 4.0 && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            nPoints3++;
            xPos3 = particle.X;
            zPos3 = particle.Z;
            eBy3 = -particle.eBy;
            graph3->SetPoint(nPoints3,xPos3,zPos3,eBy3);
        }
        if(particle.b> 0 && particle.b <2 && particle.time == 6.0 && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            nPoints5++;
            xPos5 = particle.X;
            zPos5 = particle.Z;
            eBy5 = -particle.eBy;
            graph5->SetPoint(nPoints5,xPos5,zPos5,eBy5);
        }
        
        if(particle.b> 0 && particle.b <2 && particle.time == 8.0 && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            nPoints7++;
            xPos7 = particle.X;
            zPos7 = particle.Z;
            eBy7 = -particle.eBy;
            graph7->SetPoint(nPoints7,xPos7,zPos7,eBy7);
        }
        
        if(particle.b> 0 && particle.b <2 && particle.time == 10.0 && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
            nPoints9++;
            xPos9 = particle.X;
            zPos9 = particle.Z;
            eBy9 = -particle.eBy;
            graph9->SetPoint(nPoints9,xPos9,zPos9,eBy9);
        }
      //gStyle->SetPadRightMargin(0.3);
     }
     
    TCanvas* c1 = new TCanvas("c1","UrQMD test example",800,800);
    
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelOffset(0.01,"z");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.07,"xy"); //0.5
  gStyle->SetTitleSize(0.07,"z"); //0.5
  gStyle->SetTitleSize(0.07,"title"); //0.5
  gStyle->SetPadBottomMargin(0.14); //0.12
  gStyle->SetPadLeftMargin(0.12); // 0.12
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadRightMargin(0.2);

  c1->Divide(2,3);
    c1->SetTopMargin(0.1);
    c1->SetFillColor(0);
    
    c1->cd(1);
    graph1->SetTitle("t= 0.5fm; X[fm];Z[fm];-eB_{y}/m_{#pi}^{2} ");

    graph1->Draw("contz");

    c1->cd(2);
    graph2->SetTitle("t= 2fm; X[fm];Z[fm]; -eB_{y}/m_{#pi}^{2} ");
    graph2->Draw("contz");    

    c1->cd(3);
    graph3->SetTitle("t= 4fm; X[fm];Z[fm]; -eB_{y}/m_{#pi}^{2} ");
    graph3->Draw("contz");    

    c1->cd(4);
    graph5->SetTitle("t= 6fm; X[fm];Z[fm]; -eB_{y}/m_{#pi}^{2} ");
    graph5->Draw("contz");    

    c1->cd(5);
    graph7->SetTitle("t= 8fm; X[fm];Z[fm]; -eB_{y}/m_{#pi}^{2} ");
    graph7->Draw("contz");    

    c1->cd(6);
    graph9->SetTitle("t= 10fm; X[fm];Z[fm]; -eB_{y}/m_{#pi}^{2} ");
    graph9->Draw("contz");
    
   c1->SaveAs("4GeV_ByContzRP_b0-2.pdf");
   c1->SaveAs("4GeV_ByContzRP_b0-2.png");
    
}

