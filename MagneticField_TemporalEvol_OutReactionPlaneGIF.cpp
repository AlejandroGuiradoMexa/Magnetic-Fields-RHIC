

/*
This macro read the test.f14 file from urqmd_3.4 after used the urqmdtoroot.C macro 
to create a ROOT file. The magnetic fields are calculated in the particle point aproximation
using the Lienard-Wiechert potentials. After run this macro is needed to use a ROOT function to 
create the magnetic fields evolution in time GIF, this can be done in "PlotGIF.cpp".

Author: Alejandro Guirado.
Date: 25/03/20.
*/
void StylePlots(){
    gROOT->Reset();
  gStyle->SetCanvasColor(-1);
  gStyle->SetPadColor(-1);
  gStyle->SetFrameFillColor(-1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetHistFillColor(-1);
  gStyle->SetTitleFillColor(-1);
  gStyle->SetFillColor(-1);
  gStyle->SetFillStyle(4000);
  gStyle->SetStatStyle(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetPadColor(10); //
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLabelSize(0.035,"xyz");
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTitleSize(0.04,"xy"); //0.5
  gStyle->SetTitleSize(0.04,"z"); //0.5
  gStyle->SetLabelColor(kBlack,"xyz");
}
void MagneticField_TemporalEvol_OutReactionPlaneGIF(){

  const int nTime = 20, binNumber = 50;
  const int tI = 0.0;
  const int tF = 10.0;
  const float tSteps = 0.5;

    TCanvas *canvas[nTime] = {NULL};
    TGraph2D *ByColz[nTime] = {NULL};
    char *mass_p = new char[10];
    char *legendCanvas = new char[10];

    for(int index=0; index<nTime; index++){
       ByColz[index] = new TGraph2D();
    }

    double_t nPoints[nTime] = {0.};
    double_t xPos[nTime] = {0.};
    double_t yPos[nTime] = {0.};
    double_t eBy[nTime] = {0.};
  TChain mychain("T");

    mychain.Add("10-pruebaGIFtest.root");
    Int_t nlines = 0;
    struct particula_t 
    {
      Float_t time,X,Y,Z,E,Px,Py,Pz,Pt,P,m,id,isoespin,charge,lastcoll,numbercoll,history,frezetime,frezeX,frezeY,frezeZ,frezeE,frezePx,frezePy,frezePz,b,nspec,R,PXR,eBx,eBy,eBz;
    } PARTICLE;
    
    particula_t  particle;
    mychain.SetBranchAddress("particle",&particle);
    Int_t nevent = mychain.GetEntries();
     for(Int_t j=0; j<nevent; j++){
       mychain.GetEvent(j);
       
        if(particle.b> 6 && particle.b <8 && particle.numbercoll==0 && particle.charge==1 && particle.id == 1){
      for(int time = 1; time<=nTime; time++){
          double timeUrqmd = time*0.5;
            if( timeUrqmd == particle.time){
                nPoints[time-1]++;
                xPos[time-1] = particle.X;
                yPos[time-1] = particle.Y;
                eBy[time-1] = -particle.eBy;   
                cout << particle.time <<"    " <<nPoints[time-1]
                << "   " <<  xPos[time-1] << "  "<<  eBy[time-1] <<endl;
                ByColz[time-1]->SetPoint(nPoints[time-1],xPos[time-1],yPos[time-1],eBy[time-1]);
            }
          }
        }
     }
    TLegend *legend[nTime] = {NULL};
    for(int index=0; index<nTime; index++){
       legend[index] = new TLegend(0.1,0.7,0.4,0.9);
       legend[index]->SetFillStyle(0);

    }

    TCanvas *CanvasGif = new TCanvas("CanvasGif","", 1200,800);
    double_t canvasIndex = 0.;
    for(int index=0; index<nTime; index++){
        ostringstream MassIndex,legendIndex;
        MassIndex << "canvas" << index;
        legendIndex<< "legend"<<index;
        sprintf(mass_p,"Plots/mass_p_%d.png",index);  
        canvasIndex = (index+1)*0.5; 
        cout << canvasIndex << "   "<< index <<  "  "<< (index+1)/2<<endl;                           
        sprintf(legendCanvas," t = %8.1f  [fm]",canvasIndex);                               

        legend[index]->AddEntry((TObject*)0,legendCanvas,"");
        ByColz[index]->SetMaximum(0.1);
        ByColz[index]->SetMinimum(-0.7);
        canvas[index] = new TCanvas(MassIndex.str().c_str(),"",1200,800 );
        ByColz[index]->SetTitle("Temporal evolution of MF out of reaction plane");
        ByColz[index]->GetXaxis()->SetTitle("X[fm]");
        ByColz[index]->GetYaxis()->SetTitle("Y[fm]");
        ByColz[index]->GetZaxis()->SetTitle("eBy[GeV]");
        ByColz[index]->Draw("colz");
        legend[index]->Draw();
        canvas[index]->SaveAs(mass_p);
        cout<< mass_p<<endl;
    }

    for(int index=0; index<nTime; index++){
        delete canvas[index], ByColz[index];
    }
}