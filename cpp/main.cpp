//  *********************************************************************
//  * DISCLAIMER                                                        *
//  *                                                                   *
//  * Neither the authors of this software system, nor their employing  *
//  * institutes, nor the agencies providing financial support for this *
//  * work  make  any representation or  warranty, express or implied,  *
//  * regarding  this  software system or assume any liability for its  *
//  * use.                                                              *
//  *                                                                   *
//  * This  code  implementation is the  intellectual property  of the  *
//  * OpenGATE collaboration.                                           *
//  * By copying,  distributing  or modifying the Program (or any work  *
//  * based  on  the Program)  you indicate  your  acceptance of  this  *
//  * statement, and all its terms.                                     *
//  *********************************************************************
//
//######################################################################################
//# Authors   : Nicolas A Karakatsanis, Sadek A. Nehmeh, CR Schmidtlein                #
//#                                                                                    #
//# Program   : Bin_GATE_v1.0.c  29-JUL-2010                                           #
//#                                                                                    #
//# Objective : To read the coincidences TTree from the .root file, and generates the  #
//#             corresponding Michelogram and Projection files.                        #
//#                                                                                    #
//# Input     : Monte Carlo data from GATE and egsPET                                  #
//#                                                                                    #
//# Output    : 1 Michelogram  files according to various binning definitions          #
//#           : 2 Projection  files according to various binning definitions           #
//#                                                                                    #
//######################################################################################
//#                                                                                    #
//# This file is last modified on Nov 12, 2023 by: N. Karakatsanis                     #
//#                                                                                    #
//# The data are input from a root file produced by Gate simulating extended FOV of    #
//# mCT scanner. This scanner will have 5xFOV thus the root file contains information  #
//# on every gantry. In this case there are 5 gantries with gantryID (0 -> 4).         #
//# The ring numbers are defined based on the gantryID.                                #
//#                                                                                    #
//# The virtual rings between the blocks are taken into consideration here.            #
//#                                                                                    #
//# The central FOV is taken into consideration => N_RINGS = 55                        #
//# The maximum and minimum rings should be specified if the user wishes to change     #
//# the number or the order of gantries.                                               #
//#                                                                                    #
//# The odd rings are removed....                                                      #
//#                                                                                    #
//#                                                                                    #
//#                                                                                    #
//#  NEW WAY TO RUN:                                                                   #
//#                                                                                    #
//# HOW TO COMPILE:                                                                    #
//# 1) Compile using this command line in the terminal:                                #
//#       g++ (name of file) `root-config --cflags --libs`                             #
//#                                                                                    #
//# HOW TO RUN:                                                                        #
//# 1) After compiling the code type the following command line:                       #
//#       ./a.out 'directory name of root files' 'output file name' 'yes/no'           #
//#                                                                                    #
//# NOTE: 1) To drop the odd ring put 'yes' as an argument when running, otherwise the #
//# rings will not be droped.                                                          #
//#       2) You will be prompted to enter the minimum ring number                     #
//######################################################################################



#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TH2D.h"
#include "TDirectory.h"
#include "TList.h"
#include "Rtypes.h"
#include "TChainElement.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TRandom.h"
#include <time.h>
#include <cmath>

// PETSIRD Includes
#include "protocols.h"

using namespace std ;

// ETSI PET scanner model (cylindricalPET)
#define N_STAT_RINGS            40 // ETSIPETscanner: Number of physical rings
#define N_SEG                   79 // ETSIPETscanner: Number of segments (2 N_RINGS - 1)
#define N_DET                  600 // ETSIPETscanner: Detectors per ring
#define S_WIDTH                380 // ETSIPETscanner: Number of radial sinogram bins ### CHECK ###
#define N_RSEC                  60 // ETSIPETscanner: Number of resector
#define N_MODULE                16 // ETSIPETscanner: Number of modules (2x8)
#define N_MOD_xy                 2 // ETSIPETscanner: Number of tangential modules
#define N_MOD_z                  8 // ETSIPETscanner: Number of axial modules
#define N_SUBMOD                 1 // ETSIPETscanner: Number of submodules (1x1)
#define N_SMOD_xy                1 // ETSIPETscanner: Number of tangential submodules
#define N_SMOD_z                 1 // ETSIPETscanner: Number of axial submodules
#define N_CRYSTAL               25 // ETSIPETscanner: Number of crystals (5x5)
#define N_CRY_xy                 5 // ETSIPETscanner: Number of tangential crystals
#define N_CRY_z                  5 // ETSIPETscanner: Number of axial crystals
#define MAX_D_RING              39 // ETSIPETscanner: Maximum ring difference
#define SPAN                     1 // ETSIPETscanner: Span factor or axial compression factor (can only be an odd number)
#define N_PLANES              1600 // ETSIPETscanner: Total number of sinograms
#define FOV                    128 // ETSIPETscanner: Width of the FOV (mm)
#define SLICES_PER_FOV          79 // ETSIPETscanner: Number of slices per FOV
#define USE_OFFSET               0 // ETSIPETscanner: On/Off use of offset
#define OFFSET                 139 // ETSIPETscanner: Sets initial sinogram angle
#define N_ROOTFILES             36 // Number of ROOT files to process. Set it to 0 if you wish to skip ROOT binning (e.g. to bin only normalization data)
#define N_CBM_STEPS              1 // Number of CBM steps including the initial position (1 for static, 10 for VisionSparse/VisionCheckerboard every-other-mini-block, 20 $
#define CBM_DIRECTION_sign      -1 // Set it to -1 if CBM (i.e. object) motion is to the left or lower ring indices (then shift rings to the right i.e. higher ring indice$
                                   // Otherwise set it to 1, for CBM (object) motion to the right i.e. to higher ring indices, (then shift rings to the left, i.e. lower r$

#define BIN_TO_MICH             0 // Bin prompts to a Michelogram 4-Dimensional array [ring2][ring1][phi][u] NOTE: NOT PREFERRED IF MAX_D_RING<N_RINGS -1 AS THEN ITS SIZ$
#define BIN_TO_PROJ             0 // Bin prompts to Projectiongram (Viewgram) format and segment order: -MAX.RING.DIFF, ...,-2,-1,0,+1,+2, ..., +MAX.RING.DIFF
#define BIN_TO_SINO             0 // Bin prompts to Sinogram format in the following segment order: 0,+1,-1,+2,-2, ..., +MAX.RING.DIFF, -MAX.RING.DIFF

#define TxVirtualCrystalNUM      0 // ETSIPETscanner: Set the number of virtual crystals that you wish to interleave transaxially in the gaps between TxPhysCrystalNUM physical crystals
#define AxVirtualCrystalNUM      0 // ETSIPETscanner: Set the number of virtual crystals that you wish to interleave axially in the gaps between AxPhysCrystalNUM physical crystal rings
#define TxPhysCrystalNUM        10 // ETSIPETscanner: Number of physical crystals every which you wish to add TxVirtualCrystalNUM virtual crystals transaxially
#define AxPhysCrystalNUM         5 // ETSIPETscanner: Number of physical crystal rings every which you wish to add  AxVirtualCrystalNUM virtual crystal rings

#define CALC_SENS                0 // ETSIPETscanner: Set this flag to 1 only if you wish to calculate and save the axial sensitivity profile

unsigned int ans,ans1;
//unsigned short ans,ans1;
unsigned short ans_SC, ans_RC, ans_SRC, ans_S, ans_R;

/*
  Current data constructions available for output.
  These are: Michelograms:  [ring1][ring2][phi][u]
  total counts, trues, scatter, and randoms can be independently collected
*/
//Float_t    Mich_r1r2fu[N_RINGS][N_RINGS][N_DET/2][S_WIDTH]={0};

void usage()
{
  std::cout << "Usage: bin_gate [options]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -r, --root-prefix <prefix>  Prefix for root files" << std::endl;
  std::cout << "  -p, --petsird-file <file>   PETSiRD file" << std::endl;
  std::cout << "  -s, --sino-file <file>      Sinogram file" << std::endl;
  std::cout << "  -l, --legacy <yes/no>       Legacy option" << std::endl;
}

int main(int argc, char** argv)
{

  std::string root_prefix  = std::string{};
  std::string petsird_file = std::string{};
  std::string sino_file = std::string{};
  std::string legacy_option = std::string("no");

  // Parse command line args:
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-r" || arg == "--root-prefix") {
      root_prefix = argv[++i];
    } else if (arg == "-p" || arg == "--petsird-file") {
      petsird_file = argv[++i];
    } else if (arg == "-s" || arg == "--sino-file") {
      sino_file = argv[++i];
    } else if (arg == "-l" || arg == "--legacy") {
      legacy_option = argv[++i];
    } else if (arg == "-h" || arg == "--help") {
      usage();
      return 0;
    } else {
      std::cerr << "Unknown argument: " << arg << std::endl;
      return 1;
    }
  }

  if (root_prefix.empty()) {
    std::cerr << "Missing root prefix" << std::endl;
    usage();
    return 1;
  }

  if (petsird_file.empty()) {
    std::cerr << "Missing petsird file" << std::endl;
    usage();
    return 1;
  }

  // Print arguments and exit
  std::cout << "root_prefix: " << root_prefix << std::endl;
  std::cout << "petsird_file: " << petsird_file << std::endl;
  if (!sino_file.empty()) {
    std::cout << "sino_file: " << sino_file << std::endl;
  }

  time_t initialTime = time(NULL);
  time_t finalTime, totalTime;

  time_t ttime1, ttime2, time3, time4, time5, time6;

  Double_t PI = acos(-1.0);

  //  Blur variables
  Double_t FWHM_theta                 = 0.33; //  Full Width Half Max of crystal bluring
  Double_t FWHM_z                     = 0.33; //  Full Width Half Max of crystal bluring
  Double_t sigma_theta                = FWHM_theta/(2.0*sqrt(2.0*log(2.0)));
  Double_t sigma_z                    = FWHM_z/(2.0*sqrt(2.0*log(2.0)));
  Double_t C_theta                    = 0.40;
  Double_t C_z                        = 0.80;

  Int_t N_RINGS = N_STAT_RINGS + N_CBM_STEPS - 1; //Calculation of total number of rings in the system matrix (accounting for added rings due to CBM)

  // Calculation of the total number of sinogram planes accounting for the MAX_D_RING AND N_RINGS AND ASSUMING SPAN=1
  Int_t S_PLANES = 0;
  for (int s = 0 ; s < 2*MAX_D_RING + 1 ; s++)
     S_PLANES += N_RINGS - (s+1)/2;
  Int_t S_PLANES_SPAN;

  printf("Maximum ring difference: %d \n",MAX_D_RING);
  printf("Span factor: %d \n",SPAN);

  string filedir, inputfilename ;
  string prompts_filename, norm_outputfilename, Moutputfilename, Poutputfilename, Soutputfilename;

  FILE  *Mich_r1r2fuFile, *Proj_File, *Sino_File;
  prompts_filename = sino_file;

  if (BIN_TO_MICH)
    {
      Moutputfilename = "Mich_"+ prompts_filename + ".s" ;
      cout << "Michelogram file name is = " << Moutputfilename << endl ;
      Mich_r1r2fuFile = fopen(Moutputfilename.c_str(),"wb");
    }

  if (BIN_TO_PROJ)
    {
      Poutputfilename = "Proj_"+ prompts_filename + ".s" ;
      cout << "Projection file name is = " << Poutputfilename << endl ;
      Proj_File = fopen(Poutputfilename.c_str(),"wb");
    }

  if (!sino_file.empty())
    {
      Soutputfilename = "Sino_"+ prompts_filename + ".s" ;
      cout << "Sinogram file name is = " << Soutputfilename << endl ;
      Sino_File = fopen(Soutputfilename.c_str(),"wb");
    }

  // This is an input user parameter to allow exclusion of odd ring counts from binning.
  string remODD;
  remODD = legacy_option;

  int RMIN, RMIN_CBM;
  int RMAX, RMAX_CBM;
  cout<<endl<<endl;

  //cout<<"Please Enter The Minimum Ring Number"<<endl<<endl;
  //cout<<">>>>>>>>>>>>>>>>:"<<endl;
  //cin>>RMIN;

  RMIN = 0;
  RMAX = N_STAT_RINGS;

  ofstream sensitivityOut("sensitivity_center.txt");
  float *Sensitivity;

  if (CALC_SENS)
     float Sensitivity[2*N_RINGS-1]={0};

  float ****Mich_r1r2fu = new float***[N_RINGS];
  for( int i = 0; i < N_RINGS; i++ )
    Mich_r1r2fu[i] = new float**[N_RINGS];

  for( int i = 0; i < N_RINGS; i++ )
    for( int j = 0; j < N_RINGS; j++ )
      Mich_r1r2fu[i][j] = new float*[N_DET/2];

  for( int i = 0; i < N_RINGS; i++ )
      for( int j = 0; j < N_RINGS; j++ )
          for( int k = 0; k < (N_DET/2); k++ )
              Mich_r1r2fu[i][j][k] = new float[S_WIDTH];

  for( int i = 0; i < N_RINGS; i++ )
      for( int j = 0; j < N_RINGS; j++ )
          for( int k = 0; k < (N_DET/2); k++ )
              for(int l = 0; l< S_WIDTH; l++)
                  Mich_r1r2fu[i][j][k][l] = 0.0;

  cout<<"The 4 dimensional array is allocated"<<endl;

  //#####################################################################
  //#             SINOGRAMS AND PROJECTION PLANES BINNING               #
  //#####################################################################
  Int_t                         ring1, ring2, crystal1, crystal2;
  Int_t                         phi, u;
  unsigned long long int        totalCBM_Counts_in_FOV = 0, totalCBM_Counts_binned = 0, Counts_in_FOV, Counts_binned;
  int                           flip, swap, zi, c1, c2;
  int                           N_ROOTFILES_PER_CBM_STEP, RING_SHIFT, ROOTFILE_INDEX;

  //#####################################################################
  //#              Loop over the .root file in the directory "PATH"     #
  //#####################################################################
  Int_t   Trues = 0, Scatters = 0, Randoms = 0;
  Int_t   nbytes = 0;


  //####################################################################
  //#             Declaration of leaves types - TTree Coincidences     #
  //####################################################################
  Float_t         		axialPos, rotationAngle, sinogramS, sinogramTheta;
  Char_t          		comptVolName1[40], comptVolName2[40];
  Int_t           		compton1, compton2, gantryID1, gantryID2;
  Int_t           		runID, sourceID1, sourceID2, eventID1, eventID2;
  Int_t           		layerID1, layerID2, crystalID1, crystalID2;
  Int_t           		submoduleID1, submoduleID2, moduleID1, moduleID2, rsectorID1, rsectorID2;
  Int_t           		comptonPhantom1, comptonPhantom2;
  Float_t         		energy1, energy2;
  Float_t         		globalPosX1, globalPosX2, globalPosY1, globalPosY2, globalPosZ1, globalPosZ2;
  Float_t         		sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2;
  Double_t        		time1, time2;
  unsigned long long int        nentries, totalCBM_nentries;
  Int_t           		counterA, counterB;
  counterA = 0;
  counterB = 0;

  //######################################################################################
  //#                        Set branch addresses - TTree Coincidences                   #
  //######################################################################################

  filedir = root_prefix;
  TChain *Coincidences = new TChain("Coincidences");

  ostringstream fileNumber[N_ROOTFILES];
  string fileN[N_ROOTFILES];
  N_ROOTFILES_PER_CBM_STEP = N_ROOTFILES/N_CBM_STEPS;

  for (int CBM_step_i = 0; CBM_step_i<N_CBM_STEPS; CBM_step_i++ )
    {
     cout << "Processing CBM step " << CBM_step_i << endl;
     //Initializations for each CBM step
     RING_SHIFT=CBM_step_i;
      //IN CASE OF RIGHT CBM DIRECTION (BED MOVES TOWARDS HIGHER RING INDICES), FIRST SHIFT ORIGINAL RING INDICES TO HIGHER VALUES TO AVOID NEGATIVE RING INDICES DURING CBM R$
      if (CBM_DIRECTION_sign>0)
        {
         RMIN_CBM=RMIN-CBM_DIRECTION_sign*RING_SHIFT+N_CBM_STEPS-1;
         RMAX_CBM=RMAX-CBM_DIRECTION_sign*RING_SHIFT+N_CBM_STEPS-1;
        }
      else //LEFT CBM DIRECTION (BED MOVES TOWARDS LOWER RING INDICES)
        {
         RMIN_CBM=RMIN-CBM_DIRECTION_sign*RING_SHIFT;
         RMAX_CBM=RMAX-CBM_DIRECTION_sign*RING_SHIFT;
        }
     Counts_in_FOV = 0;
     Counts_binned = 0;
     TChain *Coincidences = new TChain("Coincidences");


     for(int i = 0; i<N_ROOTFILES_PER_CBM_STEP; i++ )
        {
         ROOTFILE_INDEX=N_ROOTFILES_PER_CBM_STEP*CBM_step_i + i;
         fileNumber[ROOTFILE_INDEX] << ROOTFILE_INDEX + 1;
         fileN[ROOTFILE_INDEX] = fileNumber[ROOTFILE_INDEX].str();

         inputfilename = filedir + fileN[ROOTFILE_INDEX] +".root" ;
         //inputfilename = filedir +".root" ;

         cout << "Input file name is " << inputfilename << endl;

         //  TChain *Coincidences = new TChain("Coincidences") ;
         Coincidences->Add(inputfilename.c_str()) ;
        }

     /*
       inputfilename = filedir +".root" ;

       cout << "Input file name is " << inputfilename << endl;

       //  TChain *Coincidences = new TChain("Coincidences") ;
       Coincidences->Add(inputfilename.c_str()) ;

     */

  Coincidences->SetBranchStatus("*",0);
  Coincidences->SetBranchAddress("axialPos",&axialPos);
  Coincidences->SetBranchAddress("comptVolName1",&comptVolName1);
  Coincidences->SetBranchAddress("comptVolName2",&comptVolName2);
  Coincidences->SetBranchAddress("comptonCrystal1",&compton1);
  Coincidences->SetBranchAddress("comptonCrystal2",&compton2);
  Coincidences->SetBranchAddress("crystalID1",&crystalID1);
  Coincidences->SetBranchAddress("crystalID2",&crystalID2);
  Coincidences->SetBranchAddress("comptonPhantom1",&comptonPhantom1);
  Coincidences->SetBranchAddress("comptonPhantom2",&comptonPhantom2);
  Coincidences->SetBranchAddress("energy1",&energy1);
  Coincidences->SetBranchAddress("energy2",&energy2);
  Coincidences->SetBranchAddress("eventID1",&eventID1);
  Coincidences->SetBranchAddress("eventID2",&eventID2);
  Coincidences->SetBranchAddress("globalPosX1",&globalPosX1);
  Coincidences->SetBranchAddress("globalPosX2",&globalPosX2);
  Coincidences->SetBranchAddress("globalPosY1",&globalPosY1);
  Coincidences->SetBranchAddress("globalPosY2",&globalPosY2);
  Coincidences->SetBranchAddress("globalPosZ1",&globalPosZ1);
  Coincidences->SetBranchAddress("globalPosZ2",&globalPosZ2);
  Coincidences->SetBranchAddress("layerID1",&layerID1);
  Coincidences->SetBranchAddress("layerID2",&layerID2);
  Coincidences->SetBranchAddress("moduleID1",&moduleID1);
  Coincidences->SetBranchAddress("moduleID2",&moduleID2);
  Coincidences->SetBranchAddress("rotationAngle",&rotationAngle);
  Coincidences->SetBranchAddress("rsectorID1",&rsectorID1);
  Coincidences->SetBranchAddress("rsectorID2",&rsectorID2);
  Coincidences->SetBranchAddress("runID",&runID);
  Coincidences->SetBranchAddress("sinogramS",&sinogramS);
  Coincidences->SetBranchAddress("sinogramTheta",&sinogramTheta);
  Coincidences->SetBranchAddress("sourceID1",&sourceID1);
  Coincidences->SetBranchAddress("sourceID2",&sourceID2);
  Coincidences->SetBranchAddress("sourcePosX1",&sourcePosX1);
  Coincidences->SetBranchAddress("sourcePosX2",&sourcePosX2);
  Coincidences->SetBranchAddress("sourcePosY1",&sourcePosY1);
  Coincidences->SetBranchAddress("sourcePosY2",&sourcePosY2);
  Coincidences->SetBranchAddress("sourcePosZ1",&sourcePosZ1);
  Coincidences->SetBranchAddress("sourcePosZ2",&sourcePosZ2);
  Coincidences->SetBranchAddress("submoduleID1",&submoduleID1);
  Coincidences->SetBranchAddress("submoduleID2",&submoduleID2);
  Coincidences->SetBranchAddress("time1",&time1);
  Coincidences->SetBranchAddress("time2",&time2);
  Coincidences->SetBranchAddress("gantryID1",&gantryID1);
  Coincidences->SetBranchAddress("gantryID2",&gantryID2);

  nentries = (unsigned long long int)(Coincidences->GetEntries());

  printf("Total Number of Coincidence Events in the ROOT file:= %llu \n",nentries );

  for (unsigned long long int i = 0 ; i < nentries ; i++)
    {
      if ((i%250000)  == 0 && i!=0)  printf("... %llu ",i);
      if ((i%1000000) == 0 && i!=0)  printf("... %llu\n",i);

      nbytes += Coincidences->GetEntry(i);

      // Update the number of Trues and Randoms...
      //------------------------------------------
      if (abs(sinogramS) > FOV/2)
	continue;

      if (eventID1 == eventID2)
	{
	  if (comptonPhantom1 == 0 && comptonPhantom2 == 0) Trues++;
	  else Scatters++;
	}
      else Randoms++;

      //-----------------------------------
      //  Identify the ring#...
      //-----------------------
      ring1 = (Int_t)(crystalID1/N_CRY_xy)
	    + (Int_t)(submoduleID1/N_SMOD_xy)*N_CRY_z
	    + (Int_t)(moduleID1/N_MOD_xy)*N_SMOD_z*N_CRY_z
            + (Int_t)(gantryID1)*52;
      ring2 = (Int_t)(crystalID2/N_CRY_xy)
	    + (Int_t)(submoduleID2/N_SMOD_xy)*N_CRY_z
	    + (Int_t)(moduleID2/N_MOD_xy)*N_SMOD_z*N_CRY_z
            + (Int_t)(gantryID2)*52;

      //----------------------------------
      //   Take into account the virtual rings, change the order of the actual rings
      //---------------------------------------------------

      //-------------------------
      //  Inserting virtual crystals axially: There are AxVirtualCrystalNUM virtual crystals every AxPhysCrystalNUM physical crystals axially
      //----------------------------------
      ring1 = ring1 + AxVirtualCrystalNUM*floor(ring1/AxPhysCrystalNUM);
      ring2 = ring2 + AxVirtualCrystalNUM*floor(ring2/AxPhysCrystalNUM);


      if ( abs(ring1 - ring2) > MAX_D_RING )  continue;

      //-----------------------
      //  Identify the crystal#...
      //-----------------------------------
      crystal1 = rsectorID1 * N_MOD_xy * N_SMOD_xy * N_CRY_xy
	       + (moduleID1%N_MOD_xy) * N_SMOD_xy * N_CRY_xy
	       + (submoduleID1%N_SMOD_xy) * N_CRY_xy
	       + (crystalID1%N_CRY_xy);
      crystal2 = rsectorID2 * N_MOD_xy * N_SMOD_xy * N_CRY_xy
	       + (moduleID2%N_MOD_xy) * N_SMOD_xy * N_CRY_xy
	       + (submoduleID2%N_SMOD_xy) * N_CRY_xy
	       + (crystalID2%N_CRY_xy);

      //-------------------------
      //  Inserting virtual crystals transaxially: There are TxVirtualCrystalNUM virtual crystals every TxPhysCrystalNUM physical crystals transaxially
      //----------------------------------
      crystal1 = crystal1 + TxVirtualCrystalNUM*floor(crystal1/TxPhysCrystalNUM);
      crystal2 = crystal2 + TxVirtualCrystalNUM*floor(crystal2/TxPhysCrystalNUM);

      //-----------------------------------------------------
      //  Rotate the image correctly#...
      //--------------------------------
      if (USE_OFFSET == 1)
	{
	  crystal1 = crystal1 + OFFSET;
	  crystal2 = crystal2 + OFFSET;
	  if (crystal1 >= N_DET)  crystal1 = crystal1 - N_DET;
	  if (crystal2 >= N_DET)  crystal2 = crystal2 - N_DET;
	}
      //--------------------------------
      //  Bin the crystal ring pairs into Michelograms
      //  u - radial sinogram component
      //  phi - azimuthal sinogram component
      //  ring pairs are sorted according to c1 < c2 else flip
      //  where c1 and c2 are crystals at phi(u = S_WIDTH/2)
      //--------------------------------
      phi = ((crystal1 + crystal2 + N_DET/2)%N_DET)/2;

      if (((crystal1 + crystal2) < (3*N_DET/2)) && ((crystal1 + crystal2) >= (N_DET/2)))
	u    =  abs(crystal1 - crystal2) -  N_DET/2 + S_WIDTH/2;
      else u = -abs(crystal1 - crystal2) +  N_DET/2 + S_WIDTH/2;

      if ( u >= S_WIDTH || u < 0 ) continue;

      if (u%2 == 0)
	{
	  zi = (N_DET/2 - (crystal1 - crystal2) - 1)/2;
	  if (zi >=  N_DET/4) zi = zi - N_DET/2 + 1;
	  if (zi <= -N_DET/4) zi = zi + N_DET/2 - 1;
	}
      else
	{
	  zi = (N_DET/2 - (crystal1 - crystal2))/2;
	  if (zi >=  N_DET/4) zi = zi - N_DET/2;
	  if (zi <= -N_DET/4) zi = zi + N_DET/2;
	}

      c1 = crystal1 + zi;
      c2 = crystal2 - zi;
      if (c1 >= N_DET) c1 = c1 - N_DET;
      if (c1 < 0)      c1 = c1 + N_DET;
      if (c2 >= N_DET) c2 = c2 - N_DET;
      if (c2 < 0)      c2 = c2 + N_DET;
      if (c1 < c2) flip = 0;
      else         flip = 1;

      if (flip)
	{
	  swap  = ring1;
	  ring1 = ring2;
	  ring2 = swap;
	}

     // Update the different arrays...
     //-------------------------------
      //ALL EVENTS
      if ( ( (eventID1 == eventID2 ) ))
	{
	 //true+scatter
        if (comptonPhantom1 == 0 && comptonPhantom2 == 0)
         {
          if ( (ring1<RMAX) && (ring2<RMAX) && (ring1>RMIN-1) && (ring2>RMIN-1) )
            {
               if(remODD == "yes")
                {
                if ((ring1%2==0) && (ring2%2==0))
                // if ((ring1%2!=0) && (ring2%2!=0))
                      {
                       Mich_r1r2fu[ring2-RMIN][ring1-RMIN][phi][u] += 1.;
                       Counts_binned++;
                      }
                     else
                       Mich_r1r2fu[ring2-RMIN][ring1-RMIN][phi][u] += 0.;
                 }
               else
                {
                 Mich_r1r2fu[ring2-RMIN][ring1-RMIN][phi][u] += 1.;
                 Counts_binned++;
                }
	     //true
	   }
        }
       else
	 {
	  //scatter
	 }
       }
     else
       {
	 //random
       }
     Counts_in_FOV++;
    }

     printf("\n");
     if (N_CBM_STEPS>1)
       {
        printf("CBM STEP %d has been processed...\n",CBM_step_i);
        totalCBM_nentries=totalCBM_nentries+nentries;
        totalCBM_Counts_in_FOV=totalCBM_Counts_in_FOV+Counts_in_FOV;
        totalCBM_Counts_binned=totalCBM_Counts_binned+Counts_binned;
       }
     printf("Total Number of Coincidence Events in the ROOT file:= %llu ...\n",nentries );
     printf("Total Number of Coincidence Events in the FOV:= %llu ...\n",Counts_in_FOV );
     printf("Total Number of Coincidence Events registered in list-mode or sinogram format:= %llu ...\n",Counts_binned );
    } //END FOR-LOOP ACROSS N_CBM_STEPS

 if (N_CBM_STEPS>1)
   {
    printf("Total Number of Coincidence Events in the ROOT file (across all %d CBM steps) := %llu ...\n",N_CBM_STEPS,totalCBM_nentries );
    printf("Total Number of Coincidence Events in the FOV (across all %d CBM steps) := %llu ...\n",N_CBM_STEPS,totalCBM_Counts_in_FOV );
    printf("Total Number of Coincidence Events registered in list-mode or sinogram format (across all %d CBM steps) := %llu ...\n",N_CBM_STEPS,totalCBM_Counts_binned );
   }


  // Write the data to disk, and then close Michelogram file...
  //----------------------------------------------------------------
  //ans = fwrite(Mich_r1r2fu,4,(N_RINGS*N_RINGS*N_DET/2*S_WIDTH),Mich_r1r2fuFile);


  ttime1 = time(NULL);
  cout<<"time to finish binnig = "<<(-initialTime + ttime1)/60 <<"min"<<endl;
/// Applying Symmetries


  if (BIN_TO_MICH)
    {
      printf("Writing Michelogram File to Disk...\n");
      float DummyArray[1] = {0.0};
      for( int i = 0; i < N_RINGS; i++ )
        for( int j = 0; j < N_RINGS; j++ )
          for( int k = 0; k < (N_DET/2); k++ )
             for(int l = 0; l< S_WIDTH; l++)
               {
                DummyArray[0] = Mich_r1r2fu[i][j][k][l];
                ans = fwrite(DummyArray,4,1.0,Mich_r1r2fuFile);
               }
      fclose(Mich_r1r2fuFile);
    }

  // Calculating sensitivity profile
  if (CALC_SENS)
    {
     printf("   Calculating the sensitivity axial profile...\n");
     for(int slice = 0; slice < 2*N_RINGS-1; slice++)
       {
        for(int r2 = 0, r1 = slice; r2 < slice, r1>=0; r2++, r1--)
            for(int Phi = 0; Phi < N_DET/2; Phi++ )
                for(int U = 0; U < S_WIDTH; U++)
                   if(r1<N_RINGS && r2<N_RINGS)
                      Sensitivity[slice]+= Mich_r1r2fu[r1][r2][Phi][U];

        sensitivityOut<< slice <<"       " << Sensitivity[slice]<<endl;
       }
    }

  // Generate projection files...
  //-----------------------------
  // From Segment number -MAX_D_RING to +MAX_D_RING
  // After phi
  // After z
  // After r
  if (BIN_TO_PROJ)
    {
      printf("Writing Projection File to Disk...\n");
      Int_t S_NUM;
      for (int i = 0 ; i < 2*MAX_D_RING + 1 ; i++)
        {
          if (i <= MAX_D_RING) S_NUM = N_RINGS - MAX_D_RING + i;
          else                 S_NUM = N_RINGS + MAX_D_RING - i;

          for (int j = 0 ; j < (N_DET/2) ; j++)
            {
              float Proj[S_NUM][S_WIDTH];
              for (int k = 0 ; k < S_NUM ; k++)
                {
                  if (i <= MAX_D_RING) ring1 = k;
                  else                 ring2 = k;

                  if (i <= MAX_D_RING) ring2 = ring1 + MAX_D_RING - i;
                  else                 ring1 = ring2 - MAX_D_RING + i;

                  for (int l = 0 ; l < S_WIDTH ; l++)
                    Proj[k][l] = Mich_r1r2fu[ring2][ring1][j][l];
                }
              ans = fwrite(Proj,4,(S_NUM*S_WIDTH),Proj_File);
            }
          // printf("\n");
        }
      fclose(Proj_File);
    } //END IF

  if (!sino_file.empty())
    {
      printf("Writing Sinogram to Disk...\n");
      Int_t S_NUM;
      for (int i = 0 ; i < 2*MAX_D_RING + 1 ; i++)
        {
          S_NUM = N_RINGS - (i+1)/2;

          for (int k = 0 ; k < S_NUM ; k++)
            {
              if (!(i%2))       ring1 = k; //zero or negative ring diference ring1-ring2
              else      ring2 = k; //positive ring difference ring1-ring2

              if (!(i%2))       ring2 = ring1 +  (i/2); //zero or negative ring diference ring1-ring2
              else              ring1 = ring2 + (i/2) + 1; //positive ring difference ring1-ring2

              float Sino[N_DET/2][S_WIDTH];

              for (int j = 0 ; j < (N_DET/2) ; j++)
                for (int l = 0 ; l < S_WIDTH ; l++)
                  Sino[j][l] = Mich_r1r2fu[ring2][ring1][j][l];

              ans = fwrite(Sino,4,((N_DET/2)*S_WIDTH),Sino_File);
            }
          // printf("\n");
        }
      fclose(Sino_File);
    } //END IF

finalTime = time(NULL);
totalTime = -initialTime + finalTime;
cout<<"total calculation time in seconds = "<< totalTime <<" sec"<<endl;
cout<<"total calculation time in minutes = "<< totalTime/60 <<" min"<<endl;

  return(0);
}
