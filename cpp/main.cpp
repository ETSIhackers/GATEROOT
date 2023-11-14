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
#include "types.h"
#include "binary/protocols.h"

using namespace std ;

// ETSI PET scanner model (cylindricalPET)
#define N_STAT_RINGS            40 // ETSIPETscanner: Number of physical rings
#define N_SEG                   79 // ETSIPETscanner: Number of segments (2 N_RINGS - 1)
#define N_DET                  600 // ETSIPETscanner: Detectors per ring
#define S_WIDTH                380 // ETSIPETscanner: Number of radial sinogram bins ### CHECK ###
#define N_RSEC                  60 // ETSIPETscanner: Number of resector
#define N_RSEC_xy               60 // BIOGRAPH VISION: Number of rsectors per ring (transaxially)
#define N_RSEC_z                 1 // BIOGRAPH VISION: Regular
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


void usage()
{
  std::cout << "Usage: bin_gate [options]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -r, --root-prefix <prefix>  Prefix for root files" << std::endl;
  std::cout << "  -p, --petsird-file <file>   PETSiRD file" << std::endl;
  std::cout << "  -s, --sino-file <file>      Sinogram file" << std::endl;
  std::cout << "  -l, --legacy <yes/no>       Legacy option" << std::endl;
}

int calculate_detector_id(int gantry_id, int rsector_id, int module_id, int submodule_id, int crystal_id, int rmin = 0)
{
  int ring = (Int_t)(gantry_id)*N_RSEC_z*N_MOD_z*N_SMOD_z*N_CRY_z
        + (Int_t)(rsector_id/N_RSEC_xy)*N_MOD_z*N_SMOD_z*N_CRY_z
        + (Int_t)(module_id/N_MOD_xy)*N_SMOD_z*N_CRY_z
        + (Int_t)(submodule_id/N_SMOD_xy)*N_CRY_z
        + (Int_t)(crystal_id/N_CRY_xy);

  int crystal = (Int_t)(crystal_id%N_CRY_xy)
        + (Int_t)(submodule_id%N_SMOD_xy)*N_CRY_xy
        + (Int_t)(module_id%N_MOD_xy)*N_SMOD_xy*N_CRY_xy
        + (Int_t)(rsector_id%N_RSEC_xy)*N_MOD_xy*N_SMOD_xy*N_CRY_xy;

  return crystal + (ring-rmin)*N_DET;
}

// single ring as example
prd::ScannerInformation
get_scanner_info(float radius, int n_detectors, int n_rings)
{
  const int NUMBER_OF_TOF_BINS = 25;
  const int NUMBER_OF_ENERGY_BINS = 1;

  std::vector<float> angles;
  for (int i = 0; i < n_detectors; ++i)
  {
    angles.push_back(static_cast<float>(2 * M_PI * (1.0f*i) / n_detectors));
  }

  std::vector<prd::Detector> detectors;
  int detector_id = 0;
  for (int r =0; r < n_rings; r++)
  {
    for (auto angle : angles)
    {
      // Create a new detector
      prd::Detector d;
      d.x = radius * std::sin(angle);
      d.y = radius * std::cos(angle);
      d.z = 0.0+r; // TODO: calculate this
      d.id = detector_id++; // TODO: Make consistent with detector ID calculation above
      detectors.push_back(d);
    }
  }

  typedef yardl::NDArray<float, 1> FArray1D;
  // TOF info (in mm)
  FArray1D::shape_type tof_bin_edges_shape = { NUMBER_OF_TOF_BINS + 1 };
  FArray1D tof_bin_edges(tof_bin_edges_shape);
  for (std::size_t i = 0; i < tof_bin_edges.size(); ++i)
    tof_bin_edges[i] = (i - NUMBER_OF_TOF_BINS / 2.F) / NUMBER_OF_TOF_BINS * 2 * radius;
  FArray1D::shape_type energy_bin_edges_shape = { NUMBER_OF_ENERGY_BINS + 1 };
  FArray1D energy_bin_edges(energy_bin_edges_shape);
  for (std::size_t i = 0; i < energy_bin_edges.size(); ++i)
    energy_bin_edges[i] = 430.F + i * (650.F - 430.F) / NUMBER_OF_ENERGY_BINS;
  prd::ScannerInformation scanner_info;
  scanner_info.detectors = detectors;
  scanner_info.tof_bin_edges = tof_bin_edges;
  scanner_info.tof_resolution = 9.4F; // in mm
  scanner_info.energy_bin_edges = energy_bin_edges;
  scanner_info.energy_resolution_at_511 = .11F;    // as fraction of 511
  scanner_info.listmode_time_block_duration = 1.F; // ms
  return scanner_info;
}

uint32_t tofToIdx(float tof, const prd::ScannerInformation& scanner_info)
{
  float tof_mm = tof * 0.03;
  for (size_t i = 0; i < scanner_info.tof_bin_edges.size() - 1; ++i)
  {
    if (tof_mm >= scanner_info.tof_bin_edges[i] && tof_mm < scanner_info.tof_bin_edges[i+1])
    {
      return static_cast<uint32_t>(i);
    }
  }
  throw std::runtime_error("TOF out of range");
}

uint32_t energyToIdx(float energy, const prd::ScannerInformation& scanner_info)
{
  for (size_t i = 0; i < scanner_info.energy_bin_edges.size() - 1; ++i)
  {
    if (energy >= scanner_info.energy_bin_edges[i] && energy < scanner_info.energy_bin_edges[i+1])
    {
      return static_cast<uint32_t>(i);
    }
  }
  std::stringstream ss;
  ss << "Energy out of range: " << energy;
  throw std::runtime_error(ss.str());
}

int main(int argc, char** argv)
{

  std::string root_prefix  = std::string{};
  std::string petsird_file = std::string{};
  std::string sino_file = std::string{};
  std::string legacy_option = std::string("no");
  int number_of_root_files = 2;

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
    } else if (arg == "-n" || arg == "--number-of-root-files") {
      number_of_root_files = atoi(argv[++i]);
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

  Int_t N_RINGS = N_STAT_RINGS + N_CBM_STEPS - 1; //Calculation of total number of rings in the system matrix (accounting for added rings due to CBM)

  string filedir, inputfilename ;
  string prompts_filename, norm_outputfilename;


  // This is an input user parameter to allow exclusion of odd ring counts from binning.
  string remODD;
  remODD = legacy_option;

  //#####################################################################
  //#              Loop over the .root file in the directory "PATH"     #
  //#####################################################################
  Int_t   Trues = 0, Scatters = 0, Randoms = 0;

  //####################################################################
  //#             Declaration of leaves types - TTree Coincidences     #
  //####################################################################
  Float_t         		axialPos, rotationAngle, sinogramS, sinogramTheta;
  Char_t          		comptVolName1[255], comptVolName2[255];
  Int_t           		compton1, compton2, gantryID1, gantryID2;
  Int_t           		runID, sourceID1, sourceID2, eventID1, eventID2;
  Int_t           		layerID1, layerID2, crystalID1, crystalID2;
  Int_t           		submoduleID1, submoduleID2, moduleID1, moduleID2, rsectorID1, rsectorID2;
  Int_t           		comptonPhantom1, comptonPhantom2;
  Float_t         		energy1, energy2;
  Float_t         		globalPosX1, globalPosX2, globalPosY1, globalPosY2, globalPosZ1, globalPosZ2;
  Float_t         		sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2;
  Double_t        		time1, time2;
  unsigned long long int        nentries;

  //######################################################################################
  //#                        Set branch addresses - TTree Coincidences                   #
  //######################################################################################

  filedir = root_prefix;
  TChain *Coincidences = new TChain("Coincidences");

  for (int i = 0; i < number_of_root_files; i++) {
    std::ostringstream fileNumber;
    fileNumber << i + 1;
    inputfilename = filedir + fileNumber.str() + ".root";
    std::cout << "Input file name is " << inputfilename << std::endl;
    Coincidences->Add(inputfilename.c_str());
  }


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

  // Output PETSIRD
  prd::Header header;
  prd::ScannerInformation scanner = get_scanner_info(300.0f, N_DET, N_RINGS);
  prd::ExamInformation exam;

  header.exam = exam;
  header.scanner = scanner;

  // Write PETSiRD file
  prd::binary::PrdExperimentWriter writer(petsird_file);
  writer.WriteHeader(header);

  long current_time_block = -1;
  prd::TimeBlock time_block;
  unsigned long Counts_binned = 0;
  for (unsigned long long int i = 0 ; i < nentries ; i++)
  {
    if (i % 1000000 == 0) {
      printf("Processing event %llu of %llu\n", i, nentries);
    }

    Coincidences->GetEntry(i);
    if (eventID1 == eventID2)
    {
	    if (comptonPhantom1 == 0 && comptonPhantom2 == 0) {
        prd::CoincidenceEvent event;
        event.detector_1_id = calculate_detector_id(gantryID1, rsectorID1, moduleID1, submoduleID1, crystalID1);
        event.detector_2_id = calculate_detector_id(gantryID2, rsectorID2, moduleID2, submoduleID2, crystalID2);
        event.tof_idx = static_cast<uint32_t>(tofToIdx(1.0e12f*(time1 - time2), scanner));
        event.energy_1_idx = static_cast<uint32_t>(energyToIdx(1.0e3*energy1, scanner));
        event.energy_2_idx = static_cast<uint32_t>(energyToIdx(1.0e3*energy2, scanner));

        if (false && i%100000 == 0) {
          std::cout << "Event " << i << std::endl;
          std::cout << "  detector_1_id: " << event.detector_1_id << std::endl;
          std::cout << "  detector_2_id: " << event.detector_2_id << std::endl;
          std::cout << "  tof_idx: " << event.tof_idx << std::endl;
          std::cout << "  energy_1_idx: " << event.energy_1_idx << std::endl;
          std::cout << "  energy_2_idx: " << event.energy_2_idx << std::endl;
          std::cout << "  energy_1: " << energy1 << std::endl;
          std::cout << "  energy_2: " << energy2 << std::endl;
        }

        long this_time_block = static_cast<long>(time1*1.0e3 / scanner.listmode_time_block_duration);
        if (this_time_block != current_time_block) {
          if (current_time_block != -1) {
            writer.WriteTimeBlocks(time_block);
          }
          current_time_block = this_time_block;
          time_block = prd::TimeBlock();
          time_block.id = static_cast<uint32_t>(current_time_block);
        }
        time_block.prompt_events.push_back(event);

        Counts_binned++;
        Trues++;
      } else {
        Scatters++;
      }
	  } else {
      Randoms++;
    }
  }
  writer.WriteTimeBlocks(time_block);

  printf("Total Number of Coincidence Events in the ROOT file:= %llu ...\n",nentries );
  printf("Total Number of Coincidence Events registered in list-mode or sinogram format:= %lu ...\n", Counts_binned);

  return(0);
}
