//  *********************************************************************
//  * DISCLAIMER                                                        *
//  *                                                                   *
//  * Neither the authors of this software system, nor their employing  *
//  * institutes, nor the agencies providing financial support for this *
//  * work  make  any representation or  warranty, express or implied,  *
//  * regarding  this  software system or assume any liability for its  *
//  * use.                                                              *
//  *                                                                   *
//  * By copying,  distributing  or modifying the Program (or any work  *
//  * based  on  the Program)  you indicate  your  acceptance of  this  *
//  * statement, and all its terms.                                     *
//  *********************************************************************
//
//######################################################################################
//# Authors: Nicolas A Karakatsanis, Kris Thielemans, Michael Hansen, Hideaki Tashima  #
//#                                                                                    #
//# Date  23-NOV-2023 - July 2025                                                      #
//#                                                                                    #
//# Objective : To read the coincidences TTree from the .root file, and generate  the  #
//#             corresponding PETSIRD file.                                            #
//#                                                                                    #
//# Input     : Monte Carlo data from GATE                                             #
//#                                                                                    #
//# Output    : PETSIRD file (tested with v0.7)                                        #
//#                                                                                    #
//######################################################################################
//#                                                                                    #



#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "TROOT.h"
#include "TChain.h"
#include <cmath>

#include <nlohmann/json.hpp>

// PETSIRD Includes
#include "protocols.h"
#include "types.h"
#include "binary/protocols.h"
#include "petsird_helpers/create.h"
#include "petsird_helpers/geometry.h"

using namespace std ;

constexpr float speed_of_light_mm_per_ps = 0.299792458F;

#include "petsird_helpers.h"

struct ScannerGeometry
{
  /* future expansion
  int n_rings;
  int n_det;
  int s_width;
  */
  string model_name;
  int n_rsec_xy;
  int n_rsec_z;
  int n_mod_xy;
  int n_mod_z;
  int n_smod_xy;
  int n_smod_z;
  int n_cry_xy;
  int n_cry_z;
  int n_cry_layers;
  int cry_ax_gap;
  int cry_tx_gap;
  int smod_ax_gap;
  int smod_tx_gap;
  unt mod_ax_gap;
  int mod_tx_gap;
  int rsec_ax_gap;
  int rsec_tx_gap;
  /* future expansion
  int max_d_ring;
  */
  int number_of_TOF_bins;
  float TOF_bin_width_mm;
  int number_of_energy_bins;
  float radius;
  /* future expansion
  int tx_virtual_crystal_num;
  int ax_virtual_crystal_num;
  int tx_phys_crystal_num;
  int ax_phys_crystal_num;
  */
  float detector_x_dim, detector_y_dim, detector_z_dim; // in mm
  float energy_LLD, energy_ULD;
  float EnergyResolutionAt511;
  float TOF_resolution_mm;
  float LM_time_block_duration_ms;
  /* future expansion
  float ArcLength;
  float TxFOV;
  float TxFOV_TOF; // in ps
  float module_axial_pitch;
  */
};

void WriteScannerGeometry(const ScannerGeometry& scanner_geometry, const std::string& filename)
{
  nlohmann::json j;
  /* future expansion
  j["n_rings"] = scanner_geometry.n_rings;
  j["n_det"] = scanner_geometry.n_det;
  j["s_width"] = scanner_geometry.s_width;
  */
  j["model_name"] = scanner_geometry.model_name;
  j["n_rsec_xy"] = scanner_geometry.n_rsec_xy;
  j["n_rsec_z"] = scanner_geometry.n_rsec_z;
  j["n_mod_xy"] = scanner_geometry.n_mod_xy;
  j["n_mod_z"] = scanner_geometry.n_mod_z;
  j["n_smod_xy"] = scanner_geometry.n_smod_xy;
  j["n_smod_z"] = scanner_geometry.n_smod_z;
  j["n_cry_xy"] = scanner_geometry.n_cry_xy;
  j["n_cry_z"] = scanner_geometry.n_cry_z;
  j["n_cry_layers"] = scanner_geometry.n_cry_layers;
  j["cry_ax_gap"] = scanner_geometry.cry_ax_gap;
  j["cry_tx_gap"] = scanner_geometry.cry_tx_gap;
  j["smod_ax_gap"] = scanner_geometry.smod_ax_gap;
  j["smod_tx_gap"] = scanner_geometry.smod_tx_gap;
  j["mod_ax_gap"] = scanner_geometry.mod_ax_gap;
  j["mod_tx_gap"] = scanner_geometry.mod_tx_gap;
  j["rsec_ax_gap"] = scanner_geometry.rsec_ax_gap;
  j["rsec_tx_gap"] = scanner_geometry.rsec_tx_gap;
  j["number_of_TOF_bins"] = scanner_geometry.number_of_TOF_bins;
  j["TOF_bin_width_mm"] = scanner_geometry.TOF_bin_width_mm;
  j["TOF_resolution_mm"] = scanner_geometry.TOF_resolution_mm;
  j["number_of_energy_bins"] = scanner_geometry.number_of_energy_bins;
  j["radius"] = scanner_geometry.radius;
  /* fields for future expansion using binning
  j["tx_virtual_crystal_num"] = scanner_geometry.tx_virtual_crystal_num;
  j["ax_virtual_crystal_num"] = scanner_geometry.ax_virtual_crystal_num;
  j["tx_phys_crystal_num"] = scanner_geometry.tx_phys_crystal_num;
  j["ax_phys_crystal_num"] = scanner_geometry.ax_phys_crystal_num;
  */
  j["detector_x_dim"] = scanner_geometry.detector_x_dim;
  j["detector_y_dim"] = scanner_geometry.detector_y_dim;
  j["detector_z_dim"] = scanner_geometry.detector_z_dim;
  j["energy_LLD"] = scanner_geometry.energy_LLD;
  j["energy_ULD"] = scanner_geometry.energy_ULD;
  j["EnergyResolutionAt511"] = scanner_geometry.EnergyResolutionAt511;
  j["LM_time_block_duration_ms"] = scanner_geometry.LM_time_block_duration_ms;

  /* future expansion
  j["ArcLength"] = scanner_geometry.s_width * scanner_geometry.detector_y_dim / 2.0f;
  j["TxFOV"] = 2 * scanner_geometry.radius * sin (scanner_geometry.ArcLength / (2 * scanner_geometry.radius) );
  j["TxFOV_TOF"] = scanner_geometry.number_of_TOF_bins * scanner_geometry.TOF_bin_width_mm;
  j["module_axial_pitch"] = scanner_geometry.n_cry_z * scanner_geometry.detector_z_dim + (scanner_geometry.n_cry_z - 1) * scanner_geometry.cry_ax_gap;
  */

  std::ofstream o(filename);
  o << std::setw(4) << j << std::endl;
}

auto get_value(const nlohmann::json& j, const std::string& key)
{
  if (!j.contains(key))
    {
      std::stringstream s;
      s << "ERROR: key \"" << key << "\" is missing from the JSON file";
      throw std::runtime_error(s.str());
    }
  return j[key];
}

// Function for reading json scanner geometry
ScannerGeometry ReadScannerGeometry(const std::string& filename)
{
  std::ifstream i(filename);
  nlohmann::json j;
  i >> j;

  ScannerGeometry scanner_geometry;
  /* future expansion
  scanner_geometry.n_rings = get_value(j, "n_rings");
  scanner_geometry.n_det = get_value(j, "n_det");
  scanner_geometry.s_width = get_value(j, "s_width");
  */
  scanner_geometry.model_name = get_value(j, "model_name");
  scanner_geometry.n_rsec_xy = get_value(j, "n_rsec_xy");
  scanner_geometry.n_rsec_z = get_value(j, "n_rsec_z");
  scanner_geometry.n_mod_xy = get_value(j, "n_mod_xy");
  scanner_geometry.n_mod_z = get_value(j, "n_mod_z");
  scanner_geometry.n_smod_xy = get_value(j, "n_smod_xy");
  scanner_geometry.n_smod_z = get_value(j, "n_smod_z");
  scanner_geometry.n_cry_xy = get_value(j, "n_cry_xy");
  scanner_geometry.n_cry_z = get_value(j, "n_cry_z");
  scanner_geometry.n_cry_layers = get_value(j, "n_cry_layers");
  scanner_geometry.cry_ax_gap = get_value(j, "cry_ax_gap");
  scanner_geometry.cry_tx_gap = get_value(j, "cry_tx_gap");
  scanner_geometry.smod_ax_gap = get_value(j, "smod_ax_gap");
  scanner_geometry.smod_tx_gap = get_value(j, "smod_tx_gap");
  scanner_geometry.mod_ax_gap = get_value(j, "mod_ax_gap");
  scanner_geometry.mod_tx_gap = get_value(j, "mod_tx_gap");
  scanner_geometry.rsec_ax_gap = get_value(j, "rsec_ax_gap");
  scanner_geometry.rsec_tx_gap = get_value(j, "rsec_tx_gap");
  /* fields for future expansion using binning
  scanner_geometry.max_d_ring = get_value(j, "max_d_ring");
  */
  //scanner_geometry.number_of_TOF_bins = get_value(j, "number_of_TOF_bins");
  //scanner_geometry.TOF_bin_width_mm = get_value(j, "TOF_bin_width_mm");
  scanner_geometry.number_of_TOF_bins = get_value(j, "number_of_TOF_bins");
  scanner_geometry.TOF_bin_width_mm = get_value(j, "TOF_bin_width_mm");

  scanner_geometry.number_of_energy_bins = get_value(j, "number_of_energy_bins");
  scanner_geometry.radius = get_value(j, "radius");
  /* fields for future expansion using binning
  scanner_geometry.tx_virtual_crystal_num = get_value(j, "tx_virtual_crystal_num");
  scanner_geometry.ax_virtual_crystal_num = get_value(j, "ax_virtual_crystal_num");
  scanner_geometry.tx_phys_crystal_num = get_value(j, "tx_phys_crystal_num");
  scanner_geometry.ax_phys_crystal_num = get_value(j, "ax_phys_crystal_num");
  */
  scanner_geometry.detector_x_dim = get_value(j, "detector_x_dim");
  scanner_geometry.detector_y_dim = get_value(j, "detector_y_dim");
  scanner_geometry.detector_z_dim = get_value(j, "detector_z_dim");
  scanner_geometry.energy_LLD = get_value(j, "energy_LLD");
  scanner_geometry.energy_ULD = get_value(j, "energy_ULD");
  scanner_geometry.EnergyResolutionAt511 = get_value(j, "EnergyResolutionAt511");
  scanner_geometry.TOF_resolution_mm = get_value(j, "TOF_resolution_mm");
  scanner_geometry.LM_time_block_duration_ms = get_value(j, "LM_time_block_duration_ms");

  /* future expansion
  scanner_geometry.ArcLength = scanner_geometry.s_width * scanner_geometry.detector_y_dim / 2.0f;
  scanner_geometry.TxFOV = 2 * scanner_geometry.radius * sin (scanner_geometry.ArcLength / (2 * scanner_geometry.radius) );
  scanner_geometry.TxFOV_TOF = scanner_geometry.TOF_bin_width_mm * scanner_geometry.number_of_TOF_bins;
  scanner_geometry.module_axial_pitch = scanner_geometry.n_cry_z * scanner_geometry.detector_z_dim + (scanner_geometry.n_cry_z - 1) * scanner_geometry.cry_ax_gap;
  */
  return scanner_geometry;
}

void usage()
{
  std::cout << "Usage: root_to_petsird [options]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -r, --root-prefix <root_prefix>             Prefix of root files" << std::endl;
  std::cout << "  -s, --scanner-geometry-file <filename>      Scanner geometry file" << std::endl;
  std::cout << "  -c, --normalization-file <filename>         Normalization file" << std::endl;
  std::cout << "  -p, --petsird-file <filename>               PETSiRD file" << std::endl;
  std::cout << "  -n, --number-of-root-files <number>         Number of root files" << std::endl;
  std::cout << "  -v, --verbose                               Verbose output" << std::endl;
  std::cout << "  -h, --help                                  Print this help message" << std::endl;
}

int calculate_element_index(int module_id, int submodule_id, int crystal_id, int layer_id, const ScannerGeometry& scannerGeometry)
{
  //int N_DET = scannerGeometry.n_det;
  //int N_MOD_xy = scannerGeometry.n_mod_xy;
  //int N_MOD_z = scannerGeometry.n_mod_z;
  int N_SMOD_xy = scannerGeometry.n_smod_xy;
  int N_SMOD_z = scannerGeometry.n_smod_z;
  int N_CRY_xy = scannerGeometry.n_cry_xy;
  int N_CRY_z = scannerGeometry.n_cry_z;
  int N_CRY_layers = scannerGeometry.n_cry_layers;

  // NK code
  // return (Int_t)(layer_id)
  //            + (Int_t)(crystal_id/N_CRY_xy + (crystal_id%N_CRY_xy)*N_CRY_z)*N_CRY_layers
  //            + (Int_t)(submodule_id/N_SMOD_xy + (submodule_id%N_SMOD_xy)*N_SMOD_z)*N_CRY_layers*N_CRY_xy*N_CRY_z
  //            + (Int_t)(module_id/N_MOD_xy + (module_id%N_MOD_xy)*N_MOD_z)*N_CRY_layers*N_CRY_xy*N_CRY_z*N_SMOD_xy*N_SMOD_z;
  return (Int_t)((((module_id*N_SMOD_xy*N_SMOD_z)
                   + submodule_id)*N_CRY_xy*N_CRY_z
                   + crystal_id) * N_CRY_layers
                   + layer_id);
}

void calculate_scanner_layer_xyz_coordinates(unsigned int& rsec_z_id, unsigned int& rsec_xy_id,  
                                             unsigned int& mod_z_id, unsigned int& mod_xy_id,
                                             unsigned int& smod_z_id, unsigned int& smod_xy_id,
                                             unsigned int& cry_z_id, unsigned int& cry_xy_id,
                                             unsigned int& layer_id,
                                             unsigned int global_element_index,
                                             const ScannerGeometry& scannerGeometry)
{
  rsec_z_id = ( global_element_index / (scannerGeometry.n_rsec_xy*scannerGeometry.n_mod_xy*scannerGeometry.n_smod_xy*scannerGeometry.n_cry_xy*scannerGeometry.n_cry_layers)
                                     / (scannerGeometry.n_mod_z*scannerGeometry.n_smod_z*scannerGeometry.n_cry_z));
  rsec_xy_id = ( global_element_index % (scannerGeometry.n_rsec_xy*scannerGeometry.n_mod_xy*scannerGeometry.n_smod_xy*scannerGeometry.n_cry_xy*scannerGeometry.n_cry_layers)) 
                                      / (scannerGeometry.n_mod_xy*scannerGeometry.n_smod_xy*scannerGeometry.n_cry_xy*scannerGeometry.n_cry_layers);

  mod_z_id = ( global_element_index / (scannerGeometry.n_rsec_xy*scannerGeometry.n_mod_xy*scannerGeometry.n_smod_xy*scannerGeometry.n_cry_xy*scannerGeometry.n_cry_layers)) % scannerGeometry.n_mod_z;                          
  mod_xy_id = (global_element_index / (scannerGeometry.n_smod_xy*scannerGeometry.n_cry_xy*scannerGeometry.n_cry_layers)) % scannerGeometry.n_mod_xy;

  smod_z_id = ( global_element_index / (scannerGeometry.n_rsec_xy*scannerGeometry.n_mod_xy*scannerGeometry.n_smod_xy*scannerGeometry.n_cry_xy*scannerGeometry.n_cry_layers)) % scannerGeometry.n_smod_z;                          
  smod_xy_id = ( global_element_index / (scannerGeometry.n_cry_xy*scannerGeometry.n_cry_layers)) % scannerGeometry.n_smod_xy;

  cry_z_id = ( global_element_index / (scannerGeometry.n_rsec_xy*scannerGeometry.n_mod_xy*scannerGeometry.n_smod_xy*scannerGeometry.n_cry_xy*scannerGeometry.n_cry_layers)) % scannerGeometry.n_cry_z;                          
  cry_xy_id = ( global_element_index / (scannerGeometry.n_cry_layers)) % scannerGeometry.n_cry_xy;

  layer_id = global_element_index % scannerGeometry.n_cry_layers;
}

void calculate_scanner_layer_coordinates(unsigned int& rsec_id, unsigned int& mod_id, unsigned int& smod_id, unsigned int& cry_id
                                         unsigned int rsec_z_id, unsigned int rsec_xy_id, 
                                         unsigned int mod_z_id, unsigned int mod_xy_id, 
                                         unsigned int smod_z_id, unsigned int smod_xy_id, 
                                         unsigned int cry_z_id, unsigned int cry_xy_id, 
                                         const ScannerGeometry& scannerGeometry)
{
  //index order according to rings: first the structures of the 1st ring along all TX positions, then of the 2nd ring etc
  //rsec_id = rsec_xy_id + rsec_z_id*scannerGeometry.n_rsec_xy;
  //mod_id = mod_xy_id + mod_z_id*scannerGeometry.n_mod_xy;
  //smod_id = smod_xy_id + smod_z_id*scannerGeometry.n_smod_xy;
  //cry_id = cry_xy_id + cry_z_id*scannerGeometry.n_cry_xy;
	
  //index order according to TX position in each ring: first the structures of the 1st TX position along all rings, then of the 2nd TX position etc.
  //this complies with current module_pair_SGID_LUT conventions in the persird_helper petsird_generator example
  rsec_id = rsec_z_id + rsec_xy_id*scannerGeometry.n_rsec_z;
  mod_id = mod_z_id + mod_xy_id*scannerGeometry.n_mod_z;
  smod_id = smod_z_id + smod_xy_id*scannerGeometry.n_smod_z;
  cry_id = cry_z_id + cry_xy_id*scannerGeometry.n_cry_z;	
}

unsigned int calculate_module_index(unsigned int gantry_id, unsigned int rsector_id, const ScannerGeometry& scannerGeometry)
{
  const int N_RSEC_xy = scannerGeometry.n_rsec_xy;
  const int N_RSEC_z = scannerGeometry.n_rsec_z;

  // NK code
  // return (Int_t)(rsector_id/N_RSEC_xy + (rsector_id%N_RSEC_xy)*N_RSEC_z)
  //           + (Int_t)(gantry_id)*N_RSEC_xy*N_RSEC_z;
  return (Int_t)(gantry_id)*N_RSEC_z*N_RSEC_xy + (Int_t)rsector_id;
}

//! return a cuboid volume
petsird::BoxSolidVolume
get_crystal(const ScannerGeometry& scannerGeometry)
{
 //shift 0 in y and z
 using petsird::Coordinate;
  petsird::BoxShape crystal_shape{ Coordinate{ { -scannerGeometry.detector_x_dim/2, -scannerGeometry.detector_y_dim/2, -scannerGeometry.detector_z_dim/2 } },
                                   Coordinate{ { -scannerGeometry.detector_x_dim/2, -scannerGeometry.detector_y_dim/2, scannerGeometry.detector_z_dim/2 } },
                                   Coordinate{ { -scannerGeometry.detector_x_dim/2, scannerGeometry.detector_y_dim/2, scannerGeometry.detector_z_dim/2 } },
                                   Coordinate{ { -scannerGeometry.detector_x_dim/2, scannerGeometry.detector_y_dim/2, -scannerGeometry.detector_z_dim/2 } },
                                   Coordinate{ { scannerGeometry.detector_x_dim/2, -scannerGeometry.detector_y_dim/2, -scannerGeometry.detector_z_dim/2 } },
                                   Coordinate{ { scannerGeometry.detector_x_dim/2, -scannerGeometry.detector_y_dim/2, scannerGeometry.detector_z_dim/2 } },
                                   Coordinate{ { scannerGeometry.detector_x_dim/2, scannerGeometry.detector_y_dim/2, scannerGeometry.detector_z_dim/2 } },
                                   Coordinate{ { scannerGeometry.detector_x_dim/2, scannerGeometry.detector_y_dim/2, -scannerGeometry.detector_z_dim/2 } } };
  petsird::BoxSolidVolume crystal{ crystal_shape, /* material_id */ 1 };
  return crystal;
}


//! return a module of NUM_CRYSTALS_PER_MODULE cuboids
petsird::DetectorModule
get_detector_module(const ScannerGeometry& scannerGeometry)
{
  petsird::ReplicatedBoxSolidVolume rep_volume;
  {
    rep_volume.object = get_crystal(scannerGeometry);
    for (int rep_mod_z = 0; rep_mod_z < scannerGeometry.n_mod_z; ++rep_mod_z)
      for (int rep_mod_xy = 0; rep_mod_xy < scannerGeometry.n_mod_xy; ++rep_mod_xy)
        for (int rep_smod_z = 0; rep_smod_z < scannerGeometry.n_smod_z; ++rep_smod_z)
          for (int rep_smod_xy = 0; rep_smod_xy < scannerGeometry.n_smod_xy; ++rep_smod_xy)
            for (int rep_cry_z = 0; rep_cry_z < scannerGeometry.n_cry_z; ++rep_cry_z)
              for (int rep_cry_xy = 0; rep_cry_xy < scannerGeometry.n_cry_xy; ++rep_cry_xy)
                  for (int rep_cry_layer = 0; rep_cry_layer < scannerGeometry.n_cry_layers; ++rep_cry_layer)
                    {
                      petsird::RigidTransformation transform{ { { 1.F, 0.F, 0.F, scannerGeometry.radius + rep_cry_layer * scannerGeometry.detector_x_dim + scannerGeometry.detector_x_dim/2 },
                                                                { 0.F, 1.F, 0.F, (rep_mod_xy - scannerGeometry.n_mod_xy / 2 + ((scannerGeometry.n_mod_xy - 1) % 2 ) * 0.5F) * scannerGeometry.n_smod_xy * scannerGeometry.n_cry_xy * scannerGeometry.detector_y_dim
                                                                               + (rep_smod_xy - scannerGeometry.n_smod_xy / 2 + ((scannerGeometry.n_smod_xy - 1) % 2 ) * 0.5F) * scannerGeometry.n_cry_xy * scannerGeometry.detector_y_dim
                                                                               + (rep_cry_xy - scannerGeometry.n_cry_xy / 2 + ((scannerGeometry.n_cry_xy - 1) % 2 ) * 0.5F) * scannerGeometry.detector_y_dim },
                                                                { 0.F, 0.F, 1.F, (rep_mod_z - scannerGeometry.n_mod_z / 2 + ((scannerGeometry.n_mod_z - 1) % 2 ) * 0.5F) * scannerGeometry.n_smod_z * scannerGeometry.n_cry_z * scannerGeometry.detector_z_dim
                                                                               + (rep_smod_z - scannerGeometry.n_smod_z / 2 + ((scannerGeometry.n_smod_z - 1) % 2 ) * 0.5F) * scannerGeometry.n_cry_z * scannerGeometry.detector_z_dim
                                                                               + (rep_cry_z - scannerGeometry.n_cry_z / 2 + ((scannerGeometry.n_cry_z - 1) % 2 ) * 0.5F) * scannerGeometry.detector_z_dim } } };
                      rep_volume.transforms.push_back(transform);
                    }
  }

  petsird::DetectorModule detector_module;
  detector_module.detecting_elements = rep_volume;

  return detector_module;
}


//! return scanner build by rotating a module around the (0,0,1) axis
petsird::ScannerGeometry
get_scanner_geometry(const ScannerGeometry& scannerGeometry)
{
  petsird::ReplicatedDetectorModule rep_module;
  {
    rep_module.object = get_detector_module(scannerGeometry);
    std::vector<float> angles;
    for (unsigned int i = 0; i < scannerGeometry.n_rsec_xy; ++i)
      {
        angles.push_back(static_cast<float>((2 * M_PI * i) / scannerGeometry.n_rsec_xy));
      }

    for (auto angle : angles)
      for (int rep_rsec_z = 0; rep_rsec_z < scannerGeometry.n_rsec_z; ++rep_rsec_z)
        {
          petsird::RigidTransformation transform{ { { std::cos(angle), -std::sin(angle), 0.F, 0.F },
                                                    { std::sin(angle), std::cos(angle), 0.F, 0.F},
                                                    { 0.F, 0.F, 1.F, (rep_rsec_z - scannerGeometry.n_rsec_z / 2 + ((scannerGeometry.n_rsec_z - 1) % 2 ) * 0.5F) * scannerGeometry.n_mod_z * scannerGeometry.n_smod_z * scannerGeometry.n_cry_z * scannerGeometry.detector_z_dim} } };

          rep_module.transforms.push_back(transform);
        }
  }
  petsird::ScannerGeometry scanner_geometry;
  scanner_geometry.replicated_modules.push_back(rep_module);
  return scanner_geometry;
}

// set some example efficiencies in the ScannerInformation object.
void
set_detection_efficiencies(petsird::ScannerInformation& scanner, const ScannerGeometry& scannerGeometry)
{
  //const auto num_types_of_modules = scanner.scanner_geometry.replicated_modules.size();
  // only 1 type of module in the current scanner
  //assert(num_types_of_modules == 1);
  const petsird::TypeOfModule type_of_module{ 0 };
  const auto num_detection_bins = petsird_helpers::get_num_detection_bins(scanner, type_of_module);

  const auto& event_energy_bin_edges = scanner.event_energy_bin_edges[type_of_module];
  const auto num_event_energy_bins = event_energy_bin_edges.NumberOfBins();

  // set all detection_bin_efficiencies to 1 in this example
  if (scanner.detection_efficiencies.detection_bin_efficiencies)
    {
      auto& bin_effs = (*scanner.detection_efficiencies.detection_bin_efficiencies)[type_of_module];
      yardl::resize(bin_effs, { num_detection_bins });
      std::fill(begin(bin_effs), end(bin_effs), 1.F);
    }

  // check if the caller wants to have module-pair stuff. If not, return.
  if (!scanner.detection_efficiencies.module_pair_efficiencies_vectors)
    return;

  const auto& rep_module = scanner.scanner_geometry.replicated_modules[type_of_module];
  const auto num_modules = rep_module.transforms.size();

  // We will only use rotational symmetries (no translation along the axis yet)
  // We also assume all module-pairs are in coincidence, except those with the same angle.
  // Writing a module number as (z-position, angle):
  //   eff((z1,a1), (z2, a2)) == eff((z1,0), (z2, abs(a2-a1)))
  // or in linear indices
  //   eff(z1 + NZ * a1, z2 + NZ * a2) == eff(z1, z2 + NZ * abs(a2 - a1))
  // (coincident) SGIDs need to start from 0, so ignoring self-coincident angles
  const unsigned int num_SGIDs = scannerGeometry.n_rsec_z * scannerGeometry.n_rsec_z * (scannerGeometry.n_rsec_xy - 1);
  // SGID = z1 + NZ * (z2 + NZ * abs(a2 - a1) - 1)
  const auto NZ = scannerGeometry.n_rsec_z;

  auto& module_pair_SGID_LUT = (*scanner.detection_efficiencies.module_pair_sgidlut)[type_of_module][type_of_module];
  module_pair_SGID_LUT = yardl::NDArray<int, 2>({ num_modules, num_modules });
  for (unsigned int mod1 = 0; mod1 < num_modules; ++mod1)
    {
      for (unsigned int mod2 = 0; mod2 < num_modules; ++mod2)
        {
          const auto z1 = mod1 % NZ;
          const auto a1 = mod1 / NZ;
          const auto z2 = mod2 % NZ;
          const auto a2 = mod2 / NZ;
          if (a1 == a2)
            {
              module_pair_SGID_LUT(mod1, mod2) = -1;
            }
          else
            {
              module_pair_SGID_LUT(mod1, mod2) = z1 + NZ * (z2 + NZ * (std::abs(int(a2) - int(a1)) - 1));
            }
        }
    }
  // assert(module_pair_SGID_LUT).max() == num_SGIDs - 1);
  // initialise module_pair_efficiencies
  auto& module_pair_efficiencies_vector
      = (*scanner.detection_efficiencies.module_pair_efficiencies_vectors)[type_of_module][type_of_module];
  // assign an empty vector first, and reserve correct size
  module_pair_efficiencies_vector = petsird::ModulePairEfficienciesVector();
  module_pair_efficiencies_vector.reserve(num_SGIDs);

  const auto& detecting_elements = rep_module.object.detecting_elements;
  const auto num_det_els_in_module = detecting_elements.transforms.size();
  const auto num_detection_bins_in_module = num_det_els_in_module * num_event_energy_bins;
  for (unsigned int SGID = 0; SGID < num_SGIDs; ++SGID)
    {
      // extract first module_pair for this SGID. However, as this currently unused, it is commented out
      // const auto& module_pair = *std::find(module_pair_SGID_LUT.begin(), module_pair_SGID_LUT.end(), SGID);
      petsird::ModulePairEfficiencies module_pair_efficiencies;
      module_pair_efficiencies.values = yardl::NDArray<float, 2>({ num_detection_bins_in_module, num_detection_bins_in_module });
      // give some (non-physical) value
      module_pair_efficiencies.values.fill(1.F);
      module_pair_efficiencies.sgid = SGID;
      module_pair_efficiencies_vector.emplace_back(module_pair_efficiencies);
      assert(module_pair_efficiencies_vector.size() == unsigned(SGID + 1));
    }
}

void
SetEfficienciesFromFile(petsird::ScannerInformation& scanner, const ScannerGeometry& scannerGeometry, const string& filename)
{
  const petsird::TypeOfModule type_of_module{ 0 };

  const auto& rep_module = scanner.scanner_geometry.replicated_modules[type_of_module];
  //const auto num_modules = rep_module.transforms.size();
  //printf("num_modules=%ld\n", num_modules);
  //printf("scannerGeometry.n_rsec_z=%d\n", scannerGeometry.n_rsec_z);

  auto& module_pair_SGID_LUT = (*scanner.detection_efficiencies.module_pair_sgidlut)[type_of_module][type_of_module];

  const auto& detecting_elements = rep_module.object.detecting_elements;
  const auto num_det_els_in_module = detecting_elements.transforms.size();
  const auto& event_energy_bin_edges = scanner.event_energy_bin_edges[type_of_module];
  const auto num_event_energy_bins = event_energy_bin_edges.NumberOfBins();
  const auto num_detection_bins_in_module = num_det_els_in_module * num_event_energy_bins;

  const unsigned int num_SGIDs = scannerGeometry.n_rsec_z * scannerGeometry.n_rsec_z * (scannerGeometry.n_rsec_xy - 1);

  yardl::NDArray<float, 3> sum_components({num_SGIDs, num_detection_bins_in_module, num_detection_bins_in_module});
  yardl::NDArray<int, 3> count_components({num_SGIDs, num_detection_bins_in_module, num_detection_bins_in_module});

  //printf("scannerGeometry.n_smod_xy=%d, scannerGeometry.n_smod_z=%d\n", scannerGeometry.n_smod_xy, scannerGeometry.n_smod_z);
  //printf("scannerGeometry.n_cry_xy=%d, scannerGeometry.n_cry_z=%d, scannerGeometry.n_cry_layers=%d\n", scannerGeometry.n_cry_xy, scannerGeometry.n_cry_z, scannerGeometry.n_cry_layers);

  //getchar();

  std::ifstream fin(filename, std::ios::binary);
  if (!fin) {
    std::cerr << "Failed to open normalization file!\n" << std::endl;
    exit(1);
  }
  std::cout << "Reading of the normalization file ...";
  while (fin) {
    float value;
    uint32_t global_element_index1;
    uint32_t global_element_index2;
    if (fin.read(reinterpret_cast<char*>(&value), sizeof(float))
     && fin.read(reinterpret_cast<char*>(&global_element_index1), sizeof(uint32_t))
     && fin.read(reinterpret_cast<char*>(&global_element_index2), sizeof(uint32_t)))
    {
      if (value!=0) value=1./value; //invert norm correction factors to convert to detection efficiencies
      //printf("%f, %d, %d\n", value, global_element_index1, global_element_index2);
      unsigned int rsec_z_id1, rsec_xy_id1, mod_z_id1, mod_xy_id1, smod_z_id1, smod_xy_id1, cry_xy_id1, cry_z_id1, layer_id1;
      calculate_scanner_layer_xyz_coordinates (rsec_z_id1, rsec_xy_id1,
                                              mod_z_id1, mod_xy_id1, 
                                              smod_z_id1, smod_xy_id1,
                                              cry_z_id1, cry_xy_id1,
                                              layer_id1,
	                                      global_element_index1,
                                              scannerGeometry);
      unsigned int rsec_z_id2, rsec_xy_id2, mod_z_id2, mod_xy_id2, smod_z_id2, smod_xy_id2, cry_xy_id2, cry_z_id2, layer_id2;
      calculate_scanner_layer_xyz_coordinates (rsec_z_id2, rsec_xy_id2,
                                              mod_z_id2, mod_xy_id2,
                                              smod_z_id2, smod_xy_id2,
                                              cry_z_id2, cry_xy_id2,
                                              layer_id2,
	                                      global_element_index2,
                                              scannerGeometry);
      unsigned int rsec_id1, mod_id1, smod_id1, cry_id1;
      calculate_scanner_layer_coordinates(rsec_id1, mod_id1, smod_id1, cry_id1,
	                                  rsec_z_id1, rsec_xy_id1,
                                          mod_z_id1, mod_xy_id1, 
                                          smod_z_id1, smod_xy_id1, 
                                          cry_z_id1, cry_xy_id1, 
                                          scannerGeometry);
      unsigned int rsec_id2, mod_id2, smod_id2, cry_id2;
      calculate_scanner_layer_coordinates(rsec_id2, mod_id2, smod_id2, cry_id2,
				  	  rsec_z_id2, rsec_xy_id2,
				          mod_z_id2, mod_xy_id2, 
				          smod_z_id2, smod_xy_id2, 
				          cry_z_id2, cry_xy_id2, 
				          scannerGeometry);
      // int module_id2, submodule_id2, crystal_id2, layer_id2;
      // int mod2 = module_id2 * scannerGeometry.n_smod_xy * scannerGeometry.n_smod_z + submodule_id2;
      // int detection_bin2 = crystal_id2 * scannerGeometry.n_cry_layers + layer_id2;
      
      // int SGID = module_pair_SGID_LUT(mod1, mod2);
      // unsigned int mod_id1 = rsec_xy_id1 + scannerGeometry.n_rsec_xy * rsec_z_id1;
      // unsigned int mod_id2 = rsec_xy_id2 + scannerGeometry.n_rsec_xy * rsec_z_id2;
      int SGID = module_pair_SGID_LUT(rsec_id1, rsec_id2);
      // printf("module_id1=%d, submodule_id1=%d, crystal_id1=%d, layer_id1=%d\n", module_id1, submodule_id1, crystal_id1, layer_id1);
      // printf("module_id2=%d, submodule_id2=%d, crystal_id2=%d, layer_id2=%d\n", module_id2, submodule_id2, crystal_id2, layer_id2);
      // printf("mod1=%d, mod2=%d\n", mod1, mod2);
      //printf("SGID=%d, rsec_id1=%d, rsec_id2=%d\n", SGID, rsec_id1, rsec_id2);
      unsigned int element_id_in_module1 = calculate_element_index(mod_id1, smod_id1, cry_id1, layer_id1, scannerGeometry);
      unsigned int element_id_in_module2 = calculate_element_index(mod_id2, smod_id2, cry_id2, layer_id2, scannerGeometry);

      if (SGID>=0) {
        sum_components(SGID, element_id_in_module1, element_id_in_module2) += value;
        count_components(SGID, element_id_in_module1, element_id_in_module2) ++;
      }
    }
  }
  std::cout << "completed." << std::endl;
  
  auto& module_pair_efficiencies_vector
      = (*scanner.detection_efficiencies.module_pair_efficiencies_vectors)[type_of_module][type_of_module];
  for (int SGID=0; SGID < num_SGIDs; ++SGID) {
    for (unsigned int detection_bin1=0; detection_bin1<num_detection_bins_in_module; ++detection_bin1) {
      for (unsigned int detection_bin2=0; detection_bin2<num_detection_bins_in_module; ++detection_bin2) {
        if (count_components(SGID, detection_bin1, detection_bin2)>0) {
          module_pair_efficiencies_vector[SGID].values(detection_bin1, detection_bin2) 
            = sum_components(SGID, detection_bin1, detection_bin2) / count_components(SGID, detection_bin1, detection_bin2);
        } else {
          module_pair_efficiencies_vector[SGID].values(detection_bin1, detection_bin2) = 0;
        }
      }
    }
  }

}

petsird::ScannerInformation
get_scanner_info(const ScannerGeometry& scannerGeometry, bool store_det_efficiencies)
{
  const unsigned long NUMBER_OF_TOF_BINS = static_cast<unsigned long>(scannerGeometry.number_of_TOF_bins);
  const unsigned long NUMBER_OF_EVENT_ENERGY_BINS = static_cast<unsigned long>(scannerGeometry.number_of_energy_bins);
  const float energy_LLD = scannerGeometry.energy_LLD;
  const float energy_ULD =scannerGeometry.energy_ULD;

  petsird::ScannerInformation scanner_info;
  scanner_info.model_name = scannerGeometry.model_name; // "PETSIRD_GATEROOT"

  const auto num_types_of_modules = 1;
  // Pre-allocate various structures to have the correct size for num_types_of_modules
  // (We will still have to set descent values into each of these.)
  petsird_helpers::create::initialize_scanner_information_dimensions(scanner_info, num_types_of_modules,
                                                                     /* allocate_detection_bin_efficiencies = */ store_det_efficiencies,
                                                                     /* allocate_module_pair_efficiencies = */ store_det_efficiencies);

  // TODO scanner_info.bulk_materials

  // geometry
  scanner_info.scanner_geometry = get_scanner_geometry(scannerGeometry);

  // TOF and energy information
  {
    auto& all_tof_bin_edges = scanner_info.tof_bin_edges;
    auto& all_tof_resolutions = scanner_info.tof_resolution;
    auto& all_event_energy_bin_edges = scanner_info.event_energy_bin_edges;
    auto& all_event_energy_resolutions = scanner_info.energy_resolution_at_511;

    // only 1 type of module in the current scanner
    assert(num_types_of_modules == 1);
    const petsird::TypeOfModule type_of_module{ 0 };

    typedef yardl::NDArray<float, 1> FArray1D;
    // TOF info (in mm)
    FArray1D tof_bin_edges_arr;
    yardl::resize(tof_bin_edges_arr, { NUMBER_OF_TOF_BINS + 1 });
    for (std::size_t i = 0; i < tof_bin_edges_arr.size(); ++i)
      tof_bin_edges_arr[i] = (i - NUMBER_OF_TOF_BINS / 2.F) * scannerGeometry.TOF_bin_width_mm;
    const petsird::BinEdges tof_bin_edges{ tof_bin_edges_arr };
    all_tof_bin_edges[type_of_module][type_of_module] = tof_bin_edges;

    all_tof_resolutions[type_of_module][type_of_module] = scannerGeometry.TOF_resolution_mm;

    FArray1D event_energy_bin_edges_arr;
    yardl::resize(event_energy_bin_edges_arr, { NUMBER_OF_EVENT_ENERGY_BINS + 1 });
    for (std::size_t i = 0; i < event_energy_bin_edges_arr.size(); ++i)
      event_energy_bin_edges_arr[i] = energy_LLD + i * (energy_ULD - energy_LLD) / NUMBER_OF_EVENT_ENERGY_BINS;
    petsird::BinEdges event_energy_bin_edges{ event_energy_bin_edges_arr };
    all_event_energy_bin_edges[type_of_module] = event_energy_bin_edges;
    all_event_energy_resolutions[type_of_module] = scannerGeometry.EnergyResolutionAt511;    // as fraction of 511 (e.g. 0.11F)
  }

  if (store_det_efficiencies)
    set_detection_efficiencies(scanner_info, scannerGeometry); // initialize all valid efficiencies with the value of 1

  // TODO scanner_info.coincidence_policy = petsird::CoincidencePolicy::kRejectMultiples;
  scanner_info.delayed_coincidences_are_stored = false;
  scanner_info.triple_events_are_stored = false;
  return scanner_info;
}



uint32_t tofToIdx(double tofPos_mm, const petsird::ScannerInformation& scanner_info)
{
  constexpr petsird::TypeOfModule type_of_module{ 0 };
  const auto& tof_bin_edges = scanner_info.tof_bin_edges[type_of_module][type_of_module].edges;

  for (size_t i = 0; i < tof_bin_edges.size() - 1; ++i)
  {
    if (tofPos_mm >= tof_bin_edges[i] && tofPos_mm < tof_bin_edges[i+1])
    {
      return static_cast<uint32_t>(i);
    }
  }
  std::cout << "WARNING: TOF out of range: " << tofPos_mm << std::endl;
  //std::stringstream ss;
  //ss << "WARNING: TOF out of range: " << tofPos_mm;
  throw std::runtime_error("TOF out of range");
}

uint32_t energyToIdx(float energy, const petsird::ScannerInformation& scanner_info)
{
  constexpr petsird::TypeOfModule type_of_module{ 0 };
  const auto& energy_bin_edges = scanner_info.event_energy_bin_edges[type_of_module].edges;
  for (size_t i = 0; i < energy_bin_edges.size() - 1; ++i)
  {
    if (energy >= energy_bin_edges[i] && energy < energy_bin_edges[i+1])
    {
      return static_cast<uint32_t>(i);
    }
  }
  std::stringstream ss;
  ss << "WARNING: Energy out of range: " << energy;
  throw std::runtime_error(ss.str());
}

petsird::Coordinate mean_position(const petsird::BoxShape& box_shape)
{
  petsird::Coordinate mean;
  mean.c = {0, 0, 0};
  for (auto& corner : box_shape.corners)
    {
      mean.c += corner.c;
    }
  mean.c /= box_shape.corners.size();
  return mean;
}

auto create_new_time_block(const petsird::ScannerInformation& scanner)
{
  const auto num_module_types = scanner.scanner_geometry.replicated_modules.size();
  petsird::EventTimeBlock time_block;
  // allocate lists for prompt_events
  // warning: would have to do the same for delayeds once we store them
  time_block.prompt_events = petsird_helpers::create::construct_2D_nested_vector<petsird::ListOfCoincidenceEvents>(num_module_types, num_module_types);
  return time_block;
}

int main(int argc, char** argv)
{
  std::string root_prefix;
  std::string scanner_geometry_file;
  std::string normalization_file;
  std::string petsird_file;
  int number_of_root_files = 2;
  bool verbose = false;
  bool store_det_efficiencies = true;

  // Parse command line args:
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-r" || arg == "--root-prefix") {
      root_prefix = argv[++i];
    } else if (arg == "-s" || arg == "--scanner-geometry-file") {
      scanner_geometry_file = argv[++i];
    } else if (arg == "-c" || arg == "--normalization-file") {
      normalization_file = argv[++i];
    } else if (arg == "-p" || arg == "--petsird-file") {
      petsird_file = argv[++i];
    } else if (arg == "-n" || arg == "--number-of-root-files") {
      number_of_root_files = atoi(argv[++i]);
    } else if (arg == "-v" || arg == "--verbose") {
      verbose = true;
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

  if (normalization_file.empty()) {
    store_det_efficiencies=false;
  }
  // Print arguments and exit
  std::cout << "root_prefix: " << root_prefix << std::endl;
  std::cout << "scanner_geometry_file: " << scanner_geometry_file << std::endl;
  std::cout << "normalization file: " << normalization_file << std::endl;
  std::cout << "petsird_file: " << petsird_file << std::endl;

  // Read scanner geometry
  ScannerGeometry scannerGeometry;
  if (scanner_geometry_file.empty()) {
    std::cerr << "Need to specify scanner geometry" << std::endl;
    return 1;
  } else {
    try
      {
        scannerGeometry = ReadScannerGeometry(scanner_geometry_file);
      }
    catch (const std::runtime_error& e)
      {
        std::cerr << e.what() << std::endl;
        return 1;
      }
  }

  string filedir, inputfilename;
  Int_t   Trues = 0, Scatters = 0, Randoms = 0;

  //####################################################################
  //#             Declaration of leaves types - TTree Coincidences     #
  //####################################################################
  Float_t         			axialPos, rotationAngle, sinogramS, sinogramTheta;
  Char_t          			comptVolName1[255], comptVolName2[255];
  Int_t           			compton1, compton2, gantryID1, gantryID2;
  Int_t           			runID, sourceID1, sourceID2, eventID1, eventID2;
  Int_t           			layerID1, layerID2, crystalID1, crystalID2;
  Int_t           			submoduleID1, submoduleID2, moduleID1, moduleID2, rsectorID1, rsectorID2;
  Int_t           			comptonPhantom1, comptonPhantom2;
  Float_t         			energy1, energy2; //in MeV
  Float_t         			globalPosX1, globalPosX2, globalPosY1, globalPosY2, globalPosZ1, globalPosZ2; //in mm
  Float_t         			sourcePosX1, sourcePosX2, sourcePosY1, sourcePosY2, sourcePosZ1, sourcePosZ2; //in mm
  Double_t        			time1, time2; //in sec
  unsigned long long int   	nentries;

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
  petsird::Header header;
  header.scanner = get_scanner_info(scannerGeometry, store_det_efficiencies);

  if (store_det_efficiencies) {
    SetEfficienciesFromFile(header.scanner, scannerGeometry, normalization_file);
  }
	
  auto& scanner = header.scanner;

  if (verbose) {
    // Print scanner information
    std::cout << "Scanner information:" << std::endl;
    constexpr petsird::TypeOfModule type_of_module{ 0 };
    const auto& tof_bin_edges = header.scanner.tof_bin_edges[type_of_module][type_of_module];
    const auto num_tof_bins = tof_bin_edges.NumberOfBins();
    std::cout << "Number of TOF bins: " << num_tof_bins << std::endl;
    std::cout << "TOF bin edges (mm): " << tof_bin_edges.edges << std::endl;
    const auto TOF_resolution_mm = header.scanner.tof_resolution[type_of_module][type_of_module];
    std::cout << "TOF resolution : " << TOF_resolution_mm << " mm = " << TOF_resolution_mm / (speed_of_light_mm_per_ps / 2) << " ps\n";
    const auto& event_energy_bin_edges = header.scanner.event_energy_bin_edges[type_of_module];
    const auto num_event_energy_bins = event_energy_bin_edges.NumberOfBins();
    std::cout << "Number of energy bins: " << num_event_energy_bins << std::endl;
    std::cout << "Event energy bin edges: " << event_energy_bin_edges.edges << std::endl;
  }

  // Write PETSiRD file
  petsird::binary::PETSIRDWriter writer(petsird_file);
  writer.WriteHeader(header);

  long current_time_block = -1;
  auto time_block = create_new_time_block(header.scanner);

  const auto event_time_block_duration = scannerGeometry.LM_time_block_duration_ms; // ms
  const petsird::TypeOfModule type_of_module{ 0 };
  unsigned long Counts_binned = 0;
  for (unsigned long long int i = 0 ; i < nentries ; i++)
  {
    if (i % 1000000 == 0) {
      printf("Processing event %llu of %llu, (%f percent)\n", i, nentries, 100.0f*i/nentries);
    }

    Coincidences->GetEntry(i);
    if (eventID1 == eventID2)
    {
	    if (comptonPhantom1 == 0 && comptonPhantom2 == 0) {
        petsird::CoincidenceEvent event;
        petsird::ExpandedDetectionBin expanded_detection_bin;
        expanded_detection_bin.module_index = calculate_module_index(gantryID1, rsectorID1, scannerGeometry);
        expanded_detection_bin.element_index = calculate_element_index(moduleID1, submoduleID1, crystalID1, layerID1, scannerGeometry);
        expanded_detection_bin.energy_index = static_cast<uint32_t>(energyToIdx(1.0e3*energy1, scanner));
        event.detection_bins[0] = petsird_helpers::make_detection_bin(header.scanner, type_of_module, expanded_detection_bin);
        expanded_detection_bin.module_index = calculate_module_index(gantryID2, rsectorID2, scannerGeometry);
        expanded_detection_bin.element_index = calculate_element_index(moduleID2, submoduleID2, crystalID2, layerID2, scannerGeometry);
        expanded_detection_bin.energy_index = static_cast<uint32_t>(energyToIdx(1.0e3*energy2, scanner));
        event.detection_bins[1] = petsird_helpers::make_detection_bin(header.scanner, type_of_module, expanded_detection_bin);
        const double tofPos_mm = 1.0e12 * (time1 - time2) * speed_of_light_mm_per_ps / 2; //in mm
        event.tof_idx = tofToIdx(tofPos_mm, scanner);

        if (verbose && i%100000 == 0) {
          std::cout << "Event " << i << std::endl;
          const auto expanded_detection_bin0
            = petsird_helpers::expand_detection_bin(header.scanner, type_of_module, event.detection_bins[0]);
          const auto expanded_detection_bin1
            = petsird_helpers::expand_detection_bin(header.scanner, type_of_module, event.detection_bins[1]);
          std::cout << "    "
                    << "[ExpandedDetectionBin(module=" << expanded_detection_bin0.module_index << ", "
                    << "el=" << expanded_detection_bin0.element_index << ", "
                    << "energy_index=" << expanded_detection_bin0.energy_index
                    << "), ExpandedDetectionBin(module=" << expanded_detection_bin1.module_index << ", "
                    << "el=" << expanded_detection_bin1.element_index << ", "
                    << "energy_index=" << expanded_detection_bin1.energy_index << ")]\n";
          std::cout << "  tof_idx: " << event.tof_idx << std::endl;
          std::cout << "gantryID1=" << gantryID1 << ", rsectorID1=" << rsectorID1
                    << ", moduleID1=" << moduleID1 << ", submoduleID1=" << submoduleID1 << ", crystalID1=" << crystalID1 << '\n';
          std::cout << "gantryID2=" << gantryID2 << ", rsectorID2=" << rsectorID2
                    << ", moduleID2=" << moduleID2 << ", submoduleID2=" << submoduleID2 << ", crystalID2=" << crystalID2 << '\n';
          const auto box_shape0 = petsird_helpers::geometry::get_detecting_box(header.scanner, type_of_module, expanded_detection_bin0);
          const auto mean_pos0 = mean_position(box_shape0);
          const auto box_shape1 = petsird_helpers::geometry::get_detecting_box(header.scanner, type_of_module, expanded_detection_bin1);
          const auto mean_pos1 = mean_position(box_shape1);
          std::cout << "  pos 1           : " << mean_pos0.c[0] << ", " << mean_pos0.c[1] << ", " << mean_pos0.c[2] << "\n";
          std::cout << "  pos 2           : " << mean_pos1.c[0] << ", " << mean_pos1.c[1] << ", " << mean_pos1.c[2] << "\n";

          std::cout << "  GlobalPosition 1: " << globalPosX1 << ", " << globalPosY1 << ", " << globalPosZ1 << std::endl;
          float distance_1 = std::sqrt(std::pow(mean_pos0.c[0]-globalPosX1, 2) + std::pow(mean_pos0.c[1]-globalPosY1, 2) + std::pow(mean_pos0.c[2]-globalPosZ1, 2));
          std::cout << "  Distance 1: " << distance_1 << std::endl;
          std::cout << "  GlobalPosition 2: " << globalPosX2 << ", " << globalPosY2 << ", " << globalPosZ2 << std::endl;
          float distance_2 = std::sqrt(std::pow(mean_pos1.c[0]-globalPosX2, 2) + std::pow(mean_pos1.c[1]-globalPosY2, 2) + std::pow(mean_pos1.c[2]-globalPosZ2, 2));
          std::cout << "  Distance 2: " << distance_2 << std::endl;
        }
        long this_time_block = static_cast<long>(time1*1.0e3 / event_time_block_duration);
        if (this_time_block != current_time_block) {
          if (current_time_block != -1) {
            writer.WriteTimeBlocks(time_block);
          }
          current_time_block = this_time_block;
          time_block = create_new_time_block(header.scanner);
          time_block.time_interval.start = time1*1.0e3;
          time_block.time_interval.stop = time1*1.0e3 + event_time_block_duration;
        }
        time_block.prompt_events[type_of_module][type_of_module].push_back(event);
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
  writer.EndTimeBlocks();
  writer.Close();

  printf("Total Number of Coincidence Events in the ROOT file:= %llu ...\n",nentries );
  printf("Total Number of Coincidence Events registered in list-mode or sinogram format:= %lu ...\n", Counts_binned);  return(0);
}
