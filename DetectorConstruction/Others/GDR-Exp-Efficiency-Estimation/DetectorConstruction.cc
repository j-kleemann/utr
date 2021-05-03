/*
utr - Geant4 simulation of the UTR at HIGS
Copyright (C) 2017 the developing team (see README.md)

This file is part of utr.

utr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

utr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with utr.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "DetectorConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// Materials
#include "G4NistManager.hh"
#include "Units.hh"
#include "Materials.hh"

// Geometry
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

// Detectors
#include "HPGe_Clover.hh"
#include "HPGe_Collection.hh"
#include "LaBr_3x3.hh"


// Sensitive Detectors
#include "G4SDManager.hh"
#include "EnergyDepositionSD.hh"
#include "ParticleSD.hh"
#include "SecondarySD.hh"


#include <string>
using std::string;

#include <vector>
using std::vector;

#include <array>
using std::array;




double clover_distance = 8. * 25.4 * mm;
double labr_distance = 70 * mm;

struct DetectorPosition {
  const string id;
  const double theta;
  const double phi;
  const double distance;
  const double intrinsic_rotation_angle;
};

array<DetectorPosition, 6> clover_positions{
    DetectorPosition{"clover_1", 135    * degree,   0 * degree, clover_distance, 0.5 * pi},
    DetectorPosition{"clover_2", 135    * degree, 180 * degree, clover_distance, 1.5 * pi},
    DetectorPosition{"clover_3", 125.26 * degree,  45 * degree, clover_distance, 0.5 * pi},
    DetectorPosition{"clover_6", 125.26 * degree, 315 * degree, clover_distance, 0.},
    DetectorPosition{"clover_4", 125.26 * degree, 135 * degree, clover_distance, 0.},
    DetectorPosition{"clover_5", 125.26 * degree, 225 * degree, clover_distance, 0.5 * pi},
};

array<DetectorPosition, 4> labr_positions{
    DetectorPosition{"labr_1",  90 * degree,   0 * degree, labr_distance, 0},
    DetectorPosition{"labr_2",  90 * degree,  90 * degree, labr_distance, 0},
    DetectorPosition{"labr_3",  90 * degree, 180 * degree, labr_distance, 0},
    DetectorPosition{"labr_4",  90 * degree, 270 * degree, labr_distance, 0},
};

DetectorConstruction::DetectorConstruction() {}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume *DetectorConstruction::Construct()
{

  G4Colour white(1.0, 1.0, 1.0);
  G4Colour grey(0.5, 0.5, 0.5);
  G4Colour black(0.0, 0.0, 0.0);
  G4Colour red(1.0, 0.0, 0.0);
  G4Colour green(0.0, 1.0, 0.0);
  G4Colour blue(0.0, 0.0, 1.0);
  G4Colour cyan(0.0, 1.0, 1.0);
  G4Colour magenta(1.0, 0.0, 1.0);
  G4Colour yellow(1.0, 1.0, 0.0);
  G4Colour orange(1.0, 0.5, 0.0);
  G4Colour light_orange(1.0, 0.82, 0.36);

  G4NistManager *nist = G4NistManager::Instance();
  G4Material *air = nist->FindOrBuildMaterial("G4_AIR");


  /***************** WORLD *****************/
  G4double World_x = 3000. * mm;
  G4double World_y = 3150. * mm;
  G4double World_z = 8000. * mm;

  G4Box *World_Solid = new G4Box("World_Solid", World_x * 0.5, World_y * 0.5, World_z * 0.5);

  G4LogicalVolume *World_Logical = new G4LogicalVolume(World_Solid, air, "World_Logical", 0, 0, 0);

  G4VisAttributes world_vis = G4VisAttributes(red);
  world_vis.SetForceWireframe(true);
  World_Logical->SetVisAttributes(world_vis);

  G4VPhysicalVolume *World_Physical = new G4PVPlacement(0, G4ThreeVector(), World_Logical, "World", 0, false, 0);


  /***************** FLOOR *****************/
  G4Material *concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
  G4double Floor_Thickness = 26.8 * cm; // Estimated from drawings

  // Construct Floor between UTR and planet Earth
  G4Box *Room_Floor_Solid = new G4Box("Room_Floor_Solid", World_x * 0.5, Floor_Thickness * 0.5, World_z * 0.5);
  G4LogicalVolume *Room_Floor_Logical = new G4LogicalVolume(Room_Floor_Solid, concrete, "Room_Floor_Logical");
  Room_Floor_Logical->SetVisAttributes(grey);
  new G4PVPlacement(0, G4ThreeVector(0., (-World_y + Floor_Thickness) * 0.5, 0), Room_Floor_Logical, "Room_Floor", World_Logical, false, 0, false);


  /***************** BEAMPIPE *****************/
  G4double Beampipe_Inner_Radius = .5 * 45.2 * mm;
  G4double Beampipe_Outer_Radius = .5 * 51.1 * mm;
  G4double Beampipe_Length = 1000. * mm;
  G4double relative_air_density = 1e-2;
  G4double beampipe_exit_window_thickness = 2.*mm;

  // Beam pipe
  G4Material *plexiglass = nist->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4Tubs *Beampipe_Solid = new G4Tubs("Beampipe_Solid", 0., Beampipe_Outer_Radius, Beampipe_Length*0.5, 0., twopi);
  G4LogicalVolume *Beampipe_Logical = new G4LogicalVolume(Beampipe_Solid, plexiglass, "Beampipe_Logical");
  Beampipe_Logical->SetVisAttributes(G4Colour(1,1,1, 0.8));
  new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Beampipe_Logical, "Beampipe", World_Logical, false, 0, false);

  // Beam pipe vacuum
  G4double density_of_air = 1.225e-3 * g / cm3; // Density of air at sea level and 288K, Wikipedia
  G4double density = relative_air_density*density_of_air;
  G4Material *beampipe_vacuum = new G4Material("beampipe_vacuum", density, 4); 
  G4Element *nat_O = nist->FindOrBuildElement("O");
  G4Element *nat_N = nist->FindOrBuildElement("N");
  G4Element *nat_C = nist->FindOrBuildElement("C");
  G4Element *nat_Ar = nist->FindOrBuildElement("Ar");
  beampipe_vacuum->AddElement(nat_O, 23.1781 * perCent);
  beampipe_vacuum->AddElement(nat_N, 75.5268 * perCent);
  beampipe_vacuum->AddElement(nat_Ar, 1.2827 * perCent);
  beampipe_vacuum->AddElement(nat_C,  0.0124 * perCent);

  G4Tubs *Beampipe_Vacuum_Solid = new G4Tubs("Beampipe_Vacuum_Solid", 0., Beampipe_Inner_Radius, (Beampipe_Length-2*beampipe_exit_window_thickness)*0.5, 0., twopi);
  G4LogicalVolume *Beampipe_Vacuum_Logical = new G4LogicalVolume(Beampipe_Vacuum_Solid, beampipe_vacuum, "Beampipe_Vacuum_Logical");
  Beampipe_Vacuum_Logical->SetVisAttributes(red);
  new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Beampipe_Vacuum_Logical, "Vacuum", Beampipe_Logical, false, 0, false);

  // Target ring
  G4double Target_Ring_Length = 1.5*inch; // Estimated
  G4double Target_Ring_Outer_Radius = Beampipe_Inner_Radius;
  G4double Target_Ring_Inner_Radius = Target_Ring_Outer_Radius - 2.*mm;
  G4Tubs *Target_Ring_Solid = new G4Tubs("Target_Ring_Solid", Target_Ring_Inner_Radius, Target_Ring_Outer_Radius, Target_Ring_Length*0.5, 0., twopi);
  G4LogicalVolume *Target_Ring_Logical = new G4LogicalVolume(Target_Ring_Solid, plexiglass, "Target_Ring_Logical");
  Target_Ring_Logical->SetVisAttributes(orange);
  new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Target_Ring_Logical, "Target_Ring", Beampipe_Vacuum_Logical, false, 0, false);


  /***************** PROXY TARGET *****************/
  G4double Target_Density = 8 * g / cm3; // A bit lower than in pure form of 8.35 g/cm³
  G4double Target_Mass = 3 * g; // Actually more, because this value is only elemental weight
  G4double Target_Radius = 6 * mm;
  G4double Target_Length = Target_Mass / (Target_Density * pi * pow(Target_Radius,2));

  G4Isotope *Sm154 = new G4Isotope("154Sm", 62, 154);
  G4Element *pure_Sm154 = new G4Element("pure 154Sm", "154Sm", 1);
  pure_Sm154->AddIsotope(Sm154, 1.);
  G4Material *target_material = new G4Material("target_material", 8 * g / cm3, 2); 
  target_material->AddElement(pure_Sm154, 2); //Sm(2)O(3)
  // G4Element *nat_O = nist->FindOrBuildElement("O");
  target_material->AddElement(nat_O, 3);

  G4Tubs *Target_Solid = new G4Tubs("Target_Solid", 0., Target_Radius, Target_Length*0.5, 0., twopi);
  G4LogicalVolume *Target_Logical = new G4LogicalVolume(Target_Solid, target_material, "Target_Logical");
  Target_Logical->SetVisAttributes(orange);
  new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Target_Logical, "Target", Beampipe_Vacuum_Logical, false, 0, false);



  /***************** DETECTORS *****************/
  HPGe_Collection hpge_Collection;
  vector<HPGe_Clover> clovers;
  for (auto det_pos : clover_positions) {
    clovers.push_back(HPGe_Clover(World_Logical, det_pos.id));
    clovers[clovers.size() - 1].setProperties(hpge_Collection.HPGe_Clover_Yale);
    clovers[clovers.size() - 1].useDewar();
    clovers[clovers.size() - 1].Construct(G4ThreeVector(), det_pos.theta, det_pos.phi, det_pos.distance, det_pos.intrinsic_rotation_angle);
  }

  vector<LaBr_3x3> labrs;
  for (auto det_pos : labr_positions) {
    labrs.push_back(LaBr_3x3(World_Logical, det_pos.id));
    labrs[labrs.size() - 1].useHousing();
    labrs[labrs.size() - 1].Construct(G4ThreeVector(), det_pos.theta, det_pos.phi, det_pos.distance);
  }

  // TODO: Add some filters and wrapping?
	// labr1.useFilterCase();
	// labr1.useFilterCaseRing();
	// labr1.useHousing();
	// labr1.Add_Filter("G4_Cu", 1.*1.15*mm, 45.*mm);
	// labr1.Add_Filter("G4_Pb", 1.*1.2*mm, 45.*mm);
	// labr1.Add_Wrap("G4_Pb", 1.*1.2*mm); // Estimated
	// labr1.Construct(G4ThreeVector(), det_pos.theta, det_pos.phi, det_pos.distance);


  // #ifdef USE_ZERODEGREE
  //   zeroDegree_Setup.Construct(G4ThreeVector(0., 0., G3_Target_To_2nd_Target + ZeroDegree_To_2nd_Target));
  // #endif

  // #ifdef USE_TARGETS
  //   /***************** G3_TARGET *****************/
  //   g3_Target.Set_Containing_Volume(beampipe_Long.Get_Beampipe_Vacuum());
  //   g3_Target.Construct(G4ThreeVector(0., 0., -beampipe_Long.Get_Z_Axis_Offset_Z()));
  // #endif

  return World_Physical;
}

void DetectorConstruction::ConstructSDandField()
{
  int detIDNo=0;
  for (auto det_pos : labr_positions) {
      detIDNo++;
      auto detName = det_pos.id;
      EnergyDepositionSD *sensitiveDet = new EnergyDepositionSD(detName, detName);
      G4SDManager::GetSDMpointer()->AddNewDetector(sensitiveDet);
      sensitiveDet->SetDetectorID(detIDNo);
      SetSensitiveDetector(detName, sensitiveDet, true);
  }
  for (auto det_pos : clover_positions) {
    for (int subCrystalNo=1; subCrystalNo < 5; subCrystalNo++) {
      detIDNo++;
      auto detName = det_pos.id + "_" + std::to_string(subCrystalNo);
      EnergyDepositionSD *sensitiveDet = new EnergyDepositionSD(detName, detName);
      G4SDManager::GetSDMpointer()->AddNewDetector(sensitiveDet);
      sensitiveDet->SetDetectorID(detIDNo);
      SetSensitiveDetector(detName, sensitiveDet, true);
    }
  }
  Max_Sensitive_Detector_ID = detIDNo;

  /********* ZeroDegree detector *******/

  // #ifdef USE_ZERODEGREE
  //   EnergyDepositionSD *ZeroDegreeSD = new EnergyDepositionSD("ZeroDegree", "ZeroDegree");
  //   G4SDManager::GetSDMpointer()->AddNewDetector(ZeroDegreeSD);
  //   ZeroDegreeSD->SetDetectorID(0);
  //   SetSensitiveDetector("ZeroDegree", ZeroDegreeSD, true);
  // #endif

  // EnergyDepositionSD *HPGe1SD = new EnergyDepositionSD("HPGe1", "HPGe1");
  // G4SDManager::GetSDMpointer()->AddNewDetector(HPGe1SD);
  // HPGe1SD->SetDetectorID(1);
  // SetSensitiveDetector("HPGe1", HPGe1SD, true);
}
