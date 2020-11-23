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

// Materials
#include "G4Material.hh"
#include "G4NistManager.hh"

// Geometry
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

// Sensitive Detectors
#include "EnergyDepositionSD.hh"
#include "G4SDManager.hh"
#include "ParticleSD.hh"
#include "SecondarySD.hh"

// Units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

// Colors
#include "G4Color.hh"
#include "NamedColors.hh"

// Config for target_rotation_angle by using zerodegree_offset cmake variable
#include "utrConfig.h"


DetectorConstruction::DetectorConstruction() {}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume *DetectorConstruction::Construct() {
	/***************** Geometric and Physical Properties *****************/

	const double target_rotation_angle = zerodegree_offset*deg;
	const int target_segments = 4;
	const double target_total_length = 18*mm;
	const double target_piece_length = target_total_length/target_segments;
	const double target_radius = 10/2*mm;
	const double target_container_length = 20*mm;
	const double target_container_radius = 12/2*mm;
	const double world_x = 2 * target_container_length;
	const double world_y = 2  * 2 * target_container_radius;
	const double world_z = 2 * target_container_length;

	const double sn_density = 7.265*g/pow(cm, 3);

	const double sn112_enrichment = 95.1*perCent;
	const double sn116_enrichment = 97.8*perCent;



	G4NistManager *nist = G4NistManager::Instance();

	/***************** World Volume *****************/
	G4Material *vacuum = nist->FindOrBuildMaterial("G4_Galactic");

	G4Box *world_solid = new G4Box("world_solid", world_x/2, world_y/2, world_z/2);
	G4LogicalVolume *world_logical = new G4LogicalVolume(world_solid, vacuum, "world_logical");
	G4VPhysicalVolume *world_physical = new G4PVPlacement(0, G4ThreeVector(), world_logical, "world", 0, false, 0);
	// world_logical->SetVisAttributes(G4VisAttributes::GetInvisible());
	world_logical->SetVisAttributes(G4Color(1, 1, 1, 0.1));


	/***************** Tin Target Material *****************/
	// Construct tin isotopes
	G4Isotope *sn112 = new G4Isotope("Sn112_isotope", 50, 112);
	G4Isotope *sn114 = new G4Isotope("Sn114_isotope", 50, 114);
	G4Isotope *sn115 = new G4Isotope("Sn115_isotope", 50, 115);
	G4Isotope *sn116 = new G4Isotope("Sn116_isotope", 50, 116);
	G4Isotope *sn117 = new G4Isotope("Sn117_isotope", 50, 117);
	G4Isotope *sn118 = new G4Isotope("Sn118_isotope", 50, 118);
	G4Isotope *sn119 = new G4Isotope("Sn119_isotope", 50, 119);
	G4Isotope *sn120 = new G4Isotope("Sn120_isotope", 50, 120);
	G4Isotope *sn122 = new G4Isotope("Sn122_isotope", 50, 122);
	G4Isotope *sn124 = new G4Isotope("Sn124_isotope", 50, 124);

	// Construct natural tin element (isotopic abundances from Wikipedia)
	G4Element *sn_nat_element = new G4Element("Sn_nat_element", "Sn_nat", 10);
	sn_nat_element->AddIsotope(sn112,  0.97*perCent);
	sn_nat_element->AddIsotope(sn114,  0.65*perCent);
	sn_nat_element->AddIsotope(sn115,  0.34*perCent);
	sn_nat_element->AddIsotope(sn116, 14.53*perCent);
	sn_nat_element->AddIsotope(sn117,  7.68*perCent);
	sn_nat_element->AddIsotope(sn118, 24.23*perCent);
	sn_nat_element->AddIsotope(sn119,  8.59*perCent);
	sn_nat_element->AddIsotope(sn120, 32.59*perCent);
	sn_nat_element->AddIsotope(sn122,  4.63*perCent);
	sn_nat_element->AddIsotope(sn124,  5.79*perCent);

	// Construct pure 112-tin and pure 116-tin "elements"
	G4Element *sn112_element = new G4Element("Sn112_element", "Sn112", 1);
	sn112_element->AddIsotope(sn112, 100.*perCent);
	G4Element *sn116_element = new G4Element("Sn116_element", "Sn116", 1);
	sn116_element->AddIsotope(sn116, 100.*perCent);

	// Construct enriched target material from natural and pure "elements"
	G4Material *sn112_target_material = new G4Material("Sn112_target_material", sn_density, 2);
	sn112_target_material->AddElement( sn112_element, sn112_enrichment);
	sn112_target_material->AddElement(sn_nat_element, (1.-sn112_enrichment));

	G4Material *sn116_target_material = new G4Material("Sn116_target_material", sn_density, 2);
	sn116_target_material->AddElement( sn116_element, sn116_enrichment);
	sn116_target_material->AddElement(sn_nat_element, (1.-sn116_enrichment));


	/***************** Tin Target Container *****************/
	G4Material *pvc = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
	G4Tubs *target_container_solid = new G4Tubs("target_container_solid", 0, target_container_radius, target_container_length/2, 0, twopi);
	G4LogicalVolume *target_container_logical = new G4LogicalVolume(target_container_solid, pvc, "target_container_logical");
	G4VPhysicalVolume *target_container_physical = new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0,1,0), target_rotation_angle), G4ThreeVector(), target_container_logical, "target_container", world_logical, false, 0);
	target_container_logical->SetVisAttributes(orange);


	/***************** Tin Target Pieces *****************/
	G4Tubs *target_piece_solid = new G4Tubs("target_piece_solid", 0, target_radius, target_piece_length/2, 0, twopi);
	const G4double target_pos_offset = target_piece_length/2-target_total_length/2;

	G4LogicalVolume *target_sn116_piece_1_logical = new G4LogicalVolume(target_piece_solid, sn116_target_material, "target_sn116_piece_1_logical");
	G4VPhysicalVolume *target_sn116_piece_1_physical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0*target_piece_length+target_pos_offset), target_sn116_piece_1_logical, "target_sn116_piece_1", target_container_logical, false, 0);
	G4LogicalVolume *target_sn116_piece_2_logical = new G4LogicalVolume(target_piece_solid, sn116_target_material, "target_sn116_piece_2_logical");
	G4VPhysicalVolume *target_sn116_piece_2_physical = new G4PVPlacement(0, G4ThreeVector(0, 0, 2*target_piece_length+target_pos_offset), target_sn116_piece_2_logical, "target_sn116_piece_2", target_container_logical, false, 0);
	target_sn116_piece_1_logical->SetVisAttributes(G4Color(0.2, 0.2, 0.2));
	target_sn116_piece_2_logical->SetVisAttributes(G4Color(0.6, 0.6, 0.6));


	G4LogicalVolume *target_sn112_piece_1_logical = new G4LogicalVolume(target_piece_solid, sn112_target_material, "target_sn112_piece_1_logical");
	G4VPhysicalVolume *target_sn112_piece_1_physical = new G4PVPlacement(0, G4ThreeVector(0, 0, 1*target_piece_length+target_pos_offset), target_sn112_piece_1_logical, "target_sn112_piece_1", target_container_logical, false, 0);
	G4LogicalVolume *target_sn112_piece_2_logical = new G4LogicalVolume(target_piece_solid, sn112_target_material, "target_sn112_piece_2_logical");
	G4VPhysicalVolume *target_sn112_piece_2_physical = new G4PVPlacement(0, G4ThreeVector(0, 0, 3*target_piece_length+target_pos_offset), target_sn112_piece_2_logical, "target_sn112_piece_2", target_container_logical, false, 0);
	target_sn112_piece_1_logical->SetVisAttributes(G4Color(0.4, 0.4, 0.4));
	target_sn112_piece_2_logical->SetVisAttributes(G4Color(0.8, 0.8, 0.8));

	return world_physical;
}

void DetectorConstruction::ConstructSDandField() {
	EnergyDepositionSD *target_sn116_piece_1_energyDepositionSD = new EnergyDepositionSD("target_sn116_piece_1_sd", "target_sn116_piece_1_sd");
	G4SDManager::GetSDMpointer()->AddNewDetector(target_sn116_piece_1_energyDepositionSD);
	target_sn116_piece_1_energyDepositionSD->SetDetectorID(1);
	SetSensitiveDetector("target_sn116_piece_1_logical", target_sn116_piece_1_energyDepositionSD, true);

	EnergyDepositionSD *target_sn112_piece_1_energyDepositionSD = new EnergyDepositionSD("target_sn112_piece_1_sd", "target_sn112_piece_1_sd");
	G4SDManager::GetSDMpointer()->AddNewDetector(target_sn112_piece_1_energyDepositionSD);
	target_sn112_piece_1_energyDepositionSD->SetDetectorID(2);
	SetSensitiveDetector("target_sn112_piece_1_logical", target_sn112_piece_1_energyDepositionSD, true);

	EnergyDepositionSD *target_sn116_piece_2_energyDepositionSD = new EnergyDepositionSD("target_sn116_piece_2_sd", "target_sn116_piece_2_sd");
	G4SDManager::GetSDMpointer()->AddNewDetector(target_sn116_piece_2_energyDepositionSD);
	target_sn116_piece_2_energyDepositionSD->SetDetectorID(3);
	SetSensitiveDetector("target_sn116_piece_2_logical", target_sn116_piece_2_energyDepositionSD, true);
	
	EnergyDepositionSD *target_sn112_piece_2_energyDepositionSD = new EnergyDepositionSD("target_sn112_piece_2_sd", "target_sn112_piece_2_sd");
	G4SDManager::GetSDMpointer()->AddNewDetector(target_sn112_piece_2_energyDepositionSD);
	target_sn112_piece_2_energyDepositionSD->SetDetectorID(4);
	SetSensitiveDetector("target_sn112_piece_2_logical", target_sn112_piece_2_energyDepositionSD, true);
}
