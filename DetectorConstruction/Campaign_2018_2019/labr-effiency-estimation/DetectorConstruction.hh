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

#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH 1

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VUserDetectorConstruction.hh"

#include "utrConfig.h"

class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
	DetectorConstruction();
	~DetectorConstruction();

	virtual G4VPhysicalVolume *Construct();
	virtual void ConstructSDandField();

	void print_info() const;
    void placeLabr(
      G4LogicalVolume *World_Logical,
      G4String Detector_Name,
      G4ThreeVector global_coordinates,
      G4double LaBr_rt,
	  G4double LaBr_dy,
	  G4double LaBr_dz,
	  G4double LaBr_phi,
	  G4double LaBr_theta,
      
	  G4double LaBr_AngleX,
	  G4double LaBr_AngleY,
	  G4double LaBr_AngleZ,
      
	  G4double LaBr_Cu_Radius,
	  G4double LaBr_Cu_Thickness,
	  G4double LaBr_Pb_Radius,
	  G4double LaBr_Pb_Thickness,
	  G4double LaBr_Pb_Wrap_Thickness,
	  G4double LaBr_Pb_Wrap_Length
    );
    
    
    
    

private:
	G4double World_x;
	G4double World_y;
	G4double World_z;

	G4double G3_Target_To_2nd_Target;
	G4double Collimator_Entrance_To_G3_Target;
};

#endif
