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
#include "DetectorConstructionConfig.hh"
#include "utrConfig.h"


#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

// Materials
#include "G4NistManager.hh"
#include "Materials.hh"
#include "Units.hh"

// Geometry
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

// Detectors
#include "HPGe_Clover.hh"
#include "HPGe_Coaxial.hh"
#include "HPGe_Collection.hh"

// Sensitive Detectors
#include "EnergyDepositionSD.hh"
#include "G4SDManager.hh"
#include "ParticleSD.hh"
#include "SecondarySD.hh"

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <array>
using std::array;

#include <algorithm>

DetectorConstruction::DetectorConstruction() {}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume *DetectorConstruction::Construct() {

const G4Colour invisible(1.0, 1.0, 1.0, 0.);
const G4Colour white(1.0, 1.0, 1.0);
const G4Colour lightGrey(0.75, 0.75, 0.75);
const G4Colour grey(0.5, 0.5, 0.5);
const G4Colour transparentGrey(0.5, 0.5, 0.5, 0.5);
const G4Colour darkGrey(0.25, 0.25, 0.25);
const G4Colour black(0.0, 0.0, 0.0);
const G4Colour red(1.0, 0.0, 0.0);
const G4Colour green(0.0, 1.0, 0.0);
const G4Colour blue(0.0, 0.0, 1.0);
const G4Colour cyan(0.0, 1.0, 1.0);
const G4Colour magenta(1.0, 0.0, 1.0);
const G4Colour yellow(1.0, 1.0, 0.0);
const G4Colour orange(1.0, 0.5, 0.0);

G4NistManager *nist = G4NistManager::Instance();

// --------------- General ---------------

const auto subtractionSolidBuffer = 10. * mm; // G4SubtractionSolid doesn't always fully cut solids when the cutting solid is equal in dimension to the to be cut solid, hence, add some buffer

//------Room-------

const auto activationRoomLength = 1500. * mm; // Measured off utr drawings
const auto activationRoomHeight = 1000. * mm;
const auto activationRoomWidth = 500. * mm;
const auto activationRoomFloorThickness = 200. * mm;
const auto activationSideWallThickness = 150. * mm;

//---------THREE SCREW HOLDER-----

const auto ThreeScrewHolderBottomWidth = 1.397 * mm;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 
const auto ThreeScrewHolderBottomLength = ThreeScrewHolderBottomWidth + 5.461 * mm; // From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026
const auto ThreeScrewHolderBottomAirWidth = 5.461 * mm;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 
const auto ThreeScrewHolderTopLength = 6.8072 * mm; //From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 
//const auto ThreeScrewHolderTopWidth = 1.6002 * mm;
const auto ThreeScrewHolderGrooveWidth = 1.6002 * mm;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 
const auto ThreeScrewHolderGrooveOuterDiameter = 38.1 *mm;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 
const auto ThreeScrewHolderGrooveInnerDiamter = 36.322 * mm;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 

const auto ThreeScrewHolderLength = ThreeScrewHolderBottomLength + ThreeScrewHolderGrooveWidth + ThreeScrewHolderTopLength;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 
const auto ThreeScrewHolderOuterDiameter = 38.1 * mm;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 
const auto ThreeScrewHolderInnerDiameter = 0. * mm;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 

const auto ThreeScrewHolderAirOuterDiamter = 26.416 * mm;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 
const auto ThreeScrewHolderAirLength = 12.2682 * mm;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 

const auto SubtractionSolidBuffer1 = 0.8875 * mm;//From three screw holder technical draw: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026 

//------Activation Target -----

const auto ActivationTargetThicknessSmO = 0.190 * mm;
const auto ActivationTargetOuterDiameter = 12. * mm;
const auto ActivationTargetThicknessCeO = 0.550 * mm;

//-----Activation Target Container---

const auto ActivationTargetContainerOuterDiameter = 25.4 * mm; //From activation target container draw: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
//const auto ActivationTargetContainerInnerDiameter = 12. * mm;
const auto ActivationTargetHeight = 4.5 * mm;//From activation target container draw: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
const auto ActivationTargetContainerThickness = 1.5 * mm;//From activation target container draw: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
const auto ActivationTargetContainerLength = ActivationTargetHeight * mm + ActivationTargetThicknessSmO;
const auto ActivationTargetContainerLengthCeO = ActivationTargetHeight * mm + ActivationTargetThicknessCeO;
const auto ActivationTargetContainerAirWidth = ActivationTargetThicknessSmO;
const auto ActivationTargetContainerAirWidthCeO = ActivationTargetThicknessCeO;
const auto ActivationTargetContainerAirOuterDiameter = 25.4 * mm; //From activation target container draw: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
const auto ActivationTargetContainerAirInnerDiameter = 12. * mm; //From activation target container draw: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829

//-------Activation Target Au

const auto ActivationTargetAuOuterDiameter = 12.7 * mm;
const auto ActivationTargetAuThickness = 0.02 * mm; //not the real value just for building 

//-----Activation Target Container Lid----

const auto ActivationTargetLidOuterDiameter = 12. * mm;//From activation target container draw: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
const auto ActivationTargetLidInnerDiameter = 10. * mm;//From activation target container draw: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
const auto ActivationTargetLidInnerHeight = 2. *mm + ActivationTargetThicknessSmO;
const auto ActivationTargetLidInnerHeightCeO = 2. *mm + ActivationTargetThicknessCeO;
const auto ActivationTargetLidHeight = 3. * mm;//From activation target container draw: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
const auto ActivationTargetLidHeightCeO = 3. * mm; //From activation target container draw: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829

//------Green Holder------

const auto greenholderheight = 69.85 * mm;//From technical draw 
const auto greenholderlength = 254.0 * mm;//From technical draw 
const auto greenholderwidth = 2.921 * mm;//From technical draw 
const auto greenholderbottomlength = 88.9 * mm;//From technical draw 
const auto greenholderdistance = 82.55 * mm;//From technical draw

//------Offset------

const auto offset= 0.0508*mm;

//-----Slide----- 

auto const slidelength = 82.55 * mm;//From technical draw
auto const slideheight = 63.5 * mm;//From technical draw
auto const slidewidth = 1.58369 * mm;//From technical draw
auto const holediameter = 36.322 * mm;//From technical draw
auto const holewidth = 1.58369 * mm;//From technical draw

// ------Detector Holder-----

const auto greendetectorholderouterdiameterI = 101.6 * mm;//From technical draw
const auto greendetectorholderinnerdiameterI = 58.42 * mm;//From technical draw
const auto greendetectorholderwidthI = 1.778 * mm;//From technical draw

const auto greendetectorholderouterdiameterII = 101.6 * mm;//From technical draw
const auto greendetectorholderinnerdiameterII = 76.2 * mm; //From technical draw
const auto greendetectorholderwidthII = 17.653 * mm;//From technical draw

const auto greendetectorholderouterdiameterIII = 88.9 * mm;//From technical draw
const auto greendetectorholderinnerdiameterIII = 76.2 * mm;//From technical draw
const auto greendetectorholderwidthIII = 38.1 * mm;//From technical draw

//------Detector Planar----

const auto detectorFirstLayerOuterDiameter = 76.1 * mm;
const auto detectorFirstLayerInnerDiameter = 58.42 * mm;
const auto detectorFirstLayerThickness = 0.1 * mm;

const auto detectorWindowOuterDiameter = 58.42 * mm;
const auto detectorWindowInnerDiameter = 0. * mm;
const auto detectorWindowThickness = 0.5 * mm;

const auto detectorSecondLayerOuterDiameter = 76.1 * mm;
const auto detectorSecondLayerInnerDiameter = 75.6 * mm;
const auto detectorSecondLayerLength = 140. * mm;

const auto detectorEndCapOuterDiameter = 76.1 * mm;
const auto detectorEndCapInnerDiameter = 0. * mm;
const auto detectorEndCapThickness = 0.1 * mm;

const auto detectorCrystalOuterDiameter = 58.42 * mm;
const auto detectorCrystalInnerDiameter = 0. * mm;
const auto detectorCrystalLength = 21.* mm;

const auto detectorCrystalBerylliumWindowDistance = 5. * mm;

//-------Detector Coaxial-----

const auto detectorCoaxialFirstLayerOuterDiameter = 76.1 * mm;
const auto detectorCoaxialFirstLayerInnerDiameter = 58.42 * mm;
const auto detectorCoaxialFirstLayerThickness = 0.1 * mm;

const auto detectorCoaxialWindowOuterDiameter = 58.42 * mm;
const auto detectorCoaxialWindowInnerDiameter = 0. * mm;
const auto detectorCoaxialWindowThickness = 0.5 * mm;

const auto detectorCoaxialSecondLayerOuterDiameter = 76.1 * mm;
const auto detectorCoaxialSecondLayerInnerDiameter = 75.6 * mm;
const auto detectorCoaxialSecondLayerLength = 140. * mm;

const auto detectorCoaxialEndCapOuterDiameter = 76.1 * mm;
const auto detectorCoaxialEndCapInnerDiameter = 0. * mm;
const auto detectorCoaxialEndCapThickness = 0.1 * mm;

const auto detectorCoaxialCrystalOuterDiameter = 56. * mm;
const auto detectorCoaxialCrystalInnerDiameter = 0. * mm;
const auto detectorCoaxialCrystalLength = 53.5* mm;

const auto detectorCoaxialCrystalAirOuterDiameter = 9. * mm;
const auto detectorCoaxialCrystalAirLength = 22. * mm;


const auto detectorCrystalAluminiumWindowDistance = 5. * mm;

//-----Table-----

const auto tableheight =700. * mm;
const auto tablelength = 250. * mm;
const auto tablewidth = 2.5 * mm;

const auto tablefeetlength = 200. *mm;
const auto tablefeetheight = 300.0 * mm;
const auto tablefeetwidth = 30. *mm;

//------World-------------

const auto worldX = (activationRoomWidth + activationSideWallThickness) * 2;
const auto worldZ = (activationRoomLength + activationSideWallThickness) * 2;
const auto worldY = (activationRoomHeight) * 2;
 
auto *worldSolid = new G4Box("worldSolid", worldX / 2., worldY / 2., worldZ / 2.);
auto *worldLogical = new G4LogicalVolume(worldSolid, nist->FindOrBuildMaterial("G4_AIR"), "worldlogical");
G4VPhysicalVolume *worldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), worldLogical, "world", nullptr, false, 0);
auto worldVis = G4VisAttributes(red);
worldVis.SetForceWireframe(true);
worldLogical->SetVisAttributes(worldVis);

//-------Floor and Walls--------

auto *activationFloorSolid = new G4Box("activationFloorSolid", worldX / 2., activationRoomFloorThickness / 2., worldZ / 2.);
auto *activationFloorLogical = new G4LogicalVolume(activationFloorSolid, nist->FindOrBuildMaterial("G4_CONCRETE"), "activationFloorLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,-(activationRoomHeight + activationRoomFloorThickness / 2), 0), activationFloorLogical, "activationFloor", worldLogical, false, 0);
activationFloorLogical->SetVisAttributes(grey);


auto *activationSideWallSolid = new G4Box("activationSideWallSolid", worldX / 2., activationRoomHeight , activationSideWallThickness / 2.);
auto *activationSideWallLogical = new G4LogicalVolume(activationSideWallSolid, nist->FindOrBuildMaterial("G4_CONCRETE"), "activationSideWallLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0, 0, activationRoomLength + 75. * mm), activationSideWallLogical, "activationSideWall", worldLogical, false, 0);
auto activationSideWallVis = G4VisAttributes(grey);
activationSideWallVis.SetForceWireframe(true);
activationSideWallLogical->SetVisAttributes(activationSideWallVis);

auto *activationSideWallSolidleft = new G4Box("activationSideWallSolidleft", worldX / 2., activationRoomHeight , activationSideWallThickness / 2.);
auto *activationSideWallleftLogical = new G4LogicalVolume(activationSideWallSolidleft, nist->FindOrBuildMaterial("G4_CONCRETE"), "activationSideWallleftLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -activationRoomLength - 75. * mm), activationSideWallleftLogical, "activationSideWallleft", worldLogical, false, 0);
auto activationSideWallleftVis = G4VisAttributes(grey);
activationSideWallleftVis.SetForceWireframe(true);
activationSideWallleftLogical->SetVisAttributes(activationSideWallleftVis);

//-----Table-----

auto *tableSolid = new G4Box("tableSolid", tablelength / 2., tablewidth / 2., tableheight / 2.);
auto * tableLogical = new G4LogicalVolume(tableSolid, nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "table;Logical");
new G4PVPlacement(nullptr, G4ThreeVector(0,-greenholderheight /2. - greenholderwidth / 2. - 2 * tablewidth, 0), tableLogical, "table", worldLogical, false, 0);
tableLogical->SetVisAttributes(darkGrey);

auto *tablefeetSolid = new G4Box("tablefeetSolid", tablefeetlength / 2., tablefeetheight / 2., tablefeetwidth / 2.);
auto * tablefeetLogical = new G4LogicalVolume(tablefeetSolid, nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "tablefeetLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,-greenholderheight / 2. - greenholderwidth - 2 * tablewidth - tablelength / 2. -ThreeScrewHolderOuterDiameter / 2. -5. * mm, -tableheight / 2. + tablefeetwidth / 2. ), tablefeetLogical, "table", worldLogical, false, 0);
tablefeetLogical->SetVisAttributes(darkGrey);

auto *tablefeetrigthSolid = new G4Box("tablefeetrightSolid", tablefeetlength / 2., tablefeetheight / 2., tablefeetwidth / 2.);
auto * tablefeetrightLogical = new G4LogicalVolume(tablefeetSolid, nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "tablefeetrightLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,-greenholderheight / 2. - greenholderwidth - 2 * tablewidth - tablelength / 2. -ThreeScrewHolderOuterDiameter / 2. -5. * mm, +tableheight / 2. - tablefeetwidth / 2. ), tablefeetrightLogical, "table", worldLogical, false, 0);
tablefeetrightLogical->SetVisAttributes(darkGrey);

//------ Green Holder-----

auto *greenholderBoxLeftSolid = new G4Box("greenholderBoxLeftSolid", greenholderwidth / 2., greenholderheight / 2., greenholderlength / 2.);
auto *greenholderBoxLeftLogical = new G4LogicalVolume(greenholderBoxLeftSolid, nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "greenholderBoxLeftLogical");
new G4PVPlacement(nullptr, G4ThreeVector(greenholderdistance / 2. + greenholderwidth/2. ,0,-77.*mm), greenholderBoxLeftLogical, "greenholdetBocLeft", worldLogical, false, 0); // greenholderlength = 254.0 * mm and if the activation target is in the center og it -> 254.0 / 2 = 127. * mm. Since activation target is 50.* mm from the detector -> 127. *mm - 50.*mm = 77. * mm; So I shifting the green holder with -77.*mm (to the left)
auto greenholderBoxLeftVis = G4VisAttributes(green);
greenholderBoxLeftVis.SetForceWireframe(true);
greenholderBoxLeftLogical->SetVisAttributes(greenholderBoxLeftVis);

auto *greenholderBoxRightSolid = new G4Box("greenholderBoxRightSolid", greenholderwidth / 2., greenholderheight / 2., greenholderlength / 2.);
auto*greenholderBoxRightLogical = new G4LogicalVolume(greenholderBoxRightSolid, nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "greenholderBoxRightLogical");
new G4PVPlacement(nullptr, G4ThreeVector(-greenholderdistance / 2. - greenholderwidth/2., 0, -77.*mm), greenholderBoxRightLogical, "greenholderBoxRight", worldLogical, false, 0);
auto greenholderBoxRightVis = G4VisAttributes(green);
greenholderBoxRightVis.SetForceWireframe(true);
greenholderBoxRightLogical->SetVisAttributes(greenholderBoxRightVis);

auto *greenholderBoxBottomSolid = new G4Box("greenholderBoxBottomSolid", greenholderbottomlength / 2., greenholderwidth / 2., greenholderlength / 2.);
auto *greenholderBoxBottomLogical = new G4LogicalVolume(greenholderBoxBottomSolid, nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "greenholderBoxBottomLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,-tablewidth / 2. - greenholderheight / 2., -77.*mm), greenholderBoxBottomLogical, "greenholderBoxBottom", worldLogical, false, 0);
greenholderBoxBottomLogical->SetVisAttributes(green);

//------Green Detector Holder -----

auto *greendetectorholderISolid = new G4Tubs("greedetectorholderSolid", greendetectorholderinnerdiameterI / 2., greendetectorholderouterdiameterI / 2., greendetectorholderwidthI/ 2., 0, twopi);
auto *greendetectorholderILogical = new G4LogicalVolume(greendetectorholderISolid, nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "greendetectorholderILogical");
new G4PVPlacement(nullptr, G4ThreeVector(0, +tablewidth / 2. + greendetectorholderwidthI / 2., +50. * mm+ greendetectorholderwidthI / 2.), greendetectorholderILogical, "greendetectorholder", worldLogical, false, 0);
  greendetectorholderILogical->SetVisAttributes(green);

auto *greendetectorholderIISolid = new G4Tubs("greedetectorholderSolid", greendetectorholderinnerdiameterII / 2., greendetectorholderouterdiameterII / 2., greendetectorholderwidthII/ 2., 0, twopi);
auto *greendetectorholderIILogical = new G4LogicalVolume(greendetectorholderIISolid, nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "greendetectorholderIILogical");
new G4PVPlacement(nullptr, G4ThreeVector(0, +tablewidth / 2. + greendetectorholderwidthI / 2., +50. * mm + greendetectorholderwidthI + greendetectorholderwidthII / 2.),greendetectorholderIILogical, "greendetectorholder", worldLogical, false, 0);
greendetectorholderIILogical->SetVisAttributes(green);

auto *greendetectorholderIIISolid = new G4Tubs("greendetectorholderIII", greendetectorholderinnerdiameterIII / 2., greendetectorholderouterdiameterIII / 2., greendetectorholderwidthIII / 2., 0, twopi);
auto *greendetectorholderIIILogical = new G4LogicalVolume(greendetectorholderIIISolid, nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), "greendetectorholderIIILogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,+tablewidth / 2. + greendetectorholderwidthI / 2., +50. * mm + greendetectorholderwidthI + greendetectorholderwidthII + greendetectorholderwidthIII / 2.), greendetectorholderIIILogical, "greendetectorholder", worldLogical, false, 0);
greendetectorholderIIILogical->SetVisAttributes(green);

//------Slide----

auto *slideSolidBox = new G4Box("slideSolidBox", slidelength / 2., slideheight / 2., slidewidth / 2.);
auto *slideSolidHole = new G4Tubs("slideSolidHole", 0., holediameter / 2., holewidth / 2 + subtractionSolidBuffer, 0, twopi);
auto *slideSolid = new G4SubtractionSolid("slideSolid", slideSolidBox, slideSolidHole, nullptr, G4ThreeVector(0, 0, 0));

auto *slideLogical = new G4LogicalVolume(slideSolid, nist->FindOrBuildMaterial("G4_POLYETHYLENE"), "slideLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), slideLogical, "slide", worldLogical, false, 0);
slideLogical->SetVisAttributes(lightGrey); 

//----Three Screw Holder-----

const auto threescrewholderlidthickness = 1.6002 * mm;
const auto densityAcrylicMaterial = 1.18 * g/cm3;

auto *Acrylic = new G4Material("Acrylic", densityAcrylicMaterial, 3);
Acrylic->AddElement(nist->FindOrBuildElement("C"),5);
Acrylic->AddElement(nist->FindOrBuildElement("H"),8);
Acrylic->AddElement(nist->FindOrBuildElement("O"),2);

auto *ThreeScrewHolderSolid = new G4Tubs("ThreeScrewHolderSolid", ThreeScrewHolderInnerDiameter / 2., ThreeScrewHolderOuterDiameter / 2., ThreeScrewHolderLength / 2., 0, twopi);
auto *GrooveSolid = new G4Tubs("GrooveSolid", ThreeScrewHolderGrooveInnerDiamter / 2., ThreeScrewHolderGrooveOuterDiameter / 2. + SubtractionSolidBuffer1 / 2., ThreeScrewHolderGrooveWidth / 2., 0, twopi);

auto *ThreeScrewHolderGrooveSub = new G4SubtractionSolid("ThreeScrewHolderGrooveSub", ThreeScrewHolderSolid, GrooveSolid, nullptr, G4ThreeVector(0,0, +(ThreeScrewHolderLength / 2. - ThreeScrewHolderTopLength - ThreeScrewHolderGrooveWidth / 2.)));

auto *ThreeScrewHolderLogical = new G4LogicalVolume(ThreeScrewHolderGrooveSub,Acrylic, "ThreeScrewHolderLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0,-(ThreeScrewHolderLength / 2. - ThreeScrewHolderTopLength - ThreeScrewHolderGrooveWidth / 2.)), ThreeScrewHolderLogical, "ThreeScrewHolder", worldLogical, false, 0);
ThreeScrewHolderLogical->SetVisAttributes(grey);

auto *ThreeScrewHolderAirSolid = new G4Tubs("ThreeScrewHolderAirSolid", ThreeScrewHolderInnerDiameter / 2., ThreeScrewHolderAirOuterDiamter / 2., ThreeScrewHolderAirLength / 2., 0, twopi);
auto *ThreeScrewHolderAirLogical = new G4LogicalVolume(ThreeScrewHolderAirSolid, nist->FindOrBuildMaterial("G4_AIR"), "ThreeScrewHolderAirLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0, ThreeScrewHolderLength / 2. - ThreeScrewHolderAirLength / 2. - threescrewholderlidthickness),ThreeScrewHolderAirLogical, "threescrewholderAIR", ThreeScrewHolderLogical, false, 0);
ThreeScrewHolderAirLogical->SetVisAttributes(blue);

const auto activationTargetInnerDiameter = 12. * mm;
double targetnatSmO_Density ;
double natSmO_Container_Mass;
double targetnatCeO_Density;
double natCeO_Container_Mass;
double Au_Mass ;
double targetAu_Density;

std::cout << "DetectorConstruction: Requested TARGET is '" << TARGET << "'\n";
if (TARGET == "SmO_8") {
natSmO_Container_Mass = (2.05587-1.88185) * g;
} else if (TARGET == "SmO_9"){
 natSmO_Container_Mass = (2.07222-1.88257) * g;
} else if (TARGET == "SmO_10"){
 natSmO_Container_Mass = (2.11479-1.88234) * g;
} else if (TARGET =="SmO_11"){
 natSmO_Container_Mass = (2.07366-1.88180) * g;
} else if (TARGET == "SmO_12"){
 natSmO_Container_Mass = (2.03059-1.88428) * g;
} else if (TARGET == "SmO_13"){
 natSmO_Container_Mass = (2.09876-1.87385) * g;
} else if (TARGET == "SmO_14"){
 natSmO_Container_Mass = (2.05948-1.88342) * g;
} else if (TARGET == "CeO_1"){
 natCeO_Container_Mass = (2.35393-1.88349) * g;
} else if (TARGET == "CeO_2"){
 natCeO_Container_Mass = (2.33186-1.88597) * g;
} else if (TARGET == "CeO_3"){
natCeO_Container_Mass = (2.30733-1.88050) * g;
} else if (TARGET == "CeO_4"){
 natCeO_Container_Mass = (2.31009-1.87927) * g;
} else if (TARGET == "CeO_5"){
natCeO_Container_Mass = (2.33911-1.88289) * g;
} else if (TARGET == "CeO_6"){
 natCeO_Container_Mass = (2.35403-1.88406) * g;
} else if (TARGET == "CeO_7"){
natCeO_Container_Mass = (2.37485-1.88253) * g;
} else if (TARGET == "Au_1"){
 Au_Mass = (0.05169) * g;
} else if (TARGET == "Au_2"){
 Au_Mass = (0.04535) * g;
} else if (TARGET == "Au_3"){
 Au_Mass = (0.04798) * g;
 }else if (TARGET == "Au_4"){
 Au_Mass = (0.04818) * g;
} else if (TARGET == "Au_5"){
Au_Mass = (0.04950) * g;
} else if (TARGET == "Au_6"){
Au_Mass = (0.04987) * g;
} else if (TARGET == "Au_7"){
 Au_Mass = (0.05209) * g;
} else if (TARGET == "Au_8"){
 Au_Mass = (0.04526) * g;
} else if (TARGET == "Au_9"){
 Au_Mass = (0.04776) * g;
} else if (TARGET == "Au_10"){
 Au_Mass = (0.04534) * g;
} else if (TARGET == "Au_11"){
 Au_Mass = (0.04860) * g;
} else if (TARGET == "Au_12"){
 Au_Mass = (0.05168) * g;
} else if (TARGET == "MixSource3"){

const auto calibrationMixSourceDiameter3 = 25.4 * mm;
const auto calibrationMixSourceThickness3 = 3.734 * mm; //not sure, from excel Sean sent

auto *calibrationMixSource3Solid = new G4Tubs("calibrationMixSource3Solid", 0., calibrationMixSourceDiameter3 / 2., calibrationMixSourceThickness3 / 2., 0, twopi);
auto *calibrationMixSource3Logical = new G4LogicalVolume(calibrationMixSource3Solid, Acrylic, "calibrationMixSource");
new G4PVPlacement(nullptr,G4ThreeVector(0,0,0), calibrationMixSource3Logical,"calibrationMixSource", ThreeScrewHolderAirLogical, false, 0);
calibrationMixSource3Logical->SetVisAttributes(red);
} else if (TARGET == "MixSource4"){

const auto calibrationMixSourceDiameter4 = 25.4 * mm;
const auto calibrationMixSourceThickness4 = 3.734 * mm; //not sure, from excel Sean sent

auto *calibrationMixSource4Solid = new G4Tubs("calibrationMixSource4Solid", 0., calibrationMixSourceDiameter4 / 2., calibrationMixSourceThickness4 / 2., 0, twopi);
auto *calibrationMixSource4Logical = new G4LogicalVolume(calibrationMixSource4Solid, Acrylic, "calibrationMixSource");
new G4PVPlacement(nullptr,G4ThreeVector(0,0,0), calibrationMixSource4Logical,"calibrationMixSource", ThreeScrewHolderAirLogical, false, 0);
calibrationMixSource4Logical->SetVisAttributes(red);
} else if (TARGET == "none"){}

if (TARGET == "SmO_8" || TARGET == "SmO_9" || TARGET == "SmO_10" || TARGET == "SmO_11" || TARGET == "SmO_12" || TARGET == "SmO_13" || TARGET == "SmO_14" ){

targetnatSmO_Density = natSmO_Container_Mass / (pi / 4. * activationTargetInnerDiameter * activationTargetInnerDiameter * ActivationTargetThicknessSmO);

auto *activation154Sm_nat_container_Material = new G4Material("natSmO_Container_Material", targetnatSmO_Density, 2);
activation154Sm_nat_container_Material->AddElement(nist->FindOrBuildElement("Sm"), 2); // 2 Sm atoms in Sm(2)O(3)
activation154Sm_nat_container_Material->AddElement(nist->FindOrBuildElement("O"), 3); // 3 O atoms in Sm(2)O(3)

auto *ActivationTargetContainerSolid = new G4Tubs("ActivationTargteHolderSolid", 0., ActivationTargetContainerOuterDiameter / 2., ActivationTargetContainerLength / 2., 0, twopi);
auto *ActivationTargetContainerLogical = new G4LogicalVolume(ActivationTargetContainerSolid, nist->FindOrBuildMaterial("G4_POLYETHYLENE"), "ActivationTargetContainerLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0, -ThreeScrewHolderAirLength / 2. + ActivationTargetContainerLength /2. +ThreeScrewHolderGrooveWidth / 2. +ThreeScrewHolderBottomAirWidth - ActivationTargetContainerThickness - offset), ActivationTargetContainerLogical,"activationtargetcontainer", ThreeScrewHolderAirLogical, false, 0);
ActivationTargetContainerLogical->SetVisAttributes(green); 

auto *ActivationTargetSolid = new G4Tubs("ActivationTargetSolid", 0.,ActivationTargetOuterDiameter / 2.,ActivationTargetThicknessSmO / 2., 0, twopi);
auto *ActivationTargetSolidLogical = new G4LogicalVolume(ActivationTargetSolid, activation154Sm_nat_container_Material,"ActivationTargetLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0, -ActivationTargetContainerLength/2.+ ActivationTargetThicknessSmO/2.+ ActivationTargetContainerThickness), ActivationTargetSolidLogical, "activationtarget", ActivationTargetContainerLogical, false, 0);
ActivationTargetSolidLogical->SetVisAttributes(orange);

auto *ActivationTargetContainerAirSolid = new G4Tubs("ActivationTargetContainerAirSolid", ActivationTargetContainerAirInnerDiameter / 2., ActivationTargetContainerAirOuterDiameter / 2., ActivationTargetContainerAirWidth / 2., 0, twopi);
auto *ActivationTargetContainerAirLogical = new G4LogicalVolume(ActivationTargetContainerAirSolid, nist->FindOrBuildMaterial("G4_AIR"),"ActivationTargetContainerAirSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,0, ActivationTargetContainerLength / 2. - ActivationTargetContainerAirWidth / 2.), ActivationTargetContainerAirLogical, "ActivationTargetContainerAir", ActivationTargetContainerLogical, false, 0);
ActivationTargetContainerAirLogical->SetVisAttributes(blue);

auto *ActivationTargetLidSolid = new G4Tubs("ActivationTargetLidSolid", 0., ActivationTargetLidOuterDiameter / 2., ActivationTargetLidHeight / 2., 0, twopi);
auto *ActivationTargetLidLogical = new G4LogicalVolume(ActivationTargetLidSolid, nist->FindOrBuildMaterial("G4_POLYETHYLENE"), "activationtargetlidLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0,+ ActivationTargetContainerLength/2. - ActivationTargetLidHeight /2 ), ActivationTargetLidLogical, "activationtargtelid", ActivationTargetContainerLogical, false, 0);
ActivationTargetLidLogical->SetVisAttributes(red);

auto *ActivationTargetLidAirSolid = new G4Tubs("ActivationTargetLidAirSolid", 0., ActivationTargetLidInnerDiameter / 2., ActivationTargetLidInnerHeight / 2., 0, twopi);
auto *ActivationTargetLidAirLogical = new G4LogicalVolume(ActivationTargetLidAirSolid, nist->FindOrBuildMaterial("G4_AIR"), "activationtargetLidAir");
new G4PVPlacement(nullptr,G4ThreeVector(0,0,+ActivationTargetLidHeight / 2. - ActivationTargetLidInnerHeight / 2. ), ActivationTargetLidAirLogical, "activationtargetLIDAIR", ActivationTargetLidLogical, false, 0);
ActivationTargetLidAirLogical->SetVisAttributes(blue);

}
else if ( TARGET == "CeO_1" || TARGET == "CeO_2" || TARGET == "CeO_3" || TARGET == "CeO_4" || TARGET == "CeO_5" || TARGET == "CeO_6"|| TARGET == "CeO_7"){
  
targetnatCeO_Density = natCeO_Container_Mass / (pi / 4. * activationTargetInnerDiameter * activationTargetInnerDiameter * ActivationTargetThicknessCeO);

auto *activationCeO_nat_container_Material = new G4Material("natSmO_Container8_Material", targetnatCeO_Density, 2);
activationCeO_nat_container_Material->AddElement(nist->FindOrBuildElement("Ce"), 1); // 1 Ce atoms in Ce(1)O(2)
activationCeO_nat_container_Material->AddElement(nist->FindOrBuildElement("O"), 2); // 32 O atoms in Ce(1)O(2)

auto *ActivationTargetContainerSolid = new G4Tubs("ActivationTargteHolderSolid", 0., ActivationTargetContainerOuterDiameter / 2., ActivationTargetContainerLengthCeO / 2., 0, twopi);
auto *ActivationTargetContainerLogical = new G4LogicalVolume(ActivationTargetContainerSolid, nist->FindOrBuildMaterial("G4_POLYETHYLENE"), "ActivationTargetContainerLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0, -ThreeScrewHolderAirLength / 2. + ActivationTargetContainerLength /2. +ThreeScrewHolderGrooveWidth / 2. +ThreeScrewHolderBottomAirWidth - ActivationTargetContainerThickness - offset), ActivationTargetContainerLogical,  "activationtargetcontainer", ThreeScrewHolderAirLogical, false, 0);
ActivationTargetContainerLogical->SetVisAttributes(green); 

auto *ActivationTargetSolid = new G4Tubs("ActivationTargetSolid", 0.,ActivationTargetOuterDiameter / 2.,  ActivationTargetThicknessCeO / 2., 0, twopi);
auto *ActivationTargetSolidLogical = new G4LogicalVolume(ActivationTargetSolid, activationCeO_nat_container_Material,"ActivationTargetLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0, -ActivationTargetContainerLengthCeO /2. + ActivationTargetThicknessCeO/2. + ActivationTargetContainerThickness), ActivationTargetSolidLogical, "activationtarget", ActivationTargetContainerLogical, false, 0);
ActivationTargetSolidLogical->SetVisAttributes(orange);

auto *ActivationTargetContainerAirSolid = new G4Tubs("ActivationTargetContainerAirSolid", ActivationTargetContainerAirInnerDiameter / 2., ActivationTargetContainerAirOuterDiameter / 2., ActivationTargetContainerAirWidthCeO / 2., 0, twopi);
auto *ActivationTargetContainerAirLogical = new G4LogicalVolume(ActivationTargetContainerAirSolid, nist->FindOrBuildMaterial("G4_AIR"),"ActivationTargetContainerAirSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,0, ActivationTargetContainerLengthCeO / 2. - ActivationTargetContainerAirWidthCeO / 2.), ActivationTargetContainerAirLogical, "ActivationTargetContainerAir", ActivationTargetContainerLogical, false, 0);
ActivationTargetContainerAirLogical->SetVisAttributes(yellow);

auto *ActivationTargetLidSolid = new G4Tubs("ActivationTargetLidSolid", 0., ActivationTargetLidOuterDiameter / 2., ActivationTargetLidHeightCeO / 2., 0, twopi);
auto *ActivationTargetLidLogical = new G4LogicalVolume(ActivationTargetLidSolid, nist->FindOrBuildMaterial("G4_POLYETHYLENE"), "activationtargetlidLogical"); 
new G4PVPlacement(nullptr, G4ThreeVector(0,0,+ ActivationTargetContainerLengthCeO/2. - ActivationTargetLidHeightCeO /2), ActivationTargetLidLogical, "activationtargtelid", ActivationTargetContainerLogical, false, 0);
ActivationTargetLidLogical->SetVisAttributes(red);

auto *ActivationTargetLidAirSolid = new G4Tubs("ActivationTargetLidAirSolid", 0., ActivationTargetLidInnerDiameter / 2., ActivationTargetLidInnerHeightCeO / 2., 0, twopi);
auto *ActivationTargetLidAirLogical = new G4LogicalVolume(ActivationTargetLidAirSolid, nist->FindOrBuildMaterial("G4_AIR"), "activationtargetLidAir");
new G4PVPlacement(nullptr,G4ThreeVector(0,0,+ActivationTargetLidHeightCeO / 2. - ActivationTargetLidInnerHeightCeO / 2.), ActivationTargetLidAirLogical, "activationtargetLIDAIR", ActivationTargetLidLogical, false, 0);
ActivationTargetLidAirLogical->SetVisAttributes(blue);
}
else if (TARGET == "AU_1" || TARGET == "Au_2" || TARGET == "Au_3" || TARGET == "Au_4" || TARGET == "Au_5" || TARGET == "Au_6" || TARGET == "Au_6" || TARGET == "Au_7" || TARGET == "Au_8" || TARGET == "Au_9" || TARGET == "Au_10" || TARGET == "Au_11" || TARGET == "Au_12"){

auto *ActivationTargetAuSolid = new G4Tubs("ActivationTargetAuSolid", 0., ActivationTargetAuOuterDiameter / 2., ActivationTargetAuThickness / 2., 0, twopi);
auto *ActivationTargetAuLogical = new G4LogicalVolume(ActivationTargetAuSolid, nist->FindOrBuildMaterial("G4_Au"), "ActivationTargetAuLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0, -ThreeScrewHolderAirLength / 2. +ThreeScrewHolderGrooveWidth / 2. +ThreeScrewHolderBottomAirWidth + offset), ActivationTargetAuLogical, "ActivationTargetAu", ThreeScrewHolderAirLogical, false,0);
ActivationTargetAuLogical->SetVisAttributes(orange);
}
if (TARGET == "SmO_8" || TARGET == "SmO_9" || TARGET == "SmO_10" || TARGET == "SmO_11" || TARGET == "SmO_12" || TARGET == "SmO_13" || TARGET == "SmO_14" || TARGET == "CeO_1" || TARGET == "CeO_2" || TARGET == "CeO_3" || TARGET == "CeO_4" || TARGET == "CeO_5" || TARGET == "CeO_6"|| TARGET == "CeO_7" || TARGET == "MixSource3"){

auto *detectorFirstLayerSolid = new G4Tubs("detectorFirstLayerSolid", detectorFirstLayerInnerDiameter / 2., detectorFirstLayerOuterDiameter / 2., detectorFirstLayerThickness / 2., 0, twopi);
auto *detectorFirstLayerLogical = new G4LogicalVolume(detectorFirstLayerSolid, nist->FindOrBuildMaterial("G4_Al"), "detectorFirstLayerSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,+tablewidth / 2. + greendetectorholderwidthI / 2., +50. * mm + greendetectorholderwidthI +detectorFirstLayerThickness / 2. ), detectorFirstLayerLogical, "detectorFirstLayer", worldLogical, false, 0);
detectorFirstLayerLogical->SetVisAttributes(grey);

auto *detectorWindowSolid = new G4Tubs("detectorWindowSolid", detectorWindowInnerDiameter / 2., detectorWindowOuterDiameter / 2., detectorWindowThickness / 2., 0, twopi);
auto *detectorWindowLogical = new G4LogicalVolume(detectorWindowSolid, nist->FindOrBuildMaterial("G4_Be"), "detectorWindowSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), detectorWindowLogical, "detectorFirstLayer", detectorFirstLayerLogical, false, 0);
detectorWindowLogical->SetVisAttributes(darkGrey);

auto *detectorSecondLayerSolid = new G4Tubs("detectorSecondLayerSolid", detectorSecondLayerInnerDiameter / 2., detectorSecondLayerOuterDiameter / 2., detectorSecondLayerLength / 2., 0, twopi);
auto *detectorSecondLayerLogical = new G4LogicalVolume(detectorSecondLayerSolid, nist->FindOrBuildMaterial("G4_Al"), "detectorSecondLayerSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,+tablewidth / 2. + greendetectorholderwidthI / 2., +50. * mm + greendetectorholderwidthI + detectorFirstLayerThickness + detectorSecondLayerLength / 2.), detectorSecondLayerLogical, "detectorSecondLayer", worldLogical, false, 0);
detectorSecondLayerLogical->SetVisAttributes(grey);

auto *detectorSecondLayerVacuumSolid = new G4Tubs("detectorSecondLayerVacuumSolid", 0., detectorSecondLayerInnerDiameter / 2., detectorSecondLayerLength / 2., 0, twopi);
auto *detectorSecondLayerVacuumLogical = new G4LogicalVolume(detectorSecondLayerVacuumSolid, nist->FindOrBuildMaterial("G4_Galactic"), "detectorSecondLayerVacuumSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,0.0), detectorSecondLayerVacuumLogical, "detectorSecondLayer", detectorSecondLayerLogical, false, 0);
detectorSecondLayerVacuumLogical->SetVisAttributes(invisible);

auto *detectorEndCapSolid = new G4Tubs("detectorEndCapSolid", detectorEndCapInnerDiameter / 2., detectorEndCapOuterDiameter / 2., detectorEndCapThickness / 2., 0, twopi);
auto *detectorEndCapLogical = new G4LogicalVolume(detectorEndCapSolid, nist->FindOrBuildMaterial("G4_Al"), "detectorEndCapSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,+tablewidth / 2. + greendetectorholderwidthI / 2., +50. * mm + greendetectorholderwidthI+ detectorFirstLayerThickness+ detectorSecondLayerLength + detectorEndCapThickness / 2.), detectorEndCapLogical, "detectorEndCap", worldLogical, false, 0);
detectorEndCapLogical->SetVisAttributes(grey);

auto *detectorCrystalSolid = new G4Tubs("detectorCrystalSolid", detectorCrystalInnerDiameter / 2., detectorCrystalOuterDiameter / 2., detectorCrystalLength / 2., 0, twopi);
auto *detectorCrystalLogical = new G4LogicalVolume(detectorCrystalSolid, nist->FindOrBuildMaterial("G4_Ge"),"detectorCrystalLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0,-detectorSecondLayerLength / 2. + detectorCrystalLength / 2. + detectorCrystalBerylliumWindowDistance / 2.), detectorCrystalLogical, "detectorCrystal", detectorSecondLayerVacuumLogical, false, 0);
detectorCrystalLogical->SetVisAttributes(blue);
}
else if (TARGET == "AU_1" || TARGET == "Au_2" || TARGET == "Au_3" || TARGET == "Au_4" || TARGET == "Au_5" || TARGET == "Au_6" || TARGET == "Au_6" || TARGET == "Au_7" || TARGET == "Au_8" || TARGET == "Au_9" || TARGET == "Au_10" || TARGET == "Au_11" || TARGET == "Au_12"){

auto *detectorCoaxialFirstLayerSolid = new G4Tubs("detectorCoaxialFirstLayerSolid", detectorCoaxialFirstLayerInnerDiameter / 2., detectorCoaxialFirstLayerOuterDiameter / 2., detectorCoaxialFirstLayerThickness / 2., 0, twopi);
auto *detectorCoaxialFirstLayerLogical = new G4LogicalVolume(detectorCoaxialFirstLayerSolid, nist->FindOrBuildMaterial("G4_Al"), "detectorCoaxialFirstLayerSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,+tablewidth / 2. + greendetectorholderwidthI / 2., +50. * mm + greendetectorholderwidthI +detectorCoaxialFirstLayerThickness / 2. ), detectorCoaxialFirstLayerLogical, "detectorFirstLayer", worldLogical, false, 0);
detectorCoaxialFirstLayerLogical->SetVisAttributes(grey);

auto *detectorCoaxialWindowSolid = new G4Tubs("detectorCoaxialWindowSolid", detectorCoaxialWindowInnerDiameter / 2., detectorCoaxialWindowOuterDiameter / 2., detectorCoaxialWindowThickness / 2., 0, twopi);
auto *detectorCoaxialWindowLogical = new G4LogicalVolume(detectorCoaxialWindowSolid, nist->FindOrBuildMaterial("G4_Be"), "detectorCoaxialWindowSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), detectorCoaxialWindowLogical, "detectorFirstLayer", detectorCoaxialFirstLayerLogical, false, 0);
detectorCoaxialWindowLogical->SetVisAttributes(darkGrey);

auto *detectorCoaxialSecondLayerSolid = new G4Tubs("detectorCoaxialSecondLayerSolid", detectorCoaxialSecondLayerInnerDiameter / 2., detectorCoaxialSecondLayerOuterDiameter / 2., detectorCoaxialSecondLayerLength / 2., 0, twopi);
auto *detectorCoaxialSecondLayerLogical = new G4LogicalVolume(detectorCoaxialSecondLayerSolid, nist->FindOrBuildMaterial("G4_Al"), "detectorCoaxialSecondLayerSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,+tablewidth / 2. + greendetectorholderwidthI / 2., +50. * mm + greendetectorholderwidthI + detectorCoaxialFirstLayerThickness + detectorCoaxialSecondLayerLength / 2.), detectorCoaxialSecondLayerLogical, "detectorSecondLayer", worldLogical, false, 0);
detectorCoaxialSecondLayerLogical->SetVisAttributes(grey);

auto *detectorCoaxialSecondLayerVacuumSolid = new G4Tubs("detectorCoaxialSecondLayerVacuumSolid", 0., detectorCoaxialSecondLayerInnerDiameter / 2., detectorCoaxialSecondLayerLength / 2., 0, twopi);
auto *detectorCoaxialSecondLayerVacuumLogical = new G4LogicalVolume(detectorCoaxialSecondLayerVacuumSolid, nist->FindOrBuildMaterial("G4_Galactic"), "detectorCoaxialSecondLayerVacuumSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), detectorCoaxialSecondLayerVacuumLogical, "detectorSecondLayer", detectorCoaxialSecondLayerLogical, false, 0);
detectorCoaxialSecondLayerVacuumLogical->SetVisAttributes(invisible);

auto *detectorCoaxialEndCapSolid = new G4Tubs("detectorCoaxialEndCapSolid", detectorCoaxialEndCapInnerDiameter / 2., detectorCoaxialEndCapOuterDiameter / 2., detectorCoaxialEndCapThickness / 2., 0, twopi);
auto *detectorCoaxialEndCapLogical = new G4LogicalVolume(detectorCoaxialEndCapSolid, nist->FindOrBuildMaterial("G4_Al"), "detectorCoaxialEndCapSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,+tablewidth / 2. + greendetectorholderwidthI / 2., +50. * mm + greendetectorholderwidthI+ detectorCoaxialFirstLayerThickness+ detectorCoaxialSecondLayerLength + detectorCoaxialEndCapThickness / 2.), detectorCoaxialEndCapLogical, "detectorEndCap", worldLogical, false, 0);
detectorCoaxialEndCapLogical->SetVisAttributes(grey);

auto *detectorCoaxialCrystalSolid = new G4Tubs("detectorCoaxialCrystalSolid", detectorCoaxialCrystalInnerDiameter / 2., detectorCoaxialCrystalOuterDiameter / 2., detectorCoaxialCrystalLength / 2., 0, twopi);
auto *detectorCoaxialCrystalLogical = new G4LogicalVolume(detectorCoaxialCrystalSolid, nist->FindOrBuildMaterial("G4_Ge"),"detectorCoaxialCrystalLogical");
new G4PVPlacement(nullptr, G4ThreeVector(0,0,-detectorCoaxialSecondLayerLength / 2. + detectorCoaxialCrystalLength / 2. + detectorCrystalAluminiumWindowDistance / 2.), detectorCoaxialCrystalLogical, "detectorCrystal", detectorCoaxialSecondLayerVacuumLogical, false, 0);
detectorCoaxialCrystalLogical->SetVisAttributes(blue);

auto *detectorCoaxialAirLayerVacuumSolid = new G4Tubs("detectorCoaxialAirLayerVacuumSolid", 0., detectorCoaxialCrystalAirOuterDiameter / 2., detectorCoaxialCrystalAirLength / 2., 0, twopi);
auto *detectorCoaxialAirLayerVacuumLogical = new G4LogicalVolume(detectorCoaxialAirLayerVacuumSolid, nist->FindOrBuildMaterial("G4_Galactic"), "detectorCoaxialAirLayerVacuumSolid");
new G4PVPlacement(nullptr, G4ThreeVector(0,0., detectorCoaxialCrystalLength / 2. - detectorCoaxialCrystalAirLength / 2.), detectorCoaxialAirLayerVacuumLogical, "detectorSecondLayer", detectorCoaxialCrystalLogical, false, 0);
detectorCoaxialAirLayerVacuumLogical->SetVisAttributes(cyan);

} 

return worldPhysical;

}
void DetectorConstruction::ConstructSDandField() {
}