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

  // =============== Measurements ===============

  // --------------- General ---------------

  string targetType;
  if (TARGET == "Sm_No8" || TARGET == "Sm_No9" || TARGET == "Sm_No10" || TARGET == "Sm_No11" || TARGET == "Sm_No12" || TARGET == "Sm_No13" || TARGET == "Sm_No14") {
    targetType = "Sm";
  } else if (TARGET == "Ce_No1" || TARGET == "Ce_No2" || TARGET == "Ce_No3" || TARGET == "Ce_No4" || TARGET == "Ce_No5" || TARGET == "Ce_No6" || TARGET == "Ce_No7") {
    targetType = "Ce";
  } else if (TARGET == "Au_No1" || TARGET == "Au_No2" || TARGET == "Au_No3" || TARGET == "Au_No4" || TARGET == "Au_No5" || TARGET == "Au_No6" || TARGET == "Au_No6" || TARGET == "Au_No7" || TARGET == "Au_No8" || TARGET == "Au_No9" || TARGET == "Au_No10" || TARGET == "Au_No11" || TARGET == "Au_No12") {
    targetType = "Au";
  } else if (TARGET == "MixedSourceAtDet3" || TARGET == "MixedSourceAtDet4") {
    targetType = "MixedSource";
  } else if (TARGET == "PointSourceAtDet3" || TARGET == "PointSourceAtDet4") {
    targetType = "PointSource";
  } else {
    G4cerr << "ERROR: Unknown TARGET '" << TARGET << "' was requested for the DetectorConstruction! Aborting..." << G4endl;
    throw std::exception();
  }
  std::cout << "DetectorConstruction: Requested TARGET is '" << TARGET << "' of type '" << targetType << "'\n";

  const auto useDet4InsteadOfDet3 = (targetType == "Au" || TARGET == "MixedSourceAtDet4" || TARGET == "PointSourceAtDet4");

  const auto subtractionSolidBuffer = 10. * mm; // G4SubtractionSolid doesn't always fully cut solids when the cutting solid is equal in dimension to the to be cut solid, hence, add some buffer

  const auto acrylic = nist->FindOrBuildMaterial("G4_PLEXIGLASS"); // Acrylic and Plexiglass are the same thing
  const auto vacuum = nist->FindOrBuildMaterial("G4_Galactic");

  // --------------- World ---------------

  const auto worldXLength = 200. * mm;
  const auto worldYWidth = 200. * mm;
  const auto worldZHeight = 500. * mm;

  // --------------- Target Holder ---------------

  const auto targetHolderDet3DetectorEnclosingFaceHoleDiameter = 2.5 * inch; // From technical drawing "TargetHolder_det3"
  const auto targetHolderDet3DetectorEnclosingInnerDiameter = 3.24 * inch; // From technical drawing "TargetHolder_det3"
  const auto targetHolderDet3DetectorEnclosingOuterDiameter = 4 * inch; // From technical drawing "TargetHolder_det3"
  const auto targetHolderDet3DetectorEnclosingLength = 28.7 * mm; // Not explicitly given in the drawing "TargetHolder_det3", hence manually measured off it
  const auto targetHolderDet3DetectorEnclosingFrontFaceThickness = 0.08 * inch; // From technical drawing "TargetHolder_det3"
  const auto targetHolderDet3SlideHolderLength = 6. * inch; // From technical drawing "TargetHolder_det3"
  const auto targetHolderDet3SlideHolderHeight = 2.63 * inch; // From technical drawing "TargetHolder_det3"
  const auto targetHolderDet3SlideHolderWidth = 3.5 * inch; // From technical drawing "TargetHolder_det3"
  const auto targetHolderDet3SlideHolderWallThickness = 0.15 * inch; // From technical drawing "TargetHolder_det3"
  const auto targetHolderDet3Material = new G4Material("plastic3DPrintWithOnly20PercentInfillMaterial", 0.2 * 0.94 * g / cm3, nist->FindOrBuildMaterial("G4_POLYETHYLENE")); // Known: 3D printed with 20% infill, assume its polyethylene

  const auto targetHolderDet4DetectorEnclosingFrontFaceHoleDiameter = 2.3 * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4DetectorEnclosingInnerDiameter = 3. * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4DetectorEnclosingBackOuterDiameter = 3.5 * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4DetectorEnclosingFrontOuterDiameter = 4. * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4DetectorEnclosingBackLength = 1.5 * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4DetectorEnclosingFrontLength = 0.765 * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4DetectorEnclosingFrontFaceThickness = 0.07 * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4SlideHolderLength = 10. * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4SlideHolderHeight = 2.75 * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4SlideHolderWidth = 3.4 * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4SlideHolderWallThickness = 0.115 * inch; // From technical drawing "Gooden_holder"
  const auto targetHolderDet4Material = acrylic; // From S. Finch, private communication

  const auto targetHolderDetectorEnclosingFrontFaceHoleDiameter = useDet4InsteadOfDet3 ? targetHolderDet4DetectorEnclosingFrontFaceHoleDiameter : targetHolderDet3DetectorEnclosingFaceHoleDiameter;
  const auto targetHolderDetectorEnclosingInnerDiameter = useDet4InsteadOfDet3 ? targetHolderDet4DetectorEnclosingInnerDiameter : targetHolderDet3DetectorEnclosingInnerDiameter;
  const auto targetHolderDetectorEnclosingFrontOuterDiameter = useDet4InsteadOfDet3 ? targetHolderDet4DetectorEnclosingFrontOuterDiameter : targetHolderDet3DetectorEnclosingOuterDiameter;
  const auto targetHolderDetectorEnclosingFrontLength = useDet4InsteadOfDet3 ? targetHolderDet4DetectorEnclosingFrontLength : targetHolderDet3DetectorEnclosingLength;
  const auto targetHolderDetectorEnclosingFrontFaceThickness = useDet4InsteadOfDet3 ? targetHolderDet4DetectorEnclosingFrontFaceThickness : targetHolderDet3DetectorEnclosingFrontFaceThickness;
  const auto targetHolderSlideHolderLength = useDet4InsteadOfDet3 ? targetHolderDet4SlideHolderLength : targetHolderDet3SlideHolderLength;
  const auto targetHolderSlideHolderHeight = useDet4InsteadOfDet3 ? targetHolderDet4SlideHolderHeight : targetHolderDet3SlideHolderHeight;
  const auto targetHolderSlideHolderWidth = useDet4InsteadOfDet3 ? targetHolderDet4SlideHolderWidth : targetHolderDet3SlideHolderWidth;
  const auto targetHolderSlideHolderWallThickness = useDet4InsteadOfDet3 ? targetHolderDet4SlideHolderWallThickness : targetHolderDet3SlideHolderWallThickness;
  const auto targetHolderMaterial = useDet4InsteadOfDet3 ? targetHolderDet4Material : targetHolderDet3Material;

  // --------------- Slide ---------------

  auto const slideThickness = 0.0625 * inch; // From technical drawing "Slide"
  auto const slideHeight = std::min(2. * 1.25 * inch, useDet4InsteadOfDet3 ? (targetHolderDet4SlideHolderHeight - 2. * targetHolderDet4SlideHolderWallThickness) : (targetHolderDet3SlideHolderHeight - 2. * targetHolderDet3SlideHolderWallThickness)); // From technical drawing "Slide", but as grooves in the SlideHolder are not implemented here, just make the slideHeight as small as the "inner" height of the SlideHolder
  auto const slideWidth = std::min(3.25 * inch, useDet4InsteadOfDet3 ? (targetHolderDet4SlideHolderWidth - 2. * targetHolderDet4SlideHolderWallThickness) : (targetHolderDet3SlideHolderWidth - 2. * targetHolderDet3SlideHolderWallThickness)); // From technical drawing "Slide", but as grooves in the SlideHolder are not implemented here, just make the slideWidth as small as the inner length of the SlideHolder
  auto const slideHoleDiameter = 2. * 0.715 * inch; // From technical drawing "Slide"

  // --------------- Three-Screw Holder ---------------

  const auto threeScrewHolderInnerLength = 0.483 * inch; // From three screw holder technical drawing: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026
  const auto threeScrewLidThickness = 0.063 * inch; // From three screw holder technical drawing: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026
  const auto threeScrewHolderLength = 0.538 * inch + threeScrewLidThickness; // From three screw holder technical drawing: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026
  const auto threeScrewHolderBottomThickness = threeScrewHolderLength - threeScrewLidThickness - threeScrewHolderInnerLength; // From three screw holder technical drawing: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026
  const auto threeScrewHolderInnerDiameter = 1.040 * inch; // From three screw holder technical drawing: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026
  const auto threeScrewHolderOuterDiameter = 1.5 * inch; // From three screw holder technical drawing: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026
  const auto threeScrewHolderGrooveInnerDiameter = 1.43 * inch; // From three screw holder technical drawing: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026
  const auto threeScrewHolderGrooveWidth = 0.063 * inch; // From three screw holder technical drawing: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026
  const auto threeScrewHolderDistanceGrooveCenterToBottom = 0.270 * inch + threeScrewHolderGrooveWidth / 2.; // From three screw holder technical drawing: https://elog.ikp.physik.tu-darmstadt.de//clovershare/2026

  // --------------- Activation Target Container ---------------
  // Was used to encapsulate SmO and CeO powder, not used for Au foil disks

  const auto activationTargetContainerOuterDiameter = 25.4 * mm; // From activation target container technical drawing: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  const auto activationTargetContainerInnerDiameter = 12. * mm; // "In contrast to the attached drawings, each "Activation Target Holder" has an inner diameter of 12 mm as well instead of 14 mm due to a miscommunication" From: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  const auto activationTargetContainerLength = 4.5 * mm; // From activation target container technical drawing: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  const auto activationTargetContainerBottomThickness = 1.5 * mm; // From activation target container technical drawing: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  const auto activationTargetContainerLidLength = 3. * mm; // From activation target container technical drawing: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  const auto activationTargetContainerLidOuterDiameter = activationTargetContainerInnerDiameter; // From activation target container technical drawing: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  const auto activationTargetContainerLidRecessInnerDiameter = 10. * mm; // Assumed due to "In contrast to the attached drawings, each "Activation Target Holder" has an inner diameter of 12 mm as well instead of 14 mm due to a miscommunication" From: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  const auto activationTargetContainerLidRecessDepth = 2. * mm; // From activation target container technical drawing: https://elog.ikp.physik.tu-darmstadt.de/clovershare/829

  // --------------- Activation Targets ---------------

  // Geometric beam widening from 8mm collimator to target position:
  const auto collimatorRoomCollimatorAperture = 8. * mm;
  // From H.R. Weller et al. Prog. Part. Nucl. Phys. 62, 257 (2009): "collimator is located about 60 m away from the collision point."
  // From U. Friman-Gayer and S. Finch private communication: "It is 53 m from the beam collision point to collimator"
  const auto collisionPointToCollimator = 53. * m;
  const auto activationTargetHolderToTargetPos = 54.0 * inch; // From ELOG https://elog.ikp.physik.tu-darmstadt.de/clovershare/2025

  const auto utrFirstLeadWallToTargetPos = 60. * inch; // From nutr: Estimated
  const auto utrFirstLeadWallLength = 8. * inch; // From nutr
  const auto utrUpstreamWallToTargetPos = utrFirstLeadWallToTargetPos + utrFirstLeadWallLength;
  const auto utrUpstreamWallThickness = 162. * mm; // Measured off utr drawings
  const auto collimatorRoomSecondLeadWallToTargetPos = utrUpstreamWallToTargetPos + utrUpstreamWallThickness;
  const auto collimatorRoomSecondLeadWallLength = 16. * inch; // From nutr
  const auto collimatorRoomFirstLeadWallToTargetPos = collimatorRoomSecondLeadWallToTargetPos + collimatorRoomSecondLeadWallLength + 320. * mm; // From nutr
  const auto collimatorRoomFirstLeadWallLength = 16. * inch; // From nutr
  const auto collimatorRoomCollimatorToTargetPos = collimatorRoomFirstLeadWallToTargetPos + collimatorRoomFirstLeadWallLength + 200. * mm; // From nutr

  const auto collisionPointToActivationTargetPos = collisionPointToCollimator + collimatorRoomCollimatorToTargetPos + activationTargetHolderToTargetPos;
  const auto beamDiameterAtActivationTargetPos = collimatorRoomCollimatorAperture * collisionPointToActivationTargetPos / collisionPointToCollimator;

  const auto activationTargetAuDiameter = 0.5 * inch; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  const auto activationTargetAuThickness = 0.02 * mm; // According to the seller Goodfellow with ±15% standard tolerance, see also https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  const auto activationTargetSmThickness = 0.19 * mm; // Assumed, calculated as the thickness of a layer of solid Sm2O3 (with literature density) in the target container with its mass being the average of the Sm activation targets' masses
  const auto activationTargetCeThickness = 0.55 * mm; // Assumed, calculated as the thickness of a layer of solid CeO2 (with literature density) in the target container with its mass being the average of the Ce activation targets' masses

  auto activationTargetMass = -1. * g;
  if (TARGET == "Ce_No1") {
    activationTargetMass = (2.35393 - 1.88349) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Ce_No2") {
    activationTargetMass = (2.33186 - 1.88597) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Ce_No3") {
    activationTargetMass = (2.30733 - 1.88050) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Ce_No4") {
    activationTargetMass = (2.31009 - 1.87927) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Ce_No5") {
    activationTargetMass = (2.33911 - 1.88289) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Ce_No6") {
    activationTargetMass = (2.35403 - 1.88406) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Ce_No7") {
    activationTargetMass = (2.37485 - 1.88253) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Sm_No8") {
    activationTargetMass = (2.05587 - 1.88185) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Sm_No9") {
    activationTargetMass = (2.07222 - 1.88257) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Sm_No10") {
    activationTargetMass = (2.11479 - 1.88234) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Sm_No11") {
    activationTargetMass = (2.07366 - 1.88180) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Sm_No12") {
    activationTargetMass = (2.03059 - 1.88428) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Sm_No13") {
    activationTargetMass = (2.09876 - 1.87385) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Sm_No14") {
    activationTargetMass = (2.05948 - 1.88342) * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/829
  } else if (TARGET == "Au_No1") {
    activationTargetMass = 0.05169 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No2") {
    activationTargetMass = 0.04535 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No3") {
    activationTargetMass = 0.04798 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No4") {
    activationTargetMass = 0.04818 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No5") {
    activationTargetMass = 0.04950 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No6") {
    activationTargetMass = 0.04987 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No7") {
    activationTargetMass = 0.05209 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No8") {
    activationTargetMass = 0.04526 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No9") {
    activationTargetMass = 0.04776 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No10") {
    activationTargetMass = 0.04534 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No11") {
    activationTargetMass = 0.04860 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  } else if (TARGET == "Au_No12") {
    activationTargetMass = 0.05168 * g; // From https://elog.ikp.physik.tu-darmstadt.de/clovershare/843
  }

  G4Material *activationTargetMaterial = nullptr;
  if (targetType == "Sm") {
    const auto activationTargetSmDensity = activationTargetMass / (pi / 4. * activationTargetContainerInnerDiameter * activationTargetContainerInnerDiameter * activationTargetSmThickness);
    activationTargetMaterial = new G4Material("activationTargetSmMaterial", activationTargetSmDensity, 2);
    activationTargetMaterial->AddElement(nist->FindOrBuildElement("Sm"), 2); // 2 Sm atoms in Sm(2)O(3)
    activationTargetMaterial->AddElement(nist->FindOrBuildElement("O"), 3); // 3 O atoms in Sm(2)O(3)
  } else if (targetType == "Ce") {
    const auto activationTargetCeDensity = activationTargetMass / (pi / 4. * activationTargetContainerInnerDiameter * activationTargetContainerInnerDiameter * activationTargetCeThickness);
    activationTargetMaterial = new G4Material("activationTargetCeMaterial", activationTargetCeDensity, 2);
    activationTargetMaterial->AddElement(nist->FindOrBuildElement("Ce"), 1); // 1 Ce atoms in CeO(2)
    activationTargetMaterial->AddElement(nist->FindOrBuildElement("O"), 2); // 2 O atoms in CeO(2)
  } else if (targetType == "Au") {
    const auto activationTargetAuDensity = activationTargetMass / (pi / 4. * activationTargetAuDiameter * activationTargetAuDiameter * activationTargetAuThickness);
    activationTargetMaterial = new G4Material("activationTargetAuMaterial", activationTargetAuDensity, nist->FindOrBuildMaterial("G4_Au"));
  }

  // --------------- Mixed Source ---------------

  const auto mixedSourceDiameter = 25.4 * mm; // From data sheet: https://elog.ikp.physik.tu-darmstadt.de/clovershare/17
  const auto mixedSourceThickness = 6. * mm; // From data sheet: https://elog.ikp.physik.tu-darmstadt.de/clovershare/17
  const auto mixedSourceActiveAreaDiameter = 5. * mm; // From data sheet: https://elog.ikp.physik.tu-darmstadt.de/clovershare/17
  const auto mixedSourceActiveAreaThickness = 0.1 * mm; // Assumed

  // --------------- Detector 3 ---------------
  // Detector 3 is of planar geometry (so-called by Canberra a "broad energy germanium" detector) with a thin beryllium window, unfortunately no data sheet is available for it.

  const auto det3HousingOuterDiameter = targetHolderDet3DetectorEnclosingInnerDiameter; // Assumed
  const auto det3HousingLength = 80. * mm; // Assumed
  const auto det3HousingThickness = 1.5 * mm; // Assumed
  const auto det3HousingFaceThickness = 0.5 * mm; // Assumed
  const auto det3BerylliumWindowDiameter = targetHolderDet3DetectorEnclosingInnerDiameter - 2. * det3HousingThickness;
  const auto det3CrystalLength = 25. * mm; // Assumed, Mirion (former Canberra): "With (crystal) cross-sectional areas of 20 to 65 cm2 (2.52 to 4.55cm radius) and thickness' of 20 to 30 mm ..." - https://www.mirion.com/products/bege-broad-energy-germanium-detectors
  const auto det3CrystalDiameter = 70. * mm; // Assumed, Mirion (former Canberra): "With (crystal) cross-sectional areas of 20 to 65 cm2 (2.52 to 4.55cm radius) and thickness' of 20 to 30 mm ..." - https://www.mirion.com/products/bege-broad-energy-germanium-detectors
  const auto det3CrystalFaceToHousingGap = 5. * mm - det3HousingFaceThickness; // S. Finch: "This places the center of the sample 5 cm from the face of the HPGe detector. The face of the crystal will then be 5.5 cm from the sample center."

  // --------------- Detector 4 ---------------
  // Detector 4 is a regular coaxial HPGe detector with an aluminum window, but unfortunately no data sheet is available for it

  const auto det4HousingOuterDiameter = targetHolderDet4DetectorEnclosingInnerDiameter; // Assumed
  const auto det4HousingLength = 100. * mm; // Assumed
  const auto det4HousingThickness = 1.5 * mm; // Assumed
  const auto det4HousingFaceThickness = 1. * mm; // Assumed
  const auto det4CrystalLength = 65. * mm; // Assumed
  const auto det4CrystalDiameter = 68. * mm; // Assumed
  const auto det4CrystalFaceToHousingGap = 5. * mm - det4HousingFaceThickness; // S. Finch: "This places the center of the sample 5 cm from the face of the HPGe detector. The face of the crystal will then be 5.5 cm from the sample center."

  // --------------- z-Positions ---------------

  const auto zOffsetThreeScrewHolderInteriorCenterToWorldOrigin = 0.; // The reference point is chosen to be center of the interior of the Three-Screw Holder
  const auto activationTargetCenterToThreeScrewHolderInteriorCenterOffset = 0. * mm; // Assume for now, that the activation targets (containers) and sources were always centered in the interior of the Three-Screw Holder
  const auto detectorFaceToWorldOrigin = 5. * cm + zOffsetThreeScrewHolderInteriorCenterToWorldOrigin; // The Three-Screw Holder was always placed in the 5 cm slot of the Target Holder. This places the center of the sample (center of the interior of the Three-Screw Holder) 5 cm from the face of the HPGe detector
  const auto threeScrewHolderCenterToWorldOrigin = threeScrewHolderInnerLength / 2. + threeScrewHolderBottomThickness - threeScrewHolderLength / 2. + zOffsetThreeScrewHolderInteriorCenterToWorldOrigin; // Assume that the (thinner) bottom of th Three-Screw Holder faced the detector
  const auto slideCenterToWorldOrigin = threeScrewHolderCenterToWorldOrigin + threeScrewHolderLength / 2. - threeScrewHolderDistanceGrooveCenterToBottom;

  //  =============== Construction ===============

  //  --------------- World ---------------

  auto *worldSolid = new G4Box("worldSolid", worldXLength / 2., worldYWidth / 2., worldZHeight / 2.);
  auto *worldLogical = new G4LogicalVolume(worldSolid, nist->FindOrBuildMaterial("G4_AIR"), "worldLogical");
  G4VPhysicalVolume *worldPhysical = new G4PVPlacement(nullptr, G4ThreeVector(), worldLogical, "world", nullptr, false, 0);
  auto worldVis = G4VisAttributes(red);
  worldVis.SetForceWireframe(true);
  worldLogical->SetVisAttributes(worldVis);

  // --------------- Target Holder ---------------

  auto *targetHolderSlideHolderBottomSolid = new G4Box("targetHolderSlideHolderBottomSolid", (targetHolderSlideHolderWidth - 2. * targetHolderSlideHolderWallThickness) / 2., targetHolderSlideHolderWallThickness / 2., targetHolderSlideHolderLength / 2.);
  auto *targetHolderSlideHolderBottomLogical = new G4LogicalVolume(targetHolderSlideHolderBottomSolid, targetHolderMaterial, "targetHolderSlideHolderBottomLogical");
  new G4PVPlacement(nullptr, G4ThreeVector(0, -targetHolderSlideHolderHeight / 2. + targetHolderSlideHolderWallThickness / 2., detectorFaceToWorldOrigin - targetHolderDetectorEnclosingFrontFaceThickness - targetHolderSlideHolderLength / 2.), targetHolderSlideHolderBottomLogical, "targetHolderSlideHolderBottom", worldLogical, false, 0);
  targetHolderSlideHolderBottomLogical->SetVisAttributes(green);

  auto *targetHolderSlideHolderRightSolid = new G4Box("targetHolderSlideHolderRightSolid", targetHolderSlideHolderWallThickness / 2., targetHolderSlideHolderHeight / 2., targetHolderSlideHolderLength / 2.);
  auto *targetHolderSlideHolderRightLogical = new G4LogicalVolume(targetHolderSlideHolderRightSolid, targetHolderMaterial, "targetHolderSlideHolderRightLogical");
  new G4PVPlacement(nullptr, G4ThreeVector(targetHolderSlideHolderWidth / 2. - targetHolderSlideHolderWallThickness / 2., 0, detectorFaceToWorldOrigin - targetHolderDetectorEnclosingFrontFaceThickness - targetHolderSlideHolderLength / 2.), targetHolderSlideHolderRightLogical, "targetHolderSlideHolderRight", worldLogical, false, 0);

  auto *targetHolderSlideHolderLeftSolid = new G4Box("targetHolderSlideHolderLeftSolid", targetHolderSlideHolderWallThickness / 2., targetHolderSlideHolderHeight / 2., targetHolderSlideHolderLength / 2.);
  auto *targetHolderSlideHolderLeftLogical = new G4LogicalVolume(targetHolderSlideHolderLeftSolid, targetHolderMaterial, "targetHolderSlideHolderLeftLogical");
  new G4PVPlacement(nullptr, G4ThreeVector(-(targetHolderSlideHolderWidth / 2. - targetHolderSlideHolderWallThickness / 2.), 0, detectorFaceToWorldOrigin - targetHolderDetectorEnclosingFrontFaceThickness - targetHolderSlideHolderLength / 2.), targetHolderSlideHolderLeftLogical, "targetHolderSlideHolderLeft", worldLogical, false, 0);

  auto targetHolderSlideHolderSideVis = G4VisAttributes(green);
  targetHolderSlideHolderSideVis.SetForceWireframe(true);
  targetHolderSlideHolderRightLogical->SetVisAttributes(targetHolderSlideHolderSideVis);
  targetHolderSlideHolderLeftLogical->SetVisAttributes(targetHolderSlideHolderSideVis);

  auto *targetHolderDetectorEnclosingFrontSolid = new G4Tubs("targetHolderDetectorEnclosingFrontSolid", targetHolderDetectorEnclosingInnerDiameter / 2., targetHolderDetectorEnclosingFrontOuterDiameter / 2., targetHolderDetectorEnclosingFrontLength / 2., 0., twopi);
  auto *targetHolderDetectorEnclosingFrontLogical = new G4LogicalVolume(targetHolderDetectorEnclosingFrontSolid, targetHolderMaterial, "targetHolderDetectorEnclosingFrontLogical");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, detectorFaceToWorldOrigin - targetHolderDetectorEnclosingFrontFaceThickness + targetHolderDetectorEnclosingFrontLength / 2.), targetHolderDetectorEnclosingFrontLogical, "targetHolderDetectorEnclosingFront", worldLogical, false, 0);
  targetHolderDetectorEnclosingFrontLogical->SetVisAttributes(green);

  auto *targetHolderDetectorEnclosingFrontFaceSolid = new G4Tubs("targetHolderDetectorEnclosingFrontFaceSolid", targetHolderDetectorEnclosingFrontFaceHoleDiameter / 2., targetHolderDetectorEnclosingInnerDiameter / 2., targetHolderDetectorEnclosingFrontFaceThickness / 2., 0., twopi);
  auto *targetHolderDetectorEnclosingFrontFaceLogical = new G4LogicalVolume(targetHolderDetectorEnclosingFrontFaceSolid, targetHolderMaterial, "targetHolderDetectorEnclosingFrontFaceLogical");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, detectorFaceToWorldOrigin - targetHolderDetectorEnclosingFrontFaceThickness / 2.), targetHolderDetectorEnclosingFrontFaceLogical, "targetHolderDetectorEnclosingFrontFace", worldLogical, false, 0);
  targetHolderDetectorEnclosingFrontFaceLogical->SetVisAttributes(green);

  if (useDet4InsteadOfDet3) { // Only the Target Holder for Detector 4 has a enclosing second ring part
    auto *targetHolderDetectorEnclosingBackSolid = new G4Tubs("targetHolderDetectorEnclosingBackSolid", targetHolderDetectorEnclosingInnerDiameter / 2., targetHolderDet4DetectorEnclosingBackOuterDiameter / 2., targetHolderDet4DetectorEnclosingBackLength / 2., 0., twopi);
    auto *targetHolderDetectorEnclosingBackLogical = new G4LogicalVolume(targetHolderDetectorEnclosingBackSolid, targetHolderMaterial, "targetHolderDetectorEnclosingBackLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, detectorFaceToWorldOrigin - targetHolderDetectorEnclosingFrontFaceThickness + targetHolderDetectorEnclosingFrontLength + targetHolderDet4DetectorEnclosingBackLength / 2.), targetHolderDetectorEnclosingBackLogical, "targetHolderDetectorEnclosingBack", worldLogical, false, 0);
    targetHolderDetectorEnclosingBackLogical->SetVisAttributes(green);
  }

  // --------------- Slide ---------------

  auto *slideWithoutHoleSolid = new G4Box("slideWithoutHoleSolid", slideWidth / 2., slideHeight / 2., slideThickness / 2.);
  auto *slideHoleSolid = new G4Tubs("slideHoleSolid", 0., slideHoleDiameter / 2., slideThickness / 2. + subtractionSolidBuffer, 0., twopi);
  auto *slideSolid = new G4SubtractionSolid("slideSolid", slideWithoutHoleSolid, slideHoleSolid, nullptr, G4ThreeVector(0, 0, 0));
  auto *slideLogical = new G4LogicalVolume(slideSolid, acrylic, "slideLogical");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, slideCenterToWorldOrigin), slideLogical, "slide", worldLogical, false, 0);
  slideLogical->SetVisAttributes(lightGrey);

  // --------------- Three-Screw Holder ---------------

  auto *threeScrewHolderWithoutGrooveSolid = new G4Tubs("threeScrewHolderWithoutGrooveSolid", 0., threeScrewHolderOuterDiameter / 2., threeScrewHolderLength / 2., 0, twopi);
  auto *threeScrewHolderGrooveSolid = new G4Tubs("threeScrewHolderGrooveSolid", threeScrewHolderGrooveInnerDiameter / 2., threeScrewHolderOuterDiameter / 2. + subtractionSolidBuffer, threeScrewHolderGrooveWidth / 2., 0, twopi);
  auto *threeScrewHolderSolid = new G4SubtractionSolid("threeScrewHolderSolid", threeScrewHolderWithoutGrooveSolid, threeScrewHolderGrooveSolid, nullptr, G4ThreeVector(0, 0, threeScrewHolderLength / 2. - threeScrewHolderDistanceGrooveCenterToBottom));
  auto *threeScrewHolderLogical = new G4LogicalVolume(threeScrewHolderSolid, acrylic, "threeScrewHolderLogical");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, threeScrewHolderCenterToWorldOrigin), threeScrewHolderLogical, "threeScrewHolder", worldLogical, false, 0);
  threeScrewHolderLogical->SetVisAttributes(grey);

  auto *threeScrewHolderInteriorAirSolid = new G4Tubs("threeScrewHolderInteriorAirSolid", 0., threeScrewHolderInnerDiameter / 2., threeScrewHolderInnerLength / 2., 0, twopi);
  auto *threeScrewHolderInteriorAirLogical = new G4LogicalVolume(threeScrewHolderInteriorAirSolid, nist->FindOrBuildMaterial("G4_AIR"), "threeScrewHolderInteriorAirLogical");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zOffsetThreeScrewHolderInteriorCenterToWorldOrigin - threeScrewHolderCenterToWorldOrigin), threeScrewHolderInteriorAirLogical, "threeScrewHolderInteriorAir", threeScrewHolderLogical, false, 0);
  threeScrewHolderInteriorAirLogical->SetVisAttributes(blue);

  // --------------- Target Container for Sm and Ce ---------------

  if (targetType == "Sm" || targetType == "Ce") {
    const auto activationTargetThickness = (targetType == "Sm") ? activationTargetSmThickness : activationTargetCeThickness;
    const auto activationTargetContainerWithLidLength = std::max(activationTargetContainerLength, activationTargetContainerBottomThickness + activationTargetThickness + activationTargetContainerLidLength);
    const auto activationTargetContainerLidProtrusionLength = activationTargetContainerWithLidLength - activationTargetContainerLength;

    auto *activationTargetContainerTotalExtendedVolumeSolid = new G4Tubs("activationTargetContainerTotalExtendedVolumeSolid", 0., activationTargetContainerOuterDiameter / 2., activationTargetContainerWithLidLength / 2., 0, twopi);
    auto *activationTargetContainerNoLidRegionSolid = new G4Tubs("activationTargetContainerNoLidRegionSolid", activationTargetContainerLidOuterDiameter / 2., activationTargetContainerOuterDiameter / 2. + subtractionSolidBuffer, activationTargetContainerLidProtrusionLength, 0, twopi); // Twice as long as necessary to work as a subtractionSolidBuffer
    auto *activationTargetContainerNoLidRecessSolid = new G4SubtractionSolid("activationTargetContainerNoLidRecessSolid", activationTargetContainerTotalExtendedVolumeSolid, activationTargetContainerNoLidRegionSolid, nullptr, G4ThreeVector(0, 0, activationTargetContainerWithLidLength / 2.));
    auto *activationTargetContainerLidRecessSolid = new G4Tubs("activationTargetContainerLidRecessSolid", 0., activationTargetContainerLidRecessInnerDiameter / 2., activationTargetContainerLidRecessDepth, 0, twopi); // Twice as long as necessary to work as a subtractionSolidBuffer
    auto *activationTargetContainerSolid = new G4SubtractionSolid("activationTargetContainerSolid", activationTargetContainerNoLidRecessSolid, activationTargetContainerLidRecessSolid, nullptr, G4ThreeVector(0, 0, activationTargetContainerWithLidLength / 2.));

    auto *activationTargetContainerLogical = new G4LogicalVolume(activationTargetContainerSolid, nist->FindOrBuildMaterial("G4_POLYETHYLENE"), "activationTargetContainerLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, activationTargetCenterToThreeScrewHolderInteriorCenterOffset), activationTargetContainerLogical, "activationTargetContainer", threeScrewHolderInteriorAirLogical, false, 0);
    activationTargetContainerLogical->SetVisAttributes(white);

    // --------------- Sm/Ce Target ---------------

    auto *activationTargetSolid = new G4Tubs("activationTargetSolid", 0., activationTargetContainerInnerDiameter / 2., activationTargetThickness / 2., 0, twopi);
    auto *activationTargetLogical = new G4LogicalVolume(activationTargetSolid, activationTargetMaterial, "activationTargetLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -activationTargetContainerWithLidLength / 2. + activationTargetContainerBottomThickness + activationTargetThickness / 2.), activationTargetLogical, "activationTarget", activationTargetContainerLogical, false, 0);
    activationTargetLogical->SetVisAttributes(yellow);

    auto *activationTargetIrradiatedPartSolid = new G4Tubs("activationTargetIrradiatedPartSolid", 0., std::min(beamDiameterAtActivationTargetPos, activationTargetContainerInnerDiameter) / 2., activationTargetThickness / 2., 0, twopi);
    auto *activationTargetIrradiatedPartLogical = new G4LogicalVolume(activationTargetIrradiatedPartSolid, activationTargetMaterial, "activationTargetIrradiatedPartLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), activationTargetIrradiatedPartLogical, "activationTargetIrradiatedPart", activationTargetLogical, false, 0);
    activationTargetIrradiatedPartLogical->SetVisAttributes(orange);
  }

  // --------------- Au Target ---------------

  else if (targetType == "Au") {
    auto *activationTargetSolid = new G4Tubs("activationTargetSolid", 0., activationTargetAuDiameter / 2., activationTargetAuThickness / 2., 0, twopi);
    auto *activationTargetLogical = new G4LogicalVolume(activationTargetSolid, activationTargetMaterial, "activationTargetLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, activationTargetCenterToThreeScrewHolderInteriorCenterOffset), activationTargetLogical, "activationTarget", threeScrewHolderInteriorAirLogical, false, 0);
    activationTargetLogical->SetVisAttributes(yellow);

    auto *activationTargetIrradiatedPartSolid = new G4Tubs("activationTargetIrradiatedPartSolid", 0., std::min(beamDiameterAtActivationTargetPos, activationTargetAuDiameter) / 2., activationTargetAuThickness / 2., 0, twopi);
    auto *activationTargetIrradiatedPartLogical = new G4LogicalVolume(activationTargetIrradiatedPartSolid, activationTargetMaterial, "activationTargetIrradiatedPartLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), activationTargetIrradiatedPartLogical, "activationTargetIrradiatedPart", activationTargetLogical, false, 0);
    activationTargetIrradiatedPartLogical->SetVisAttributes(orange);
  }

  // --------------- Mixed Source ---------------
  // According to S. Finch, a correction is necessary for the self-attenuation of the mixed calibration source, as the calibration certificate/data sheet is for the source before the company encapsulates it in its container.
  // Hence the calibration certificate does not accurately represent the source as it is delivered and this attenuation effect has to be taken into account here by constructing the source container in GEANT4.

  else if (targetType == "MixedSource") {
    auto *mixedSourceContainerSolid = new G4Tubs("mixedSourceContainerSolid", 0., mixedSourceDiameter / 2., mixedSourceThickness / 2., 0, twopi);
    auto *mixedSourceContainerLogical = new G4LogicalVolume(mixedSourceContainerSolid, acrylic, "mixedSourceContainerLogical"); // Assume that the source container is acrylic, because this is the case for every other source
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, activationTargetCenterToThreeScrewHolderInteriorCenterOffset), mixedSourceContainerLogical, "mixedSourceContainer", threeScrewHolderInteriorAirLogical, false, 0);
    mixedSourceContainerLogical->SetVisAttributes(white);

    auto *mixedSourceActivePartSolid = new G4Tubs("mixedSourceActivePartSolid", 0., mixedSourceActiveAreaDiameter / 2., mixedSourceActiveAreaThickness / 2., 0, twopi);
    auto *mixedSourceActivePartLogical = new G4LogicalVolume(mixedSourceActivePartSolid, vacuum, "mixedSourceActivePartLogical"); // The self attenuation of the actual activated layer is not simulated, because this should be considered in the calibration certificate activity, hence, use vacuum
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), mixedSourceActivePartLogical, "mixedSourceActivePart", mixedSourceContainerLogical, false, 0);
    mixedSourceActivePartLogical->SetVisAttributes(orange);
  }

  if (!useDet4InsteadOfDet3) {
    // --------------- Detector 3 ---------------

    auto *det3HousingSolid = new G4Tubs("det3HousingSolid", 0, det3HousingOuterDiameter / 2., det3HousingLength / 2., 0., twopi);
    auto *det3HousingLogical = new G4LogicalVolume(det3HousingSolid, nist->FindOrBuildMaterial("G4_Al"), "det3HousingLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, detectorFaceToWorldOrigin + det3HousingLength / 2.), det3HousingLogical, "det3Housing", worldLogical, false, 0);
    det3HousingLogical->SetVisAttributes(darkGrey);

    const auto det3VacuumLength = det3HousingLength - det3HousingFaceThickness - det3HousingThickness;
    auto *det3VacuumSolid = new G4Tubs("det3VacuumSolid", 0, (det3HousingOuterDiameter - 2. * det3HousingThickness) / 2., det3VacuumLength / 2., 0., twopi);
    auto *det3VacuumLogical = new G4LogicalVolume(det3VacuumSolid, vacuum, "det3VacuumLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (det3HousingFaceThickness - det3HousingThickness) / 2.), det3VacuumLogical, "det3Vacuum", det3HousingLogical, false, 0);
    det3VacuumLogical->SetVisAttributes(blue);

    Sensitive_Detector_Logical_Volume_Name = "det3CrystalLogical";
    auto *det3CrystalSolid = new G4Tubs("det3CrystalSolid", 0, det3CrystalDiameter / 2., det3CrystalLength / 2., 0., twopi);
    auto *det3CrystalLogical = new G4LogicalVolume(det3CrystalSolid, nist->FindOrBuildMaterial("G4_Ge"), Sensitive_Detector_Logical_Volume_Name);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -det3VacuumLength / 2. + det3CrystalFaceToHousingGap + det3CrystalLength / 2.), det3CrystalLogical, "det4Crystal", det3VacuumLogical, false, 0);
    det3CrystalLogical->SetVisAttributes(red);

    auto *det3BerylliumWindowSolid = new G4Tubs("det3BerylliumWindowSolid", 0, det3BerylliumWindowDiameter / 2., det3HousingFaceThickness / 2., 0., twopi);
    auto *det3BerylliumWindowLogical = new G4LogicalVolume(det3BerylliumWindowSolid, nist->FindOrBuildMaterial("G4_Be"), "det3BerylliumWindowLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -det3HousingLength / 2. + det3HousingFaceThickness / 2.), det3BerylliumWindowLogical, "det3BerylliumWindow", det3HousingLogical, false, 0);
    det3BerylliumWindowLogical->SetVisAttributes(lightGrey);

  } else {
    // --------------- Detector 4 ---------------

    auto *det4HousingSolid = new G4Tubs("det4HousingSolid", 0, det4HousingOuterDiameter / 2., det4HousingLength / 2., 0., twopi);
    auto *det4HousingLogical = new G4LogicalVolume(det4HousingSolid, nist->FindOrBuildMaterial("G4_Al"), "det4HousingLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, detectorFaceToWorldOrigin + det4HousingLength / 2.), det4HousingLogical, "det4Housing", worldLogical, false, 0);
    det4HousingLogical->SetVisAttributes(darkGrey);

    const auto det4VacuumLength = det4HousingLength - det4HousingFaceThickness - det4HousingThickness;
    auto *det4VacuumSolid = new G4Tubs("det4VacuumSolid", 0, (det4HousingOuterDiameter - 2. * det4HousingThickness) / 2., det4VacuumLength / 2., 0., twopi);
    auto *det4VacuumLogical = new G4LogicalVolume(det4VacuumSolid, vacuum, "det4VacuumLogical");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (det4HousingFaceThickness - det4HousingThickness) / 2.), det4VacuumLogical, "det4Vacuum", det4HousingLogical, false, 0);
    det4VacuumLogical->SetVisAttributes(blue);

    Sensitive_Detector_Logical_Volume_Name = "det4CrystalLogical";
    auto *det4CrystalSolid = new G4Tubs("det4CrystalSolid", 0, det4CrystalDiameter / 2., det4CrystalLength / 2., 0., twopi);
    auto *det4CrystalLogical = new G4LogicalVolume(det4CrystalSolid, nist->FindOrBuildMaterial("G4_Ge"), Sensitive_Detector_Logical_Volume_Name);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -det4VacuumLength / 2. + det4CrystalFaceToHousingGap + det4CrystalLength / 2.), det4CrystalLogical, "det4Crystal", det4VacuumLogical, false, 0);
    det4CrystalLogical->SetVisAttributes(red);
  }

  return worldPhysical;
}

void DetectorConstruction::ConstructSDandField() {
  EnergyDepositionSD *SensitiveDetector = new EnergyDepositionSD(Sensitive_Detector_Logical_Volume_Name, Sensitive_Detector_Logical_Volume_Name);
  G4SDManager::GetSDMpointer()->AddNewDetector(SensitiveDetector);
  SensitiveDetector->SetDetectorID(0);
  SetSensitiveDetector(Sensitive_Detector_Logical_Volume_Name, SensitiveDetector, true);

  Max_Sensitive_Detector_ID = 0; // Necessary for EVENT_EVENTWISE output mode
}
