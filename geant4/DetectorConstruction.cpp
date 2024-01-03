#include "DetectorConstruction.h"

#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"

#include "LKG4RunManager.h"
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
//#include "G4UniformMagFieldMessenger.hh"
#include "G4UserLimits.hh"
#include "G4DataInterpolation.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"

using namespace CLHEP;

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction()
{
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  auto runManager = (LKG4RunManager *) G4RunManager::GetRunManager();
  auto par = runManager -> GetPar();
  auto nist = G4NistManager::Instance();

  //Elements
  // G4Elements( name, symbol, Z number, A )
  G4Element* elH = new G4Element("H","H", 1., 1.00794*g/mole);
  G4Element* elHe = new G4Element("He", "He", 2., 4.002602*g/mole);
  G4Element* elC = new G4Element("C", "C", 6., 12.011*g/mole);
  G4Element* elN = new G4Element("N", "N", 7., 14.00674*g/mole);
  G4Element* elO = new G4Element("O", "O", 8., 15.9994*g/mole);
  G4Element* elSi= new G4Element("Si","Si", 14., 28.0855*g/mole);
  G4Element* elCs= new G4Element("Cs","Cs", 55., 132.9055*g/mole);
  G4Element* elI= new G4Element("I","I", 53., 126.9045*g/mole);
  G4Element* elAl= new G4Element("Al","Al", 13., 26.98*g/mole);
 

  //--------fromKEBI/ATTPCDetectorConstruction.cc--------
  G4double Temperature = par -> GetParDouble("temperature");
  G4double STPTemperature = 273.15;
  G4double labTemperature = (STPTemperature + Temperature)*kelvin;
  G4double Pressure = par -> GetParDouble("pressure");
  TString detMatName = par -> GetParString("detMatName");
  G4double iC4H10Ratio = par -> GetParDouble("iC4H10Ratio");


  //-----------------gas definition_from_KEBI----------------------

  G4double densityArGas = 0.001782 *g/cm3*STPTemperature/labTemperature;
  G4double densityMethane = 0.000717 *g/cm3*STPTemperature/labTemperature;
  G4double densityHeGas = 0.000179 *g/cm3*STPTemperature/labTemperature;
  G4double densityIsobutaneGas = 0.0026756 *g/cm3*STPTemperature/labTemperature;
  G4double densityVacuum = 1e-15*g/cm3;

  G4Material* ArGas = new G4Material("ArgonGas", 18, 39.948*g/mole, densityArGas, kStateGas, labTemperature);
  G4Material* Methane = nist -> FindOrBuildMaterial("G4_METHANE");
  G4Material* MethaneGas = new G4Material("MethaneGas", densityMethane, Methane, kStateGas, labTemperature);
  G4Material* HeGas = nist -> FindOrBuildMaterial("G4_He");
  G4Material* IsobutaneGas = new G4Material("iC4H10", densityIsobutaneGas, 2, kStateGas, labTemperature);
  IsobutaneGas -> AddElement(elC, 4);
  IsobutaneGas -> AddElement(elH, 10);

  G4Material *matGas = nullptr;

    if (detMatName == "p10"){
     G4double densityGas = (.9*densityArGas + .1*densityMethane)*Pressure/760.;
     matGas = new G4Material("matP10", densityGas, 2, kStateGas, labTemperature);
     matGas -> AddMaterial(ArGas, 90 *perCent);
     matGas -> AddMaterial(MethaneGas, 10 *perCent);
    }

    else if(detMatName == "4He"){
     G4double densityGas = densityHeGas*Pressure/760.;
     matGas = new G4Material("mat4He", densityGas, 1, kStateGas, labTemperature);
     matGas -> AddMaterial(HeGas, 100 *perCent);
    }

    else if(detMatName == "iC4H10"){
     G4double densityGas = densityIsobutaneGas*Pressure/760.;
     matGas = new G4Material("matiC4H10", densityGas, 1, kStateGas, labTemperature);
     matGas -> AddMaterial(IsobutaneGas, 100 *perCent);
    }

    else if(detMatName == "4He_iC4H10"){
     G4double densityGas = ((densityIsobutaneGas*iC4H10Ratio*0.01)+(densityHeGas*(100.-iC4H10Ratio)*0.01))*Pressure/760.;
     matGas = new G4Material("matHe4_iC4H10", densityGas, 2, kStateGas, labTemperature);
     matGas -> AddMaterial(IsobutaneGas, iC4H10Ratio *perCent);
     matGas -> AddMaterial(HeGas, (100.-iC4H10Ratio) *perCent);
    }

  G4Material *Vacuum_TPC = new G4Material("Vacuum_TPC", densityVacuum, 1, kStateGas, labTemperature);
  Vacuum_TPC -> AddMaterial(HeGas, 100 *perCent);
 //-------------------------------------------------------------------------------------------


  G4double density = universe_mean_density;
  // G4Material( name, density, # of consisting atom, state )
  G4Material* Vacuum = new G4Material("Vacuum", density, 2);
  Vacuum->AddElement(elN, .7);
  Vacuum->AddElement(elO, .3);
  //N(70%), O(30%)
  
  // World
  G4double world_size = 5*m;

  G4Box* solidWorld = new G4Box("World", 0.5*world_size, 0.5*world_size, 0.5*world_size);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Vacuum, "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(0,0,0), logicWorld, "World", 0, false, -1, true);
//(회전,두는 위치,무엇을,이름,?,?,?,?)
//G4Box -> logicWorld와 physWorld를 모두 활용하여 선언

  //Water Material
  G4Material* He4_gas = new G4Material("He4_gas",0.00017*g/cm3,1);
  He4_gas->AddElement(elHe,1);
  //Water->AddElement(elO,1);

  //Water Detector
  G4double box1_size=126.0*cm;

  G4Box* GasBox1 = new G4Box("GasBox1", 0.12*box1_size, 0.12*box1_size, 0.5*box1_size);
  G4LogicalVolume* logicGasBox1 = new G4LogicalVolume(GasBox1, Vacuum, "GasBox1");
  G4VPhysicalVolume* physGasBox1 = new G4PVPlacement(0, G4ThreeVector(0, 0, 62.5*cm), logicGasBox1, "GasBox1", logicWorld, false, 1, true);

  G4Material* C2H4 = new G4Material("C2H4",1.18e-3*g/cm3,2);
  C2H4->AddElement(elC,.33);
  C2H4->AddElement(elH,.67);

  G4Material* Si14_solid = new G4Material("Si14_solid",2.33*g/cm3,1);
  Si14_solid->AddElement(elSi,1);

  G4Material* CsI_solid = new G4Material("CsI_solid",4.51*g/cm3,2);
  CsI_solid->AddElement(elCs,.5);
  CsI_solid->AddElement(elI,.5);

  G4Material* Al13_solid = new G4Material("Al13_solid",2.699*g/cm3,1);
  Al13_solid->AddElement(elAl,1);

  G4double C2H4_size = 7.0*cm;
  G4double C2H4_thickness = 0.2*cm;
  G4Box* Target = new G4Box("Target", 0.5*C2H4_size, 0.5*C2H4_size, 0.5*C2H4_thickness);
  G4LogicalVolume* logicTarget = new G4LogicalVolume(Target, C2H4, "Target");
  G4VPhysicalVolume* physTarget = new G4PVPlacement(0, G4ThreeVector(0, 0, -62.4*cm), logicTarget, "Target", logicGasBox1, false, 10, true);

  G4double Si_size = 7.5*cm;
  G4double Si_thickness = 0.1*cm;
  G4Box* Board_Si1 = new G4Box("Board_Si1", 0.27*Si_size, 0.5*Si_size, 0.5*Si_thickness);
  G4LogicalVolume* logicBoard_Si1 = new G4LogicalVolume(Board_Si1, Si14_solid, "Board_Si1");
  //G4VPhysicalVolume* physBoard_Si1 = new G4PVPlacement(0, G4ThreeVector(6.3*cm, 0, 71.585*cm), logicBoard_Si1, "Board_Si1", logicGasBox1, false, 4, true);
  G4VPhysicalVolume* physBoard_Si1 = new G4PVPlacement(0, G4ThreeVector(6.3*cm, 0, 9.085*cm), logicBoard_Si1, "Board_Si1", logicGasBox1, false, 4, true);

  //G4double Si_thickness2 = 1.5*mm;
  G4Box* Board_Si2 = new G4Box("Board_Si2", 0.27*Si_size, 0.5*Si_size, 0.5*Si_thickness);
  G4LogicalVolume* logicBoard_Si2 = new G4LogicalVolume(Board_Si2, Si14_solid, "Board_Si2");
  //G4VPhysicalVolume* physBoard_Si2 = new G4PVPlacement(0, G4ThreeVector(-6.3*cm, 0, 71.585*cm), logicBoard_Si2, "Board_Si2", logicGasBox1, false, 5, true);
  G4VPhysicalVolume* physBoard_Si2 = new G4PVPlacement(0, G4ThreeVector(-6.3*cm, 0, 9.085*cm), logicBoard_Si2, "Board_Si2", logicGasBox1, false, 5, true);


  Si_size = 7.5*cm;
  Si_thickness = 0.1*cm;

  G4double PositionThetaInDegree = 40;
  G4double PositionThetaInRad = PositionThetaInDegree/180.*3.14159;

  G4RotationMatrix* RotateDetector = new G4RotationMatrix();
  RotateDetector->rotateY(-PositionThetaInRad);

    //G4ThreeVector OriginalPosition(17.5*cm, 0, 67.5*cm);
    G4ThreeVector OriginalPosition(12.5*cm, 0, 67.5*cm);
    G4ThreeVector RevisedPosition = OriginalPosition.rotateY(PositionThetaInRad);
    RevisedPosition.setZ(RevisedPosition.z());

  //G4double box2_size=50*cm;
  G4double box2_size=40*cm;

  G4Box* GasBox2 = new G4Box("GasBox2", 0.5*box2_size, 0.5*box2_size, 0.5*box2_size);
  G4LogicalVolume* logicGasBox2 = new G4LogicalVolume(GasBox2, Vacuum, "GasBox2");
  G4VPhysicalVolume* physGasBox2 = new G4PVPlacement(RotateDetector, RevisedPosition, logicGasBox2, "GasBox2", logicWorld, false, 2, true);

  //RotateDetector = new G4RotationMatrix();
  //RotateDetector->rotateX(PositionThetaInRad);

  G4Box* Board_Si3 = new G4Box("Board_Si3", 0.5*Si_size, 0.27*Si_size, 0.5*Si_thickness);
  G4LogicalVolume* logicBoard_Si3 = new G4LogicalVolume(Board_Si3, Si14_solid, "Board_Si3");
  G4VPhysicalVolume* physBoard_Si3 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(-11.1*cm,2.77*cm,4.01*cm), logicBoard_Si3, "Board_Si3", logicGasBox2, false, 6, true);

  G4Box* Board_Si5 = new G4Box("Board_Si5", 0.5*Si_size, 0.27*Si_size, 0.5*Si_thickness);
  G4LogicalVolume* logicBoard_Si5 = new G4LogicalVolume(Board_Si5, Si14_solid, "Board_Si5");
  G4VPhysicalVolume* physBoard_Si5 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(-11.1*cm,-2.77*cm,4.01*cm), logicBoard_Si5, "Board_Si5", logicGasBox2, false, 7, true);


  G4double CsI_size = 5.0*cm;
  G4double CsI_thickness = 5.0*cm;

  G4Box* Board_CsI1 = new G4Box("Board_CsI1", 1*CsI_size, 0.5*CsI_size, 0.5*CsI_thickness);
  G4LogicalVolume* logicBoard_CsI1 = new G4LogicalVolume(Board_CsI1, CsI_solid, "Board_CsI1");
  G4VPhysicalVolume* physBoard_CsI1 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(-12.225*cm,2.98*cm,7.865*cm), logicBoard_CsI1, "Board_CsI1", logicGasBox2, false, 8, true);

  G4Box* Board_CsI3 = new G4Box("Board_CsI3", 1*CsI_size, 0.5*CsI_size, 0.5*CsI_thickness);
  G4LogicalVolume* logicBoard_CsI3 = new G4LogicalVolume(Board_CsI3, CsI_solid, "Board_CsI3");
  G4VPhysicalVolume* physBoard_CsI3 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(-12.225*cm,-2.98*cm,7.865*cm), logicBoard_CsI1, "Board_CsI3", logicGasBox2, false, 19, true);


  //Al13_solid
  G4double Al_size1 = 5.1*cm;
  G4double Al_thickness1 = 10.1*cm;

  G4double Al_size2 = 5.0*cm;
  G4double Al_thickness2 = 5.0*cm;

  //Al_outerBox
  G4Box* Box_Al1 = new G4Box("Box_Al1", 0.5*Al_thickness1, 0.5*Al_size1, 0.5*Al_size1);
  G4LogicalVolume* logicBox_Al1 = new G4LogicalVolume(Box_Al1, Al13_solid, "Box_Al1");

  G4Box* Box_Al3 = new G4Box("Box_Al3", 0.5*Al_thickness1, 0.5*Al_size1, 0.5*Al_size1);
  G4LogicalVolume* logicBox_Al3 = new G4LogicalVolume(Box_Al3, Al13_solid, "Box_Al3");

  //Al_innerBox
  G4Box* Box_inAl1 = new G4Box("Box_inAl1", 0.5*Al_size2, 0.5*Al_size2, 0.5*Al_thickness2);
  G4LogicalVolume* logicBox_inAl1 = new G4LogicalVolume(Box_inAl1, Vacuum, "Box_inAl1");

  G4Box* Box_inAl3 = new G4Box("Box_inAl3", 0.5*Al_size2, 0.5*Al_size2, 0.5*Al_thickness2);
  G4LogicalVolume* logicBox_inAl3 = new G4LogicalVolume(Box_inAl3, Vacuum, "Box_inAl3");
 
  //Create subtraction solid
  G4VSolid* solidBox_Al1 = logicBox_Al1 -> GetSolid();
  G4VSolid* solidBox_inAl1 = logicBox_inAl1 -> GetSolid();
  G4VSolid* solidBox_Al3 = logicBox_Al3 -> GetSolid();
  G4VSolid* solidBox_inAl3 = logicBox_inAl3 -> GetSolid();

  //G4SubtractionSolid* subtractionSolid_Al1 = new G4SubtractionSolid("subtractionSolid_Al1", solidBox_Al1, solidBox_inAl1);
  //G4SubtractionSolid* subtractionSolid_Al3 = new G4SubtractionSolid("subtractionSolid_Al3", solidBox_Al3, solidBox_inAl3);

  G4VSolid* subtractionSolid_Al1 = new G4SubtractionSolid("subtractionSolid_Al1", solidBox_Al1, solidBox_inAl1);
  G4VSolid* subtractionSolid_Al3 = new G4SubtractionSolid("subtractionSolid_Al3", solidBox_Al3, solidBox_inAl3);

  G4LogicalVolume* logicsubtractionSolid_Al1 = new G4LogicalVolume(subtractionSolid_Al1, Al13_solid, "subtractionSolid_Al1");
  G4LogicalVolume* logicsubtractionSolid_Al3 = new G4LogicalVolume(subtractionSolid_Al3, Al13_solid, "subtractionSolid_Al3");

  //Place the subtration solid in the world volume
  G4VPhysicalVolume* physBox_Al1 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(-12.225*cm,2.98*cm,7.865*cm), logicsubtractionSolid_Al1, "subtractionSolid_Al1", logicGasBox2, false, 29, true);
  G4VPhysicalVolume* physBox_Al3 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(-12.225*cm,-2.98*cm,7.865*cm), logicsubtractionSolid_Al3, "subtractionSolid_Al3", logicGasBox2, false, 31, true);

{
  G4VisAttributes * coloredSolid_Al1 = new G4VisAttributes(G4Colour(G4Colour::Yellow()));
  coloredSolid_Al1 -> SetForceWireframe(true);
  logicsubtractionSolid_Al1 -> SetVisAttributes(coloredSolid_Al1);

  G4VisAttributes * coloredSolid_Al3 = new G4VisAttributes(G4Colour(G4Colour::Yellow()));
  coloredSolid_Al3 -> SetForceWireframe(true);
  logicsubtractionSolid_Al3 -> SetVisAttributes(coloredSolid_Al3);
}


  G4double TPC_size = 10.0*cm;
  G4double TPC_thickness = 15.0*cm;

  G4Box* TPC1 = new G4Box("TPC1", 0.5*TPC_size, 0.5*TPC_size, 0.5*TPC_thickness);
  G4LogicalVolume* logicTPC1 = new G4LogicalVolume(TPC1, matGas, "TPC1");
  {
        G4VisAttributes* attTPC = new G4VisAttributes(G4Colour(G4Colour::Gray()));
        attTPC -> SetForceWireframe(true);
        logicTPC1 -> SetVisAttributes(attTPC);
   }
  logicTPC1 -> SetUserLimits(new G4UserLimits(0.1*nm, 0.1*nm));
  G4VPhysicalVolume* physTPC1 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(-12.5*cm,0,-12.5*cm), logicTPC1, "TPC1", logicGasBox2, false, 15, true);
  runManager -> SetSensitiveDetector(physTPC1);
  //--------------------------------------------------------24.01.02
  RotateDetector = new G4RotationMatrix();
  RotateDetector->rotateY(PositionThetaInRad);

    //OriginalPosition = G4ThreeVector(-17.5*cm, 0, 67.5*cm);
    OriginalPosition = G4ThreeVector(-12.5*cm, 0, 67.5*cm);
    RevisedPosition = OriginalPosition.rotateY(-PositionThetaInRad);
    RevisedPosition.setZ(RevisedPosition.z());

  G4Box* GasBox3 = new G4Box("GasBox3", 0.5*box2_size, 0.5*box2_size, 0.5*box2_size);
  G4LogicalVolume* logicGasBox3 = new G4LogicalVolume(GasBox3, Vacuum, "GasBox3");
  G4VPhysicalVolume* physGasBox3 = new G4PVPlacement(RotateDetector, RevisedPosition, logicGasBox3, "GasBox3", logicWorld, false, 3, true);


//  RotateDetector = new G4RotationMatrix();
//  RotateDetector->rotateY(-PositionThetaInRad);

    //OriginalPosition = G4ThreeVector((-12.5*cm, 0, 67.5*cm);
    //RevisedPosition = OriginalPosition.rotateY(PositionThetaInRad);
    //RevisedPosition.setZ(RevisedPosition.z());

  G4Box* Board_Si4 = new G4Box("Board_Si4", 0.5*Si_size, 0.27*Si_size, 0.5*Si_thickness);
  G4LogicalVolume* logicBoard_Si4= new G4LogicalVolume(Board_Si4, Si14_solid, "Board_Si4");
  //G4VPhysicalVolume* physBoard_Si4 = new G4PVPlacement(RotateDetector, RevisedPosition, logicBoard_Si4, "Board_Si4", logicGasBox3, false, 7, true);
  G4VPhysicalVolume* physBoard_Si4 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(11.1*cm,2.77*cm,4.01*cm), logicBoard_Si4, "Board_Si4", logicGasBox3, false, 7, true);

  G4Box* Board_Si6 = new G4Box("Board_Si6", 0.5*Si_size, 0.27*Si_size, 0.5*Si_thickness);
  G4LogicalVolume* logicBoard_Si6= new G4LogicalVolume(Board_Si6, Si14_solid, "Board_Si6");
  G4VPhysicalVolume* physBoard_Si6 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(11.1*cm,-2.77*cm,4.01*cm), logicBoard_Si4, "Board_Si6", logicGasBox3, false, 18, true);


//  RotateDetector = new G4RotationMatrix();
//  RotateDetector->rotateX(-PositionThetaInRad);

    //OriginalPosition = G4ThreeVector(0, 0, 75.365*cm);
    //RevisedPosition = OriginalPosition.rotateX(PositionThetaInRad);
    //RevisedPosition.setZ(RevisedPosition.z());

  CsI_size = 5.0*cm;
  CsI_thickness = 5.0*cm;

  G4Box* Board_CsI2 = new G4Box("Board_CsI2", 1*CsI_size, 0.5*CsI_size, 0.5*CsI_thickness);
  G4LogicalVolume* logicBoard_CsI2 = new G4LogicalVolume(Board_CsI2, CsI_solid, "Board_CsI2");
  G4VPhysicalVolume* physBoard_CsI2 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(12.225*cm,2.98*cm,7.865*cm), logicBoard_CsI2, "Board_CsI2", logicGasBox3, false, 9, true);

  G4Box* Board_CsI4 = new G4Box("Board_CsI4", 1*CsI_size, 0.5*CsI_size, 0.5*CsI_thickness);
  G4LogicalVolume* logicBoard_CsI4 = new G4LogicalVolume(Board_CsI4, CsI_solid, "Board_CsI4");
  G4VPhysicalVolume* physBoard_CsI4 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(12.225*cm,-2.98*cm,7.865*cm), logicBoard_CsI2, "Board_CsI4", logicGasBox3, false, 20, true);


  //Al13_solid
  //G4double Al_size1 = 5.1*cm;
  //G4double Al_thickness1 = 10.1*cm;

  //G4double Al_size2 = 5.0*cm;
  //G4double Al_thickness2 = 5.0*cm;

  //Al_outerBox
  G4Box* Box_Al2 = new G4Box("Box_Al2", 0.5*Al_thickness1, 0.5*Al_size1, 0.5*Al_size1);
  G4LogicalVolume* logicBox_Al2 = new G4LogicalVolume(Box_Al2, Al13_solid, "Box_Al2");

  G4Box* Box_Al4 = new G4Box("Box_Al4", 0.5*Al_thickness1, 0.5*Al_size1, 0.5*Al_size1);
  G4LogicalVolume* logicBox_Al4 = new G4LogicalVolume(Box_Al4, Al13_solid, "Box_Al4");

  //Al_innerBox
  G4Box* Box_inAl2 = new G4Box("Box_inAl2", 0.5*Al_size2, 0.5*Al_size2, 0.5*Al_thickness2);
  G4LogicalVolume* logicBox_inAl2 = new G4LogicalVolume(Box_inAl2, Vacuum, "Box_inAl2");

  G4Box* Box_inAl4 = new G4Box("Box_inAl4", 0.5*Al_size2, 0.5*Al_size2, 0.5*Al_thickness2);
  G4LogicalVolume* logicBox_inAl4 = new G4LogicalVolume(Box_inAl4, Vacuum, "Box_inAl4");

  //Create subtraction solid
  G4VSolid* solidBox_Al2 = logicBox_Al2 -> GetSolid();
  G4VSolid* solidBox_inAl2 = logicBox_inAl2 -> GetSolid();
  G4VSolid* solidBox_Al4 = logicBox_Al4 -> GetSolid();
  G4VSolid* solidBox_inAl4 = logicBox_inAl4 -> GetSolid();

  //G4SubtractionSolid* subtractionSolid_Al2 = new G4SubtractionSolid("subtractionSolid_Al2", solidBox_Al2, solidBox_inAl2);
  //G4SubtractionSolid* subtractionSolid_Al4 = new G4SubtractionSolid("subtractionSolid_Al4", solidBox_Al4, solidBox_inAl4);

  G4VSolid* subtractionSolid_Al2 = new G4SubtractionSolid("subtractionSolid_Al2", solidBox_Al2, solidBox_inAl2);
  G4VSolid* subtractionSolid_Al4 = new G4SubtractionSolid("subtractionSolid_Al4", solidBox_Al4, solidBox_inAl4);

  G4LogicalVolume* logicsubtractionSolid_Al2 = new G4LogicalVolume(subtractionSolid_Al2, Al13_solid, "subtractionSolid_Al2");
  G4LogicalVolume* logicsubtractionSolid_Al4 = new G4LogicalVolume(subtractionSolid_Al4, Al13_solid, "subtractionSolid_Al4");

{
  G4VisAttributes * coloredSolid_Al2 = new G4VisAttributes(G4Colour(G4Colour::Yellow()));
  coloredSolid_Al2 -> SetForceWireframe(true);
  logicsubtractionSolid_Al2 -> SetVisAttributes(coloredSolid_Al2);

  G4VisAttributes * coloredSolid_Al4 = new G4VisAttributes(G4Colour(G4Colour::Yellow()));
  coloredSolid_Al4 -> SetForceWireframe(true);
  logicsubtractionSolid_Al4 -> SetVisAttributes(coloredSolid_Al4);
}

  //Place the subtration solid in the world volume
  G4VPhysicalVolume* physBox_Al2 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(12.225*cm,2.98*cm,7.865*cm), logicsubtractionSolid_Al2, "subtractionSolid_Al1", logicGasBox3, false, 30, true);
  G4VPhysicalVolume* physBox_Al4 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(12.225*cm,-2.98*cm,7.865*cm), logicsubtractionSolid_Al4, "subtractionSolid_Al3", logicGasBox3, false, 32, true);


  G4Box* TPC2 = new G4Box("TPC2", 0.5*TPC_size, 0.5*TPC_size, 0.5*TPC_thickness);
  G4LogicalVolume* logicTPC2 = new G4LogicalVolume(TPC2, matGas, "TPC2");
  {                                                                                                                                               
        G4VisAttributes* attTPC = new G4VisAttributes(G4Colour(G4Colour::Gray()));
        attTPC -> SetForceWireframe(true);
        logicTPC1 -> SetVisAttributes(attTPC);
   }
  logicTPC2 -> SetUserLimits(new G4UserLimits(0.1*nm, 0.1*nm));
  G4VPhysicalVolume* physTPC2 = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(12.5*cm,0,-12.5*cm), logicTPC2, "TPC2", logicGasBox3, false, 16, true);
  runManager -> SetSensitiveDetector(physTPC2);

  G4double BDC_size=16.0*cm;

  G4Box* BDC = new G4Box("BDC", 0.5*BDC_size, 0.5*BDC_size, 0.3*BDC_size);
  G4LogicalVolume* logicBDC = new G4LogicalVolume(BDC, Vacuum, "BDC");
  G4VPhysicalVolume* physBDC = new G4PVPlacement(0, G4ThreeVector(0, 0, 208.3*cm), logicBDC, "BDC", logicWorld, false, 11, true);


  G4double SC_size1=40.0*cm;
  G4double SC_thickness1=200*um;

  G4Box* SC1 = new G4Box("SC1", 0.5*SC_size1, 0.5*SC_size1, 0.5*SC_thickness1);
  G4LogicalVolume* logicSC1 = new G4LogicalVolume(SC1, He4_gas, "SC1");
  G4VPhysicalVolume* physSC1 = new G4PVPlacement(0, G4ThreeVector(0, 0, -12.3*cm), logicSC1, "SC1", logicWorld, false, 12, true);

  G4double SC_size2=50.0*cm;
  G4double SC_thickness2=0.3*cm;

  G4Box* SC2 = new G4Box("SC2", 0.5*SC_size2, 0.5*SC_size2, 0.5*SC_thickness1);
  G4LogicalVolume* logicSC2 = new G4LogicalVolume(SC2, He4_gas, "SC2");
  G4VPhysicalVolume* physSC2 = new G4PVPlacement(0, G4ThreeVector(0, 0, 196.3*cm), logicSC2, "SC2", logicWorld, false, 13, true);

  G4Box* SC3 = new G4Box("SC3", 0.5*SC_size2, 0.5*SC_size2, 0.5*SC_thickness2);
  G4LogicalVolume* logicSC3 = new G4LogicalVolume(SC3, He4_gas, "SC3");
  G4VPhysicalVolume* physSC3 = new G4PVPlacement(0, G4ThreeVector(0, 0, 220.3*cm), logicSC3, "SC3", logicWorld, false, 14, true);

  return physWorld;
}
