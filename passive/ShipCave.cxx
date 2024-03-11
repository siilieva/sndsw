/* Generated by Together */
#include "ShipCave.h"

#include "ShipGeoCave.h"                // for ShipGeoCave

void ShipCave::ConstructGeometry()
{
  FairGeoLoader* loader=FairGeoLoader::Instance();
  FairGeoInterface* GeoInterface =loader->getGeoInterface();
  ShipGeoCave* MGeo=new ShipGeoCave();
  MGeo->setGeomFile(GetGeometryFileName());
  GeoInterface->addGeoModule(MGeo);
  Bool_t rc = GeoInterface->readSet(MGeo);
  if ( rc ) { MGeo->create(loader->getGeoBuilder()); }
 
}
ShipCave::ShipCave()
:FairModule()
{
}

ShipCave::~ShipCave()
{

}
ShipCave::ShipCave(const char* name,  const char* Title)
  : FairModule(name ,Title)
{
  world[0] = 0;
  world[1] = 0;
  world[2] = 0;
}