#ifndef _PHYSICSOBJECTS_
#define _PHYSICSOBJECTS_

#include "TreeUtilities.h"
#include "Logging.h"

class PhysicsObject
{
 private:
  /* data */
 public:
  float Px = 0;
  float Py = 0;
  float Pz = 0;
  PhysicsObject(){};
  explicit PhysicsObject(float PxFromTree, float PyFromTree, float PzFromTree); // explicit means this constructor can be inherited by the Gammas, Jets and Pi0s
  ~PhysicsObject(){};
};

PhysicsObject::PhysicsObject(float PxFromTree, float PyFromTree, float PzFromTree)
{
  Px = PxFromTree;
  Py = PyFromTree;
  Pz = PzFromTree;
}

class IsoGamma : public PhysicsObject
{
 public:
  using PhysicsObject::PhysicsObject;
  ~IsoGamma(){};
  float E = 0;
};

void saveClustersFromEventInVector(TreeBuffer tree, std::vector<IsoGamma> IsoGammas)
{
  for (int iCluster = 0; iCluster < (int) tree.Cluster_Px->size(); iCluster++) {
    IsoGamma isoGamma(tree.Cluster_Px->at(iCluster), tree.Cluster_Px->at(iCluster), tree.Cluster_Px->at(iCluster));
    isoGamma.E = tree.Cluster_E->at(iCluster);
  }
}

class Jet : public PhysicsObject
{
 public:
  using PhysicsObject::PhysicsObject;
  ~Jet(){};
};

class Pi0 : public PhysicsObject
{
 public:
  using PhysicsObject::PhysicsObject;
  ~Pi0(){};
};

#endif // _PHYSICSOBJECTS_