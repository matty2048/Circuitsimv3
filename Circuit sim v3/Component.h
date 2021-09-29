#pragma once
#include <Eigen/Dense>
class Solver;
class Component
{
public:
    Solver* solver;
    int NodeAIndex; //the Nodes Index in the nodes vector if -1 is connected to 0v
    int NodeBIndex;

   virtual Eigen::MatrixXf BuildRHSMatrix(); //builds the A matrix for the given component
   virtual Eigen::MatrixXf BuildLHSMatrix();// builds the Z matrix for a given component
   virtual void FinishConstruction();
};