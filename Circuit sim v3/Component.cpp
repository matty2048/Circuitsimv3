#include "Component.h"

Eigen::MatrixXf Component::BuildRHSMatrix()
{
    return Eigen::MatrixXf();
}

Eigen::MatrixXf Component::BuildLHSMatrix()
{
    return Eigen::MatrixXf();
}

void Component::FinishConstruction()
{
}
