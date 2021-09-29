// Circuit sim v3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>


#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <type_traits>
#include "Component.h"

std::ostream& operator<<(std::ostream& os, const Eigen::MatrixXf& dt)
{

    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    os << dt.format(CleanFmt);
    return os;
}

class Solver
{
public:


    void AddNode(int num = 1)
    {
        n = +num;
    }
    template<typename T, typename ... Params>
    void AddComponent(Params&& ... params)
    {
        T* N = new T(std::forward<Params>(params)...);
        N->solver = this;
        N->FinishConstruction();
        Component* comp = dynamic_cast<Component*>(N);
        Components.push_back(comp);
    }

    void Solve()
    {
        std::cout << n << std::endl;
        std::cout << m << std::endl;

        for (auto i : Components)
        {
            i->BuildLHSMatrix();
        }
    }

    int n = -1, m; //n is the number of Nodes 
                   //m is the number of independant voltage sources

   // std::vector<Node> Nodes; //
    std::vector<Component*> Components;
};


struct Resistor : public Component
{
    Resistor(int NodeA, int NodeB, float R)
    {
        NodeAIndex = NodeA;
        NodeBIndex = NodeB;
        G = 1 / R;
    }
    Eigen::MatrixXf BuildRHSMatrix()
    {
        Eigen::MatrixXf ZMatrix;
        return ZMatrix;
    }
    Eigen::MatrixXf BuildLHSMatrix()  override
    {
        Eigen::MatrixXf AMatrix = Eigen::MatrixXf(solver->m + solver->n, solver->m + solver->n);
        Eigen::MatrixXf GMatrix = Eigen::MatrixXf(solver->n, solver->n);
        
        AMatrix.setZero();
        GMatrix.setZero();
        if (NodeAIndex != -1)
        {
            GMatrix(NodeAIndex, NodeAIndex) = G;
            if (NodeBIndex != -1)
            {
                GMatrix(NodeAIndex, NodeBIndex) = -G;
                GMatrix(NodeBIndex, NodeAIndex) = -G;
            }
        }
        if (NodeBIndex != -1)
        {
            GMatrix(NodeBIndex, NodeBIndex) = G;
        }
        for(int i =0; i<solver->n; i++)
        for (int j = 0; j < solver->n; j++)
        {
            AMatrix(i, j) = GMatrix(i, j);
        }
        
        return AMatrix;
    }
    float G;
};
struct VoltageSource : public Component
{
    VoltageSource(int NodeA, int NodeB, float V) : V(V)
    {
        NodeAIndex = NodeA;
        NodeBIndex = NodeB;
    }
    Eigen::MatrixXf BuildLHSMatrix()
    {
        Eigen::MatrixXf AMatrix = Eigen::MatrixXf(solver->m + solver->n, solver->m + solver->n);
        Eigen::MatrixXf BMatrix = Eigen::MatrixXf(solver->m,solver->n);
        Eigen::MatrixXf CMatrix = Eigen::MatrixXf(solver->n,solver->m);

        BMatrix.setZero();
        CMatrix.setZero();

        if (NodeAIndex != -1)
        {
            BMatrix(ID, NodeAIndex) = 1;
        }
        if (NodeBIndex != -1)
        {
            BMatrix(ID, NodeBIndex) = -1;
        }

        std::cout << BMatrix << std::endl;


        return AMatrix;
    }
    void FinishConstruction()
    {
        ID = this->solver->m;
        this->solver->m += 1;
    }
    float V;
    int ID;
};



int main()
{
    Solver s;
    s.AddNode(2); //Adds the Number of non 0v Nodes
    s.AddComponent<VoltageSource>(-1,0,3);
    s.AddComponent<Resistor>(-1, 1, 3);
    s.AddComponent<Resistor>(1, -1);

    s.Solve();
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
