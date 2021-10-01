// Circuit sim v3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>


#include <iostream>
#include <glm/glm.hpp>
#include <Eigen/Dense>

#include <implot.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <gl/glew.h>
#include <GLFW/glfw3.h>

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
    bool first = true;
    void Solve()
    {
        t += dt;
        if (first)
        {
            Last = Eigen::MatrixXf(m + n, 1);
            Last.setZero();
            first = false;
        }
       Eigen::MatrixXf AMatrix = Eigen::MatrixXf(n + m, n + m);
       Eigen::MatrixXf ZMatrix = Eigen::MatrixXf(m + n, 1);
       Eigen::MatrixXf XMatrix = Eigen::MatrixXf(m + n, 1);
      
       AMatrix.setZero();
       ZMatrix.setZero();
       XMatrix.setZero();


       for (auto i : Components)
       {
           AMatrix += i->BuildLHSMatrix();
           ZMatrix += i->BuildRHSMatrix();
       }
       
       XMatrix = AMatrix.inverse() * ZMatrix;
       std::cout << AMatrix << std::endl;
       std::cout << std::endl << std::endl;
       std::cout << XMatrix << std::endl;
       std::cout << std::endl;
      // std::cout << std::endl;
       Last = XMatrix;
    }

    int n = -1, m; //n is the number of Nodes 
                   //m is the number of independant voltage sources
    int c;
    float dt = 0.1;
    float t = 0;
   // std::vector<Node> Nodes; //
    Eigen::MatrixXf Last;
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
       Eigen::MatrixXf ZMatrix = Eigen::MatrixXf(solver->n + solver->m,1);
       ZMatrix.setZero();
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
struct CurrentSource : public Component
{
    CurrentSource(int NodeA, int NodeB, float I) :I(I)
    {
        NodeAIndex = NodeA;
        NodeBIndex = NodeB;
    }

    Eigen::MatrixXf BuildLHSMatrix()
    {
        Eigen::MatrixXf AMatrix = Eigen::MatrixXf(solver->n + solver->m, solver->n + solver->m);
        AMatrix.setZero();
        return AMatrix;
    }

    Eigen::MatrixXf BuildRHSMatrix()
    {
        Eigen::MatrixXf ZMatrix = Eigen::MatrixXf(solver->n + solver->m, 1);
        ZMatrix.setZero();
        if(NodeAIndex != -1)
        ZMatrix(NodeAIndex) += I;
        if (NodeBIndex != -1)
            ZMatrix(NodeBIndex) -= I;

        return ZMatrix;

    }
    int I;
};
struct VoltageSource : public Component
{
    VoltageSource(int NodeA, int NodeB, float V) : V(V)
    {
        NodeAIndex = NodeA;
        NodeBIndex = NodeB;
    }
    void FinishConstruction()
    {
        ID = this->solver->m;
        this->solver->m += 1;
    }
    Eigen::MatrixXf BuildRHSMatrix()
    {
      Eigen::MatrixXf ZMatrix = Eigen::MatrixXf(solver->m + solver->n,1);
      ZMatrix.setZero();
      ZMatrix(solver->n + ID,0) = V;

      return ZMatrix;
    }
    Eigen::MatrixXf BuildLHSMatrix()
    {
        Eigen::MatrixXf AMatrix = Eigen::MatrixXf(solver->m + solver->n, solver->m + solver->n);
        Eigen::MatrixXf BMatrix = Eigen::MatrixXf(solver->m,solver->n);
        Eigen::MatrixXf CMatrix = Eigen::MatrixXf(solver->n,solver->m);
        
        AMatrix.setZero();
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
        
        //std::cout << BMatrix << std::endl;
        
        CMatrix = BMatrix.transpose();
        
        for (int i = solver->n; i < solver->n + solver->m; i++)
        {
            for (int j = 0; j < solver->n; j++)
            {
                AMatrix(i, j) = BMatrix(i - solver->n, j);
            }
        }
        
        for (int j = solver->n; j < solver->n + solver->m; j++)
        {
            for (int i = 0; i < solver->n; i++)
            {
                AMatrix(i, j) = CMatrix(i, j-solver->n);
            }
        }
        
        //std::cout << AMatrix << std::endl;
        return AMatrix;
    }
    float V;
    int ID;
};
struct ACSource : public Component
{
    ACSource(int NodeA, int NodeB, float(*V)(float) ) : V(V)
    {
        NodeAIndex = NodeA;
        NodeBIndex = NodeB;
    }
    Eigen::MatrixXf BuildRHSMatrix()
    {
        Eigen::MatrixXf ZMatrix = Eigen::MatrixXf(solver->m + solver->n, 1);
        ZMatrix.setZero();
        ZMatrix(solver->n + ID, 0) = V(solver->t);

        return ZMatrix;
    }
    Eigen::MatrixXf BuildLHSMatrix()
    {
        Eigen::MatrixXf AMatrix = Eigen::MatrixXf(solver->m + solver->n, solver->m + solver->n);
        Eigen::MatrixXf BMatrix = Eigen::MatrixXf(solver->m, solver->n);
        Eigen::MatrixXf CMatrix = Eigen::MatrixXf(solver->n, solver->m);

        AMatrix.setZero();
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

        //std::cout << BMatrix << std::endl;

        CMatrix = BMatrix.transpose();

        for (int i = solver->n; i < solver->n + solver->m; i++)
        {
            for (int j = 0; j < solver->n; j++)
            {
                AMatrix(i, j) = BMatrix(i - solver->n, j);
            }
        }

        for (int j = solver->n; j < solver->n + solver->m; j++)
        {
            for (int i = 0; i < solver->n; i++)
            {
                AMatrix(i, j) = CMatrix(i, j - solver->n);
            }
        }

        //std::cout << AMatrix << std::endl;
        return AMatrix;
    }
    void FinishConstruction()
    {
        ID = this->solver->m;
        this->solver->m += 1;
    }
    float(*V)(float);
    int ID;
};
struct Inductor : public Component
{
    Inductor(int NodeA, int NodeB, float L): L(L)
    {
        NodeAIndex = NodeA;
        NodeBIndex = NodeB;
         
    }
    void FinishConstruction()
    {
        ID = this->solver->m;
        Req = L/solver->dt;
        solver->m += 1;
    }
    Eigen::MatrixXf BuildLHSMatrix()
    {
        Eigen::MatrixXf AMatrix = Eigen::MatrixXf(solver->m + solver->n, solver->m + solver->n);
        AMatrix.setZero();
        AMatrix(solver->n + ID,solver->n + ID) -= Req;
        if (NodeAIndex != -1) {
            AMatrix(solver->n + ID, NodeAIndex) = 1;
            AMatrix(NodeAIndex, solver->n + ID) = 1;
        }
        if (NodeBIndex != -1) {
            AMatrix(solver->n + ID, NodeBIndex) = -1;
            AMatrix(NodeBIndex, solver->n + ID) = -1;
        }
      //  std::cout << AMatrix << std::endl;
      //  std::cout << std::endl;
        return AMatrix;
    }

    Eigen::MatrixXf BuildRHSMatrix()
    {
        Eigen::MatrixXf ZMatrix = Eigen::MatrixXf(solver->m + solver->n, 1);
        ZMatrix.setZero(); 
        V = -Req * solver->Last(solver->n + ID, 0);
        ZMatrix(solver->n + ID, 0) = V;
         
//        std::cout << ZMatrix << std::endl;
        return ZMatrix;
    }
    float V;
    float Req;
    float L;
    int ID;
};
struct Capacitor : public Component
{
    Capacitor(int NodeA, int NodeB, float C):C(C)
    {
        NodeAIndex = NodeA;
        NodeBIndex = NodeB;
    }
    void FinishConstruction()
    {
        Geq = C/solver->dt;
        
    }
    Eigen::MatrixXf BuildLHSMatrix()
    {
        Eigen::MatrixXf AMatrix = Eigen::MatrixXf(solver->m + solver->n, solver->m + solver->n);
        Eigen::MatrixXf GMatrix = Eigen::MatrixXf(solver->n, solver->n);

        AMatrix.setZero();
        GMatrix.setZero();
        if (NodeAIndex != -1)
        {
            GMatrix(NodeAIndex, NodeAIndex) = Geq;
            if (NodeBIndex != -1)
            {
                GMatrix(NodeAIndex, NodeBIndex) = -Geq;
                GMatrix(NodeBIndex, NodeAIndex) = -Geq;
            }
        }
        if (NodeBIndex != -1)
        {
            GMatrix(NodeBIndex, NodeBIndex) = Geq;
        }
        for (int i = 0; i < solver->n; i++)
            for (int j = 0; j < solver->n; j++)
            {
                AMatrix(i, j) = GMatrix(i, j);
            }

        return AMatrix;

    }
    Eigen::MatrixXf BuildRHSMatrix()
    {
        Eigen::MatrixXf ZMatrix = Eigen::MatrixXf(solver->n + solver->m, 1);
        ZMatrix.setZero();

        if (NodeAIndex != -1 && NodeBIndex != -1)
        {
            ZMatrix(NodeAIndex,0) += Geq * solver->Last(NodeAIndex, 0) - solver->Last(NodeBIndex, 0);
            ZMatrix(NodeBIndex,0) += Geq * solver->Last(NodeAIndex, 0) - solver->Last(NodeBIndex, 0);
        }
        else
            if (NodeAIndex != -1)
                ZMatrix(NodeAIndex,0) += Geq * solver->Last(NodeAIndex, 0);
        else
        if (NodeBIndex != -1)
            ZMatrix(NodeBIndex,0) -= Geq * solver->Last(NodeBIndex, 0);
        return ZMatrix;

    }
    float C;
    float Geq;
};
struct OpAmp : public Component
{
    OpAmp(int NodeA, int NodeB, int NodeC) // A is + B is - C is output
    {
        NodeAIndex = NodeA;
        NodeBIndex = NodeB;
        NodeCIndex = NodeC;
    }
    void FinishConstruction()
    {
        ID = this->solver->m;
        this->solver->m += 1;
    }
    Eigen::MatrixXf BuildRHSMatrix()
    {
        Eigen::MatrixXf ZMatrix = Eigen::MatrixXf(solver->m + solver->n, 1);
        ZMatrix.setZero();
       

        return ZMatrix;
    }
    Eigen::MatrixXf BuildLHSMatrix()
    {
        Eigen::MatrixXf AMatrix = Eigen::MatrixXf(solver->m + solver->n, solver->m + solver->n);
        Eigen::MatrixXf BMatrix = Eigen::MatrixXf(solver->m, solver->n);
        Eigen::MatrixXf CMatrix = Eigen::MatrixXf(solver->n, solver->m);

        AMatrix.setZero();
        BMatrix.setZero();
        CMatrix.setZero();

        if (NodeCIndex != -1)
        {
            CMatrix(NodeCIndex,ID) = 1;
        }

        //std::cout << BMatrix << std::endl;

        //CMatrix = BMatrix.transpose();
        if(NodeAIndex != -1)
        BMatrix(ID,NodeAIndex) = 1;
        if(NodeBIndex != -1)
        BMatrix(ID,NodeBIndex) = -1;


        for (int i = solver->n; i < solver->n + solver->m; i++)
        {
            for (int j = 0; j < solver->n; j++)
            {
                AMatrix(i, j) = BMatrix(i - solver->n, j);
            }
        }

        for (int j = solver->n; j < solver->n + solver->m; j++)
        {
            for (int i = 0; i < solver->n; i++)
            {
                AMatrix(i, j) = CMatrix(i, j - solver->n);
            }
        }

        //std::cout << AMatrix << std::endl;
        return AMatrix;
    }

    int ID;
    int NodeCIndex;
};


int main()
{


    if (!glfwInit())
        return -1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);  //initializes GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(1000, 1000, "ViewPort", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    //Creates GLFW window for drawing to
    glfwMakeContextCurrent(window);

    if (glewInit() != GLEW_OK) ///initializes glew
    {
        std::cout << "error with glew :c" << std::endl;
        std::cin;
        return -1;
    }

    glfwSwapInterval(1);


    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);       //enables glError + binds message callback on error




    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    // Setup Platform/Renderer bindings
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 460");
    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    ImPlot::CreateContext();


    Solver s;
    s.AddNode(3); //Adds the Number of non 0v Nodes
    
    //note to self make this suck less
    //also need to add dynamic stuffs later
    s.AddComponent<ACSource>(0, -1, [](float t) {return 2 * sin(t); });
    s.AddComponent<Resistor>(-1, 1, 1000);
    s.AddComponent<Resistor>(1, 2, 10000);
    s.AddComponent<OpAmp>(0, 1, 2);

    float T = 0;
    std::vector<float> Tx;
    std::vector<float> Vx;
    Tx.emplace_back(T);
    Tx.reserve(2000);
    Vx.reserve(2000);

    while (!glfwWindowShouldClose(window))
    {


        glfwPollEvents();
        glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);

        // feed inputs to dear imgui, start new frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // rendering our geometries
        s.Solve();
        Vx.emplace_back(s.Last(2,0));
        T += s.dt;
        Tx.emplace_back(T);
        if (ImPlot::BeginPlot("graph")) {
            ImPlot::PlotLine("v2", Tx.data(), Vx.data(), Tx.size());
            ImPlot::EndPlot();
        }
        //// render your GUI
        ImGui::Begin("Demo window");
        ImGui::Button("Hello!");
        ImGui::End();


        // Render dear imgui into screen
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glfwSwapBuffers(window);
    }

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
