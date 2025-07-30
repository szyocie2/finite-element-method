#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>


using namespace std;

#define NPC 9  //4, 9, 16
#define npc 3  //2, 3, 4

struct GlobalData {
    double SimulationTime;
    double SimulationStepTime;
    double Conductivity;
    double Alfa;
    double Tot;
    double InitialTemp;
    double Density;
    double SpecificHeat;
    int NodesNumber;
    int ElementsNumber;
};

struct Node {
    double x, y; //wspó³rzêdne wêz³a w uk³adzie globalnym
    bool BC;     //flaga okreœlaj¹ca, czy wêze³ spe³nia warunek brzegowy (true/false)
};

struct Jakobian {
    double J[2][2];   // Macierz Jakobiego
    double J1[2][2];  // Odwrócona Macierz Jakobiego
    double detJ;      // Wyznacznik Jakobianu
};

struct ElemUniv {
    double Ne[NPC][4];                  // dN/dksi dla ka¿dego punktu ca³kowania i ka¿dej funkcji kszta³tu
    double Nn[NPC][4];                  // dN/deta dla ka¿dego punktu ca³kowania i ka¿dej funkcji kszta³tu
    double Nc[NPC][4];                  // wartoœci funkcji kszta³tu, potrzebne do obliczeñ macierzy C
    double multiplied_weights[NPC];     // pomno¿one wagi
    double weights[npc];                // wagi
    double NPC_gaussPoints[NPC][2];     // Punkty Gaussa dla ka¿dego punktu ca³kowania NPC
    double npc_gaussPoints[npc * 4][2]; // Punkty Gaussa dla ka¿dego punktu ca³kowania npc krawêdzi elementu

    struct Surface {
        double N[npc][4]; // Macierz funkcji kszta³tu dla 4 œcianek, ka¿da ma 4 funkcje kszta³tu
    } surface[4];         // 4 powierzchnie dla ka¿dej œciany
};

struct Element {
    int ID[4] = {};             // identyfikatory wêz³ów elementu
    Jakobian Jacobian[NPC] = {};// dane Jakobiego dla ka¿dego punktu ca³kowania
    double dN_dx[NPC][4] = {};  // dN/dx dla ka¿dej funkcji kszta³tu w ka¿dym punkcie ca³kowania
    double dN_dy[NPC][4] = {};  // dN/dy dla ka¿dej funkcji kszta³tu w ka¿dym punkcie ca³kowania
    double H[4][4] = {};        // finalna macierz H
    double C[4][4] = {};        // finalna macierz C
    double P[4] = {};
};

struct Grid {
    Node* nodes;
    Element* elements;

    Grid(int nodesNumber, int elementsNumber) {
        nodes = new Node[nodesNumber];
        elements = new Element[elementsNumber];
    }

    ~Grid() {
        delete[] nodes;
        delete[] elements;
    }
};

// WskaŸnik na macierz globaln¹ H (dynamically allocated)
double** HG = nullptr;
double** CG = nullptr;
double* PG = nullptr;
double* t1 = nullptr;

void wypiszH(const Element&, int);
void wypiszG(int);
void gaussTemperatura(int, double);

void initializeGaussPointsAndWeights(ElemUniv& elemUniv) {

#if NPC == 4
    double NPC_gaussPoints[4][2] = {
        {-0.57735026919, -0.57735026919},
        { 0.57735026919, -0.57735026919},
        { 0.57735026919,  0.57735026919},
        {-0.57735026919,  0.57735026919}
    };
    double multiplied_weights[4] = { 1.0, 1.0, 1.0, 1.0 };

#elif NPC == 9
    double NPC_gaussPoints[9][2] = {
        {-0.77459666924, -0.77459666924},
        { 0.0,           -0.77459666924},
        { 0.77459666924, -0.77459666924},
        {-0.77459666924,  0.0},
        { 0.0,            0.0},
        { 0.77459666924,  0.0},
        {-0.77459666924,  0.77459666924},
        { 0.0,            0.77459666924},
        { 0.77459666924,  0.77459666924}
    };
    double multiplied_weights[9] = { 0.30864197531, 0.49382716049, 0.30864197531, 
                                     0.49382716049, 0.79012345679, 0.49382716049, 
                                     0.30864197531, 0.49382716049, 0.30864197531 };

#elif NPC == 16
    double g1 = -0.8611363116, g2 = -0.3399810436, g3 = 0.3399810436, g4 = 0.8611363116;
    double w1 = 0.3478548451,  w2 = 0.6521451549;
    double NPC_gaussPoints[16][2] = {
        {g1, g1}, {g2, g1}, {g3, g1}, {g4, g1},
        {g1, g2}, {g2, g2}, {g3, g2}, {g4, g2},
        {g1, g3}, {g2, g3}, {g3, g3}, {g4, g3},
        {g1, g4}, {g2, g4}, {g3, g4}, {g4, g4}
    };
    double multiplied_weights[16] = {   w1 * w1, w2 * w1, w2 * w1, w1 * w1,
                                        w1 * w2, w2 * w2, w2 * w2, w1 * w2,
                                        w1 * w2, w2 * w2, w2 * w2, w1 * w2,
                                        w1 * w1, w2 * w1, w2 * w1, w1 * w1 };

#else
#error "Invalid number of integration points. Please set NPC to 4, 9, or 16."
#endif

    // Kopiowanie do elemUniv
    for (int i = 0; i < NPC; ++i) {
        elemUniv.NPC_gaussPoints[i][0] = NPC_gaussPoints[i][0];
        elemUniv.NPC_gaussPoints[i][1] = NPC_gaussPoints[i][1];
        elemUniv.multiplied_weights[i] = multiplied_weights[i];
    }

#if npc == 2
    double npc_gaussPoints[8][2] = {
        { -0.57735026919, -1}, {  0.57735026919, -1},//dolna krawêdŸ elementu
        {  1, -0.57735026919}, {  1,  0.57735026919},//prawa
        {  0.57735026919,  1}, { -0.57735026919,  1},//górna
        { -1,  0.57735026919}, { -1, -0.57735026919},//lewa
    };
    double weights[2] = { 1.0, 1.0 };

#elif npc == 3
    double npc_gaussPoints[12][2] = {
        { -0.77459666924, -1}, {  0, -1}, {  0.77459666924, -1},//dolna krawêdŸ elementu
        {  1, -0.77459666924}, {  1,  0}, {  1,  0.77459666924},//prawa
        {  0.77459666924,  1}, {  0,  1}, { -0.77459666924,  1},//górna
        { -1,  0.77459666924}, { -1,  0}, { -1, -0.77459666924} //lewa
    };
    double weights[3] = { 0.55555555555,0.88888888888, 0.55555555555 };

#elif npc == 4
    double npc_gaussPoints[16][2] = {
        { -0.8611363116, -1}, { -0.3399810436, -1}, {  0.3399810436, -1}, {  0.8611363116, -1},//dolna krawêdŸ elementu
        {  1, -0.8611363116}, {  1, -0.3399810436}, {  1,  0.3399810436}, {  1,  0.8611363116},//prawa
        {  0.8611363116,  1}, {  0.3399810436,  1}, { -0.3399810436,  1}, { -0.8611363116,  1},//górna
        { -1,  0.8611363116}, { -1,  0.3399810436}, { -1, -0.3399810436}, { -1, -0.8611363116} //lewa
    };
    double weights[4] = { 0.347855, 0.6552145, 0.6552145, 0.347855 };

#else
#error "Invalid number of integration points. Please set npc to 2, 3, or 4."
#endif

    // Kopiowanie do elemUniv
    for (int i = 0; i < npc * 4; i++) {
        elemUniv.npc_gaussPoints[i][0] = npc_gaussPoints[i][0];
        elemUniv.npc_gaussPoints[i][1] = npc_gaussPoints[i][1];
    }
    for (int i = 0; i < npc; i++) {
        elemUniv.weights[i] = weights[i];
    }
}

void ObliczPochodneFunkcjiKsztaltu(ElemUniv& elemUniv) {

    initializeGaussPointsAndWeights(elemUniv);

    for (int ip = 0; ip < NPC; ++ip) {//obliczamy pochodne funkcji ksztaltu wzgledem ksi (Ne) oraz eta (Nn) oraz dla macierzy C
        double xi = elemUniv.NPC_gaussPoints[ip][0];
        double eta = elemUniv.NPC_gaussPoints[ip][1];

        elemUniv.Ne[ip][0] = -0.25 * (1 - eta);
        elemUniv.Ne[ip][1] = 0.25 * (1 - eta);
        elemUniv.Ne[ip][2] = 0.25 * (1 + eta);
        elemUniv.Ne[ip][3] = -0.25 * (1 + eta);

        elemUniv.Nn[ip][0] = -0.25 * (1 - xi);
        elemUniv.Nn[ip][1] = -0.25 * (1 + xi);
        elemUniv.Nn[ip][2] = 0.25 * (1 + xi);
        elemUniv.Nn[ip][3] = 0.25 * (1 - xi);

        elemUniv.Nc[ip][0] = 0.25 * (1 - xi) * (1 - eta);
        elemUniv.Nc[ip][1] = 0.25 * (1 + xi) * (1 - eta);
        elemUniv.Nc[ip][2] = 0.25 * (1 + xi) * (1 + eta);
        elemUniv.Nc[ip][3] = 0.25 * (1 - xi) * (1 + eta);
    }

    for (int s = 0; s < 4; ++s) { // 4 powierzchnie
        for (int ip = 0; ip < npc; ++ip) {//obliczamy pochodne funkcji ksztaltu wzgledem ksi (Ne) oraz eta (Nn)
            double xi = elemUniv.npc_gaussPoints[npc * s + ip][0];
            double eta = elemUniv.npc_gaussPoints[npc * s + ip][1];

            elemUniv.surface[s].N[ip][0] = 0.25 * (1 - xi) * (1 - eta); // N1
            elemUniv.surface[s].N[ip][1] = 0.25 * (1 + xi) * (1 - eta); // N2
            elemUniv.surface[s].N[ip][2] = 0.25 * (1 + xi) * (1 + eta); // N3
            elemUniv.surface[s].N[ip][3] = 0.25 * (1 - xi) * (1 + eta); // N4
        }
    }
}

void AgregujDoMacierzyGlobalnej(Grid& grid, int elemIndex) {
    Element& element = grid.elements[elemIndex];

    // Agregacja lokalnej macierzy H do macierzy globalnej HG
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int global_i = element.ID[i] - 1;
            int global_j = element.ID[j] - 1;

            HG[global_i][global_j] += element.H[i][j];
        }
    }

    //Agregacja lokalnej macierzy C do macierzy globalnej CG
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int global_i = element.ID[i] - 1;
            int global_j = element.ID[j] - 1;

            CG[global_i][global_j] += element.C[i][j];
        }
    }

    // Agregacja lokalnego wektoru P do wektoru globalnego PG
    for (int i = 0; i < 4; ++i) {
        int global_i = element.ID[i] - 1;

        PG[global_i] += element.P[i];
    }
}

void ObliczJakobian(Grid& grid, int elemIndex, const ElemUniv& elemUniv, double conductivity, double Density, double SpecificHeat) {

    Element& element = grid.elements[elemIndex];
    for (int ip = 0; ip < NPC; ++ip) {
        double dXdXi = 0.0, dXdEta = 0.0, dYdXi = 0.0, dYdEta = 0.0;

        for (int k = 0; k < 4; ++k) {//wyznaczamy jakobian przekszta³cenia
            int nodeIndex = element.ID[k] - 1; // -1 zeby indeks by³ od zera
            Node& node = grid.nodes[nodeIndex];

            dXdXi += elemUniv.Ne[ip][k] * node.x;
            dXdEta += elemUniv.Nn[ip][k] * node.x;
            dYdXi += elemUniv.Ne[ip][k] * node.y;
            dYdEta += elemUniv.Nn[ip][k] * node.y;
        }

        //otrzymany jakobian dla danego punktu ca³kowania
        element.Jacobian[ip].J[0][0] = dXdXi;
        element.Jacobian[ip].J[0][1] = dXdEta;
        element.Jacobian[ip].J[1][0] = dYdXi;
        element.Jacobian[ip].J[1][1] = dYdEta;

        //licze wyznaczik detJ
        element.Jacobian[ip].detJ = dXdXi * dYdEta - dXdEta * dYdXi;

        if (element.Jacobian[ip].detJ != 0) {
            double invDetJ = 1.0 / element.Jacobian[ip].detJ; //odwrócony wyznacznik jakobiego

            //odwracam Jakobian
            element.Jacobian[ip].J1[0][0] = dYdEta * invDetJ;
            element.Jacobian[ip].J1[0][1] = -dXdEta * invDetJ;
            element.Jacobian[ip].J1[1][0] = -dYdXi * invDetJ;
            element.Jacobian[ip].J1[1][1] = dXdXi * invDetJ;

            double H_pc[4][4] = {};  // Tymczasowa macierz H dla tego punktu integracji

            for (int i = 0; i < 4; ++i) {// wyliczanie pochodnych dN/dx oraz dN/dy
                element.dN_dx[ip][i] = element.Jacobian[ip].J1[0][0] * elemUniv.Ne[ip][i] + element.Jacobian[ip].J1[1][0] * elemUniv.Nn[ip][i];
                element.dN_dy[ip][i] = element.Jacobian[ip].J1[0][1] * elemUniv.Ne[ip][i] + element.Jacobian[ip].J1[1][1] * elemUniv.Nn[ip][i];
            }

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    //wartosci iloczynów pochodnych dN/dx oraz dN/dy w danych wspo³rzêdnych
                    double dNdx_i_dNdx_j = element.dN_dx[ip][i] * element.dN_dx[ip][j];
                    double dNdy_i_dNdy_j = element.dN_dy[ip][i] * element.dN_dy[ip][j];

                    // Iloczyn sumy iloczynów pochodnych funkcji kszta³tu z przewodnoœci¹, wagami i wyznacznikiem Jakobianu
                    H_pc[i][j] = conductivity * (dNdx_i_dNdx_j + dNdy_i_dNdy_j) * element.Jacobian[ip].detJ * elemUniv.multiplied_weights[ip];

                    // Sumuj do finalnej macierzy H
                    element.H[i][j] += H_pc[i][j];
                }
            }

            //OBLICZANIE C
            double Hclocal[4][4] = {};

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    Hclocal[i][j] += Density * SpecificHeat * elemUniv.Nc[ip][i] * elemUniv.Nc[ip][j] * element.Jacobian[ip].detJ * elemUniv.multiplied_weights[ip];
                    element.C[i][j] += Hclocal[i][j];
                }
            }
            //cout << "\nSprawdzenie C:" << ip + 1 << "\n";
            //for (int i = 0; i < 4; ++i) {
            //    for (int j = 0; j < 4; ++j) {
            //        cout << setw(10) << Hclocal[i][j] << " ";
            //    }
            //    cout << "\n";
            //}
            //cout << "============================================\n";
        }
        else {
            cout << "  Inverse Jacobian does not exist (detJ = 0)." << endl;
        }
    }
}

void ObliczBC(Grid& grid, int elemIndex, const ElemUniv& elemUniv, double alfa, double tot) {
    Element& element = grid.elements[elemIndex];

    //wypiszH(grid.elements[elemIndex], elemIndex);
    for (int s = 0; s < 4; ++s) { // Iteracja po krawêdziach
        double Hbclocal[4][4] = {};

        int node1 = element.ID[s] - 1;               // Pierwszy wêze³ na krawêdzi
        int node2 = element.ID[(s + 1) % 4] - 1;     // Drugi wêze³ na krawêdzi

        if (grid.nodes[node1].BC && grid.nodes[node2].BC) { // Sprawdzenie BC
            double dx = grid.nodes[node2].x - grid.nodes[node1].x;
            double dy = grid.nodes[node2].y - grid.nodes[node1].y;
            double detJ_surface = sqrt(dx * dx + dy * dy) / 2.0; //D³ugoœæ krawêdzi / 2 - d³. znormalizowanego uk³adu

            for (int ip = 0; ip < npc; ip++) {
                // Wartoœci funkcji kszta³tu na krawêdzi
                double N[4] = {
                    elemUniv.surface[s].N[ip][0],
                    elemUniv.surface[s].N[ip][1],
                    elemUniv.surface[s].N[ip][2],
                    elemUniv.surface[s].N[ip][3]
                };

                //Iteracja po funkcjach kszta³tu (lokalne wêz³y)
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        Hbclocal[i][j] += alfa * N[i] * N[j] * detJ_surface * elemUniv.weights[ip];
                    }
                    element.P[i] += alfa * N[i] * tot * detJ_surface * elemUniv.weights[ip];
                }
            }

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    element.H[i][j] += Hbclocal[i][j];
                }
            }
        }
        //cout << "\nSprawdzenie Hbc:" << s + 1 << "\n";
        //for (int i = 0; i < 4; ++i) {
        //    for (int j = 0; j < 4; ++j) {
        //        cout << setw(10) << Hbclocal[i][j] << " ";
        //    }
        //    cout << "\n";
        //}
        //cout << "============================================\n";
    }
    //for (int i = 0; i < 4; ++i) {
    //    cout << element.P[i] << " ";
    //}
    //cout << endl;
}

int main() {

    //odczyt z pliku
    ifstream file("Test2_4_4_MixGrid.txt"); //Test1_4_4 Test2_4_4_MixGrid Test3_31_31_kwadrat
    if (!file) {
        cout << "Nie mo¿na otworzyæ pliku!" << endl;
        return 1;
    }

    GlobalData globalData;
    string line;

    file >> line >> globalData.SimulationTime;
    file >> line >> globalData.SimulationStepTime;
    file >> line >> globalData.Conductivity;
    file >> line >> globalData.Alfa;
    file >> line >> globalData.Tot;
    file >> line >> globalData.InitialTemp;
    file >> line >> globalData.Density;
    file >> line >> globalData.SpecificHeat;
    file >> line >> line >> globalData.NodesNumber;
    file >> line >> line >> globalData.ElementsNumber;

    Grid grid(globalData.NodesNumber, globalData.ElementsNumber);

    while (getline(file, line)) {
        if (line.find("*Node") != string::npos) {
            for (int i = 0; i < globalData.NodesNumber; i++) {
                int nodeId;
                char comma;
                file >> nodeId >> comma >> grid.nodes[i].x >> comma >> grid.nodes[i].y;
                grid.nodes[i].BC = false; // Domyœlnie ustawiamy BC na false
            }
            break;
        }
    }

    while (getline(file, line)) {
        if (line.find("*Element") != string::npos) {
            for (int i = 0; i < globalData.ElementsNumber; ++i) {
                int elementId;
                char comma;
                file >> elementId >> comma >> grid.elements[i].ID[0] >> comma >> grid.elements[i].ID[1] >> comma >> grid.elements[i].ID[2] >> comma >> grid.elements[i].ID[3];
            }
            break;
        }
    }

    while (getline(file, line)) {
        if (line.find("*BC") != string::npos) {
            getline(file, line); // Wczytanie linii zawieraj¹cej BC
            stringstream ss(line);
            int nodeId;

            // Przetwarzanie wêz³ów z warunkami brzegowymi
            while (ss >> nodeId) {
                if (nodeId > 0 && nodeId <= globalData.NodesNumber) {
                    grid.nodes[nodeId - 1].BC = true; // Ustawienie BC na true dla wêz³a
                }
                if (ss.peek() == ',') ss.ignore(); // Ignorowanie przecinków
            }
            break;
        }
    }

    file.close();
    //koniec odczytu z pliku

    // Debug: Wyœwietlenie stanu wêz³ów
    //for (int i = 0; i < globalData.NodesNumber; i++) {
    //    cout << "Node " << i + 1 << ", BC=" << grid.nodes[i].BC << endl;
    //}

    // Dynamiczna alokacja macierzy globalnej HG
    HG = new double* [globalData.NodesNumber];
    for (int i = 0; i < globalData.NodesNumber; ++i) {
        HG[i] = new double[globalData.NodesNumber] {};
    }

    // Dynamiczna alokacja macierzy globalnej CG
    CG = new double* [globalData.NodesNumber];
    for (int i = 0; i < globalData.NodesNumber; ++i) {
        CG[i] = new double[globalData.NodesNumber] {};
    }

    // Dynamiczna alokacja macierzy globalnej BCG
    PG = new double [globalData.NodesNumber] {};
    t1 = new double [globalData.NodesNumber] {};

    ElemUniv elemUniv;
    ObliczPochodneFunkcjiKsztaltu(elemUniv);

    double deltaT = globalData.SimulationTime / globalData.SimulationStepTime;

    double* t0 = new double[globalData.NodesNumber];

    // Inicjalizacja temperatury pocz¹tkowej
    for (int i = 0; i < globalData.NodesNumber; i++) {
        t0[i] = globalData.InitialTemp;
    }

    for (int i = 0; i < deltaT; i++) {
        // Resetowanie macierzy H, C i P dla ka¿dego elementu
        for (int i = 0; i < globalData.ElementsNumber; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    grid.elements[i].H[j][k] = 0.0;
                    grid.elements[i].C[j][k] = 0.0;
                }
                grid.elements[i].P[j] = 0.0;
            }
        }

        // Resetowanie macierzy globalnych HG, CG oraz wektora PG
        for (int i = 0; i < globalData.NodesNumber; ++i) {
            for (int j = 0; j < globalData.NodesNumber; ++j) {
                HG[i][j] = 0.0;
                CG[i][j] = 0.0;
            }
            PG[i] = 0.0;
        }

        for (int i = 0; i < globalData.ElementsNumber; i++) {
            ObliczJakobian(grid, i, elemUniv, globalData.Conductivity, globalData.Density, globalData.SpecificHeat);
            ObliczBC(grid, i, elemUniv, globalData.Alfa, globalData.Tot);
            //wypiszH(grid.elements[i], i);
            AgregujDoMacierzyGlobalnej(grid, i);
        }

        for (int i = 0; i < globalData.NodesNumber; ++i) {
            for (int j = 0; j < globalData.NodesNumber; ++j) {
                CG[i][j] = CG[i][j] / globalData.SimulationStepTime;
            }
        }

        for (int i = 0; i < globalData.NodesNumber; ++i) {
            for (int j = 0; j < globalData.NodesNumber; ++j) {
                HG[i][j] += CG[i][j];
            }
        }

        for (int i = 0; i < globalData.NodesNumber; i++) {
            for (int j = 0; j < globalData.NodesNumber; j++) {
                PG[i] += CG[i][j] * t0[j];
            }
        }
        gaussTemperatura(globalData.NodesNumber, i);

        for (int i = 0; i < globalData.NodesNumber; i++) {
            t0[i] = t1[i];
        }
    }
    cout << endl;
    //wypiszG(globalData.NodesNumber);  // Wyœwietlanie macierzy globalnych

    // Zwolnienie pamiêci
    for (int i = 0; i < globalData.NodesNumber; ++i) {
        delete[] HG[i];
    }
    delete[] HG;

    for (int i = 0; i < globalData.NodesNumber; ++i) {
        delete[] CG[i];
    }
    delete[] CG;

    delete[] PG;
    delete[] t1;
    delete[] t0;

    return 0;
}

void gaussTemperatura(int NodesNumber, double numerTemp) {

    double** HGreversed = nullptr;
    HGreversed = new double* [NodesNumber];
    for (int i = 0; i < NodesNumber; ++i) {
        HGreversed[i] = new double[NodesNumber + 1] {};
    }

    for (int i = 0; i < NodesNumber; ++i) {
        for (int j = 0; j < NodesNumber + 1; ++j) {
            if (j < NodesNumber) {
                HGreversed[i][j] = HG[i][j];
            }
            else {
                HGreversed[i][j] = PG[i];
            }
        }
    }

    //sprowadzanie macierzy do postaci trójk¹tnej górnej
    for (int k = 0; k < NodesNumber; k++) {
        for (int i = k + 1; i < NodesNumber; i++) {
            for (int j = NodesNumber; j >= 0; j--) {
                HGreversed[i][j] -= HGreversed[k][j] * (HGreversed[i][k] / HGreversed[k][k]);
            }
        }
    }

    //postêpowanie wsteczne
    for (int i = NodesNumber - 1; i >= 0; i--) {
        t1[i] = HGreversed[i][NodesNumber];
        for (int j = i + 1; j < NodesNumber; j++) {
            t1[i] -= HGreversed[i][j] * t1[j];
        }
        t1[i] /= HGreversed[i][i];
    }

    //cout << defaultfloat << "\n\nObliczona temperatura t " << numerTemp + 1 <<  " dla elementow metoda Gaussa:\n";

    double max = 0;
    double min = 1200;

    for (int i = 0; i < NodesNumber; i++) {
        if (t1[i] > max)max = t1[i];
        if (t1[i] < min)min = t1[i];
    }

    cout << fixed << setprecision(10); // Dok³adnoœæ do 3 miejsc
    cout << endl << min << "\t" << max;


    //cout << "\nMacierz zamieniona pod Gaussa:\n";
    //for (int i = 0; i < NodesNumber; ++i) {
    //    for (int j = 0; j < NodesNumber + 1; ++j) {
    //        cout << setw(10) << HGreversed[i][j] << " ";
    //    }
    //    cout << "\n";
    //}
}

void wypiszH(const Element& element, int i) {

    cout << "Macierz H dla " << i + 1 << " elementu:\n";
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout << setw(10) << element.H[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "BC ";
    for (int i = 0; i < 4; ++i) {
        cout << element.P[i] << " ";
    }
    cout << "\n============================================\n";
}

void wypiszG(int NodesNumber) {
    cout << "\nMacierz globalna H:\n";
    for (int i = 0; i < NodesNumber; ++i) {
        for (int j = 0; j < NodesNumber; ++j) {
            cout << setw(10) << HG[i][j] << " ";
        }
        cout << "\n";
    }

    cout << "BC ";//wektor globalny P
    for (int i = 0; i < NodesNumber; ++i) {
        cout << setw(10) << PG[i] << " ";
    }

    cout << "\nMacierz globalna C:\n";
    for (int i = 0; i < NodesNumber; ++i) {
        for (int j = 0; j < NodesNumber; ++j) {
            cout << setw(10) << CG[i][j] << " ";
        }
        cout << "\n";
    }
}