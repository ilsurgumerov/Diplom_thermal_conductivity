//#include <iostream>
//#include <fstream>
//#include <string>
//#include <cmath>
//#include <Eigen/Dense>
//
//using namespace std;
//using namespace Eigen;
//
//typedef double T;
//
//T M = 7e3;
//
//struct point {
//    T x, y, z;
//
//    point(T X = 0, T Y = 0, T Z = 0) {
//        x = X;
//        y = Y;
//        z = Z;
//    }
//};
//
//point operator*(const T k, point p) {
//    p.x *= k;
//    p.y *= k;
//    p.z *= k;
//    return p;
//}
//
//point operator+(point p1, point p2) {
//    p1.x += p2.x;
//    p1.y += p2.y;
//    p1.z += p2.z;
//    return p1;
//}
//
//point operator-(point p1, point p2) {
//    p1.x -= p2.x;
//    p1.y -= p2.y;
//    p1.z -= p2.z;
//    return p1;
//}
//
//ostream& operator<< (ostream& out, const point& p)
//{
//    out << p.x << " " << p.y << " " << p.z;
//    return out;
//}
//
//istream& operator>> (istream& in, point& p)
//{
//    in >> p.x;
//    in >> p.y;
//    in >> p.z;
//    return in;
//}
//
//struct elem {
//    int* nodesNum;
//    int num;
//};
//
//// 2 or 4
//int ReadNum(const char* fileName)
//{
//    int num;
//    ifstream ifile;
//    ifile.open(fileName);
//    if (!ifile.is_open())
//    {
//        cout << "Error! Could not open file!" << endl;
//        return -1;
//    }
//    ifile >> num;
//    ifile.close();
//    return num;
//}
//
//int ReadRequestFile(const char* fileName, int num, point* coords, int* elemNum, int& type) {
//
//    ifstream ifile;
//    ifile.open(fileName);
//    if (!ifile.is_open())
//    {
//        cout << "Error! Could not open file!" << endl;
//        return -1;
//    }
//    ifile >> num;       //2 or 4
//    for (int i = 0; i < num; ++i)
//    {
//        ifile >> coords[i].x;
//        ifile >> coords[i].y;
//        ifile >> coords[i].z;
//    }
//
//    if (num == 2)
//        ifile >> elemNum[0];
//    else
//    {
//        ifile >> elemNum[0];
//        ifile >> elemNum[1];
//    }
//
//    ifile >> type;
//    ifile.close();
//    return 0;
//}
//
//void Mesh2p(const point coord0, const point coord1, int np, point* nodes, int i0) {
//    for (int i = 0; i < np; ++i)
//    {
//        T xi = -1 + 2. * i / (np - 1);//local coords
//
//        T N1 = (1 - xi) / 2.;
//        T N2 = (1 + xi) / 2.;
//        nodes[i + i0] = N1 * coord0 + N2 * coord1;
//
//    }
//}
//
//void MakeMesh2p(const int elemNum, const point* coords, const int type, int& NE, int& NP, int& NC, point*& nodes, elem*& elems, int**& CPN, int*& CPNsizes)
//{
//
//    NE = elemNum;
//    if (type == 1) {
//        NP = elemNum + 1;
//    }
//    else {
//        NP = 2 * elemNum + 1;
//    }
//
//    NC = 2;
//    nodes = new point[NP];
//
//
//    Mesh2p(coords[0], coords[1], NP, nodes, 0);
//
//
//
//    elems = new elem[NE];
//    if (type == 1) {
//        for (int i = 0; i < NE; ++i)
//        {
//            elems[i].num = 2;
//            elems[i].nodesNum = new int[elems[i].num];
//            elems[i].nodesNum[0] = i + 1;
//            elems[i].nodesNum[1] = i + 2;
//        }
//    }
//    else {
//        for (int i = 0; i < NE; ++i)
//        {
//            elems[i].num = 3;
//            elems[i].nodesNum = new int[elems[i].num];
//            elems[i].nodesNum[0] = 2 * i + 1;
//            elems[i].nodesNum[1] = 2 * i + 2;
//            elems[i].nodesNum[2] = 2 * i + 3;
//        }
//    }
//
//    CPNsizes = new int[2];
//    CPNsizes[0] = 1; CPNsizes[1] = 1;
//
//    CPN = new int* [2];
//    CPN[0] = new int[1]; CPN[0][0] = 1;
//    CPN[1] = new int[1]; CPN[1][0] = NP;
//}
//
//void Mesh4pType1(const point* coords, const int np0, const int np1, point* nodes) {
//
//    point a0 = coords[3] - coords[0];
//    point a1 = coords[2] - coords[1];
//    for (int i = 0; i < np0; ++i) {
//        T Li = i * 1.0 / (np0 - 1);
//        point b0 = coords[0] + Li * a0;
//        point b1 = coords[1] + Li * a1;
//        Mesh2p(b0, b1, np1, nodes, i * np1);
//    }
//}
//
//void MeshCorrToType4(const point* coords, const int np0, const int np1, point* nodes) {
//
//    point a0 = coords[3] - coords[0];
//    point a1 = coords[2] - coords[1];
//    int k = 0;
//    for (int i = 0; i < np0 - 1; ++i) {
//        T Li = (i + 0.5) * 1.0 / (np0 - 1);
//        point b0 = coords[0] + Li * a0;
//        point b1 = coords[1] + Li * a1;
//        for (int j = 0; j < np1 - 1; ++j) {
//            T Lj = (j + 0.5) * 1.0 / (np1 - 1);
//            nodes[np0 * np1 + k] = b0 + Lj * (b1 - b0);
//            k++;
//        }
//    }
//}
//
//void MakeMesh4p(const int* elemNum, const point* coords, const int type, int& NE, int& NP, int& NC, point*& nodes, elem*& elems, int**& CPN, int*& CPNsizes, Matrix<int, Dynamic, Dynamic>& bcInds) {
//
//
//    const int np0 = elemNum[0] + 1;
//    const int np1 = elemNum[1] + 1;
//    if (type == 1) {
//        NE = elemNum[0] * elemNum[1];
//    }
//    if (type == 2 || type == 3) {
//        NE = 2 * elemNum[0] * elemNum[1];
//    }
//    if (type == 4) {
//        NE = 4 * elemNum[0] * elemNum[1];
//    }
//
//    if (type != 4) {
//        NP = np0 * np1;
//    }
//    else {
//        NP = np0 * np1 + elemNum[0] * elemNum[1];
//    }
//
//    nodes = new point[NP];
//
//    Mesh4pType1(coords, np0, np1, nodes);
//
//
//
//    NC = 1;
//    CPNsizes = new int[1];
//    CPNsizes[0] = 2 * (elemNum[0] + elemNum[1]);
//    CPN = new int* [1];
//    CPN[0] = new int[CPNsizes[0]];
//
//    bcInds = MatrixXi::Zero(4, CPNsizes[0]);
//
//    for (int i = 0; i < np1; ++i) {
//        CPN[0][i] = i + 1;
//        bcInds(0, i) = i + 1;
//    }
//    bcInds(1, 0) = np1;
//
//    for (int i = 1; i < np0; ++i) {
//        CPN[0][i + np1 - 1] = i * np1 + np1;
//        bcInds(1, i) = i * np1 + np1;
//    }
//    bcInds(2, 0) = (np0 - 1) * np1 + np1;
//
//    for (int i = 1; i < np1; ++i) {
//        CPN[0][i + np1 + np0 - 2] = np0 * np1 - i;
//        bcInds(2, i) = np0 * np1 - i;
//    }
//    bcInds(3, 0) = np0 * np1 - np1 + 1;
//
//    for (int i = 1; i < np0 - 1; ++i) {
//        CPN[0][i + 2 * np1 + np0 - 3] = (np0 - i - 1) * np1 + 1;
//        bcInds(3, i) = (np0 - i - 1) * np1 + 1;
//    }
//    bcInds(3, np0 - 1) = 1;
//
//    elems = new elem[NE];
//
//    auto node1{
//        [](int i, int j, int np1) {
//            return i * np1 + j + 1;
//        }
//    };
//    auto node2{
//        [](int i, int j, int np1) {
//            return i * np1 + j + 2;
//        }
//    };
//    auto node3{
//        [](int i, int j, int np1) {
//            return (i + 1) * np1 + j + 2;
//        }
//    };
//    auto node4{
//        [](int i, int j, int np1) {
//            return (i + 1) * np1 + j + 1;
//        }
//    };
//    auto nodec{
//        [](int i0, int k) {
//            return i0 + k;
//        }
//    };
//
//    if (type == 1) {
//        for (int i = 0; i < elemNum[0]; ++i)
//        {
//            for (int j = 0; j < elemNum[1]; ++j)
//            {
//                int ind = i * elemNum[1] + j;
//                elems[ind].num = 4;
//                elems[ind].nodesNum = new int[4];
//                elems[ind].nodesNum[0] = node1(i, j, np1);
//                elems[ind].nodesNum[1] = node2(i, j, np1);
//                elems[ind].nodesNum[2] = node3(i, j, np1);
//                elems[ind].nodesNum[3] = node4(i, j, np1);
//
//            }
//
//        }
//    }
//
//
//    if (type == 2) {
//        for (int i = 0; i < elemNum[0]; ++i)
//        {
//            for (int j = 0; j < elemNum[1]; ++j)
//            {
//                int ind = i * 2 * elemNum[1] + 2 * j;
//                elems[ind].num = 3;
//                elems[ind].nodesNum = new int[3];
//                elems[ind].nodesNum[0] = node1(i, j, np1);
//                elems[ind].nodesNum[1] = node2(i, j, np1);
//                elems[ind].nodesNum[2] = node3(i, j, np1);
//
//                elems[ind + 1].num = 3;
//                elems[ind + 1].nodesNum = new int[3];
//                elems[ind + 1].nodesNum[0] = node3(i, j, np1);
//                elems[ind + 1].nodesNum[1] = node4(i, j, np1);
//                elems[ind + 1].nodesNum[2] = node1(i, j, np1);
//
//            }
//
//        }
//    }
//
//    if (type == 3) {
//
//        for (int i = 0; i < elemNum[0]; ++i)
//        {
//            for (int j = 0; j < elemNum[1]; ++j)
//            {
//                int ind = i * 2 * elemNum[1] + 2 * j;
//                elems[ind].num = 3;
//                elems[ind].nodesNum = new int[3];
//                elems[ind].nodesNum[0] = node1(i, j, np1);
//                elems[ind].nodesNum[1] = node2(i, j, np1);
//                elems[ind].nodesNum[2] = node4(i, j, np1);
//
//                elems[ind + 1].num = 3;
//                elems[ind + 1].nodesNum = new int[3];
//                elems[ind + 1].nodesNum[0] = node2(i, j, np1);
//                elems[ind + 1].nodesNum[1] = node3(i, j, np1);
//                elems[ind + 1].nodesNum[2] = node4(i, j, np1);
//
//            }
//
//        }
//    }
//
//    if (type == 4) {
//
//        MeshCorrToType4(coords, np0, np1, nodes);
//        int k = 1;
//        for (int i = 0; i < elemNum[0]; ++i)
//        {
//            for (int j = 0; j < elemNum[1]; ++j)
//            {
//                int ind = i * 4 * elemNum[1] + 4 * j;
//                elems[ind].num = 3;
//                elems[ind].nodesNum = new int[3];
//                elems[ind].nodesNum[0] = node1(i, j, np1);
//                elems[ind].nodesNum[1] = node2(i, j, np1);
//                elems[ind].nodesNum[2] = nodec(np0 * np1, k);
//
//                elems[ind + 1].num = 3;
//                elems[ind + 1].nodesNum = new int[3];
//                elems[ind + 1].nodesNum[0] = node2(i, j, np1);
//                elems[ind + 1].nodesNum[1] = node3(i, j, np1);
//                elems[ind + 1].nodesNum[2] = nodec(np0 * np1, k);
//
//                elems[ind + 2].num = 3;
//                elems[ind + 2].nodesNum = new int[3];
//                elems[ind + 2].nodesNum[0] = node3(i, j, np1);
//                elems[ind + 2].nodesNum[1] = node4(i, j, np1);
//                elems[ind + 2].nodesNum[2] = nodec(np0 * np1, k);
//
//                elems[ind + 3].num = 3;
//                elems[ind + 3].nodesNum = new int[3];
//                elems[ind + 3].nodesNum[0] = node4(i, j, np1);
//                elems[ind + 3].nodesNum[1] = node1(i, j, np1);
//                elems[ind + 3].nodesNum[2] = nodec(np0 * np1, k);
//                k++;
//            }
//
//        }
//
//    }
//
//}
//
//int WriteMeshFile(const char* name, const int NE, const int NP, const int NC, point* nodes, elem* elems, int** CPN, int* CPNsizes) {
//    ofstream ofile;
//    ofile.open(name);
//    if (!ofile.is_open())
//    {
//        cout << "Error! Could not open file!" << endl;
//        return -1;
//    }
//
//    ofile << NE << " " << NP << " " << NC << endl;
//    for (int i = 0; i < NE; ++i)
//    {
//        ofile << i + 1 << " " << elems[i].num << " ";
//        for (int j = 0; j < elems[i].num; ++j) {
//            ofile << elems[i].nodesNum[j] << " ";
//        }
//        ofile << endl;
//    }
//
//    for (int i = 0; i < NP; ++i)
//    {
//        ofile << i + 1 << " " << nodes[i] << endl;
//    }
//
//
//
//    for (int i = 0; i < NC; ++i) {
//        ofile << CPNsizes[i] << " ";
//    }
//    ofile << endl;
//    for (int i = 0; i < NC; ++i) {
//        for (int j = 0; j < CPNsizes[i]; ++j) {
//            ofile << CPN[i][j] << endl;
//        }
//    }
//
//    ofile.close();
//    return 0;
//}
//
//void CreateMeshFile(const char* RequestFileName, const char* MeshFileName, Matrix<int, Dynamic, Dynamic>& bcInds) {
//
//    const int num = ReadNum(RequestFileName);
//    int type;                           //type = 1
//    point* coords = new point[num];
//    int elemNum[2];
//
//    ReadRequestFile(RequestFileName, num, coords, elemNum, type);
//
//    int NE, NP, NC;
//    point* nodes = nullptr;
//    elem* elems = nullptr;
//    int** CPN = nullptr;
//    int* CPNsizes = nullptr;
//    //MakeMesh2p(elemNum[0], coords, type, NE, NP, NC, nodes, elems, CPN, CPNsizes);
//    MakeMesh4p(elemNum, coords, type, NE, NP, NC, nodes, elems, CPN, CPNsizes, bcInds);
//    WriteMeshFile(MeshFileName, NE, NP, NC, nodes, elems, CPN, CPNsizes);
//
//    delete[] coords;
//    delete[] nodes;
//    for (int i = 0; i < NE; ++i) {
//        delete[] elems[i].nodesNum;
//    }
//    delete[] elems;
//    for (int i = 0; i < NC; ++i) {
//        delete[] CPN[i];
//    }
//    delete[] CPN;
//    delete[] CPNsizes;
//}
//
//int ReadMeshFile(const char* fileName, int& NE, int& NP, int& NC, point*& nodes, elem*& elems, int**& CPN, int*& CPNsizes) {
//    ifstream ifile;
//    ifile.open(fileName);
//    if (!ifile.is_open())
//    {
//        cout << "Error! Could not open file!" << endl;
//        return -1;
//    }
//
//    ifile >> NE;
//    ifile >> NP;
//    ifile >> NC;
//
//    nodes = new point[NP];
//    elems = new elem[NE];
//
//    for (int i = 0; i < NE; ++i)
//    {
//        int elemInd;
//        ifile >> elemInd;
//        elemInd--;
//        ifile >> elems[elemInd].num;
//        elems[elemInd].nodesNum = new int[elems[elemInd].num];
//        for (int j = 0; j < elems[elemInd].num; ++j) {
//            ifile >> elems[elemInd].nodesNum[j];
//        }
//    }
//
//    for (int i = 0; i < NP; ++i)
//    {
//        int nodeInd;
//        ifile >> nodeInd;
//        nodeInd--;
//        ifile >> nodes[nodeInd];
//    }
//
//    CPNsizes = new int[NC];
//    for (int i = 0; i < NC; ++i) {
//
//        ifile >> CPNsizes[i];
//    }
//
//    CPN = new int* [NC];
//    for (int i = 0; i < NC; ++i) {
//        CPN[i] = new int[CPNsizes[i]];
//        for (int j = 0; j < CPNsizes[i]; ++j) {
//            ifile >> CPN[i][j];
//        }
//    }
//
//    ifile.close();
//    return 0;
//}
//
//
//struct BoundaryCondition {
//    int nodeInd;
//    int bcType; //1 or 2
//    T value; // u0 or q0
//
//};
//
//int MakeBoundaryConditions(BoundaryCondition* bc, int bcSize, Matrix<int, Dynamic, Dynamic> bcInds, point* nodes, int side, T(*f)(T), int type) {
//
//    int count = 0;
//    while (bc[count].bcType != 0 && count < bcSize) count++;
//    for (int j = 0; j < bcSize; ++j) {
//        if (bcInds(side - 1, j) != 0) {
//            bc[count].bcType = type;
//            int ind = bc[count].nodeInd = bcInds(side - 1, j) - 1;
//            if (side - 1 == 0 || side - 1 == 2) bc[count].value = f(nodes[ind].x);
//            else bc[count].value = f(nodes[ind].y);
//            count++;
//        }
//        else break;
//    }
//    return 0;
//}
//
//int WriteBoundaryFile(const char* fileName, BoundaryCondition* bc, int bcSize) {
//    ofstream ofile;
//    ofile.open(fileName);
//    if (!ofile.is_open())
//    {
//        cout << "Error! Could not open file!" << endl;
//        return -1;
//    }
//    int count = 0;
//    while (bc[count].bcType != 0 && count < bcSize) count++;
//    ofile << count << endl;
//    for (int i = 0; i < count; ++i) {
//        ofile << bc[i].nodeInd + 1 << " ";
//        ofile << bc[i].bcType << " ";//1 2
//        ofile << bc[i].value << endl;
//    }
//
//    ofile.close();
//    return 0;
//}
//
//int ReadBoundaryFile(const char* fileName, BoundaryCondition* bc, int& bcSize) {
//    ifstream ifile;
//    ifile.open(fileName);
//    if (!ifile.is_open())
//    {
//        cout << "Error! Could not open file!" << endl;
//        return -1;
//    }
//    ifile >> bcSize;
//    for (int i = 0; i < bcSize; ++i) {
//        ifile >> bc[i].nodeInd; bc[i].nodeInd--;
//        ifile >> bc[i].bcType;//1 2
//        ifile >> bc[i].value;
//    }
//    ifile.close();
//    return 0;
//}
//
//T lambda(T x, T y) {
//    return 1;
//}
//
//T q(T x, T y) {
//    return +2 * sin(x + y);
//}
//
//void Tabulate(T(*f)(T, T), point* nodes, const int NP, T* fArray) {
//    for (int i = 0; i < NP; ++i) {
//        fArray[i] = f(nodes[i].x, nodes[i].y);
//    }
//}
//
//
//T f1(T x) {
//    return -cos(x);//2+sin(x);
//}
//
//T f2(T y) {
//    return 2 + sin(1 + y);// 2;
//}
//
//T f3(T x) {
//    return cos(1 + x);// 2 + sin(x + 1);
//}
//
//T f4(T y) {
//    return 2 + sin(y);//5 + 1+ y * y;
//}
//
//void Add(elem el, Matrix<double, 3, 3> localMatrix, Vector<double, 3> localRHS, Matrix<double, Dynamic, Dynamic>& systemMatrix, Vector<double, Dynamic>& systemRHS) {
//
//    int i = el.nodesNum[0] - 1;
//    int j = el.nodesNum[1] - 1;
//    int k = el.nodesNum[2] - 1;
//
//    systemMatrix(i, i) += localMatrix(0, 0);
//    systemMatrix(i, j) += localMatrix(0, 1);
//    systemMatrix(j, i) += localMatrix(1, 0);
//    systemMatrix(j, j) += localMatrix(1, 1);
//    systemMatrix(i, k) += localMatrix(0, 2);
//    systemMatrix(j, k) += localMatrix(1, 2);
//    systemMatrix(k, i) += localMatrix(2, 0);
//    systemMatrix(k, j) += localMatrix(2, 1);
//    systemMatrix(k, k) += localMatrix(2, 2);
//
//    systemRHS(i) += localRHS(0);
//    systemRHS(j) += localRHS(1);
//    systemRHS(k) += localRHS(2);
//
//}
//
//
//void BC1(int nodeInd, T u0, Matrix<double, Dynamic, Dynamic>& systemMatrix, Vector<double, Dynamic>& systemRHS) {
//    systemMatrix(nodeInd, nodeInd) += M;
//    systemRHS(nodeInd) += u0 * M;
//}
//
//void BC2(int nodeInd1, int nodeInd2, T q0, T len, Matrix<double, Dynamic, Dynamic>& systemMatrix, Vector<double, Dynamic>& systemRHS) {
//    systemRHS(nodeInd1) += q0 * len / 2.0;
//    systemRHS(nodeInd2) += q0 * len / 2.0;
//}
//
//
//T norm(point A, point B) {
//    return sqrt((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y));
//}
//
//void ApplyBoundaryValues(BoundaryCondition* bc, int bcSize, point* nodes, T* lambdaArray, Matrix<double, Dynamic, Dynamic>& systemMatrix, Vector<double, Dynamic>& systemRHS) {
//
//    for (int i = 0; i < bcSize - 1; ++i) {
//        int nodeInd1 = bc[i].nodeInd;
//        int nodeInd2 = bc[i + 1].nodeInd;
//        switch (bc[i].bcType) {
//        case 1:
//            BC1(nodeInd1, bc[i].value, systemMatrix, systemRHS);
//            break;
//        case 2:
//            if (bc[i + 1].bcType == 2) {
//                T len = norm(nodes[nodeInd1 - 1], nodes[nodeInd2 - 1]);
//                BC2(nodeInd1, nodeInd2, bc[i].value, len, systemMatrix, systemRHS);
//            }
//            break;
//        }
//    }
//    if (bc[bcSize - 1].bcType == 1) BC1(bc[bcSize - 1].nodeInd, bc[bcSize - 1].value, systemMatrix, systemRHS);
//
//}
//
//
//
//int WriteResultFile(const char* fileName, Vector<double, Dynamic> result) {
//    ofstream ofile;
//    ofile.open(fileName);
//    if (!ofile.is_open())
//    {
//        cout << "Error! Could not open file!" << endl;
//        return -1;
//    }
//    for (int i = 0; i < result.size(); ++i) {
//        ofile << result(i) << endl;
//    }
//    return 0;
//}
//
//
//
//void Calculate(elem el, point* nodes, T* lambdaArray, T* qArray, Matrix<double, 3, 3>& localMatrix, Vector<double, 3>& localRHS) {
//
//    int nodeInd1 = el.nodesNum[0] - 1;
//    int nodeInd2 = el.nodesNum[1] - 1;
//    int nodeInd3 = el.nodesNum[2] - 1;
//    T le = norm(nodes[nodeInd3], nodes[nodeInd1]);
//
//    T delta = (nodes[nodeInd2].x - nodes[nodeInd1].x) * (nodes[nodeInd3].y - nodes[nodeInd1].y) - (nodes[nodeInd3].x - nodes[nodeInd1].x) * (nodes[nodeInd2].y - nodes[nodeInd1].y);
//    T ke = (lambdaArray[nodeInd1] + lambdaArray[nodeInd2] + lambdaArray[nodeInd3]) / 3.0;
//
//    T b1 = nodes[nodeInd2].y - nodes[nodeInd3].y;
//    T b2 = nodes[nodeInd3].y - nodes[nodeInd1].y;
//    T b3 = nodes[nodeInd1].y - nodes[nodeInd2].y;
//    T c1 = nodes[nodeInd3].x - nodes[nodeInd2].x;
//    T c2 = nodes[nodeInd1].x - nodes[nodeInd3].x;
//    T c3 = nodes[nodeInd2].x - nodes[nodeInd1].x;
//    T a1 = nodes[nodeInd2].x * nodes[nodeInd3].y - nodes[nodeInd3].x * nodes[nodeInd2].y;
//    T a2 = nodes[nodeInd3].x * nodes[nodeInd1].y - nodes[nodeInd1].x * nodes[nodeInd3].y;
//    T a3 = nodes[nodeInd1].x * nodes[nodeInd2].y - nodes[nodeInd2].x * nodes[nodeInd1].y;
//
//    T coeff = ke / (2 * delta);
//    localMatrix = Matrix<double, 3, 3>{ {b1 * b1 + c1 * c1, b1 * b2 + c1 * c2, b1 * b3 + c1 * c3},\
//                                        {b1* b2 + c1 * c2, b2* b2 + c2 * c2, b2* b3 + c2 * c3},\
//                                        {b1* b3 + c1 * c3, b2* b3 + c2 * c3, b3* b3 + c3 * c3} };
//    localMatrix *= coeff;
//
//    coeff = delta / 6. / 4.;
//
//    T q0 = qArray[nodeInd1]; T q1 = qArray[nodeInd2]; T q2 = qArray[nodeInd3];
//    localRHS = Vector<double, 3>{ 2 * q0 + q1 + q2,\
//                                  q0 + 2 * q1 + q2,\
//                                  q0 + q1 + 2 * q2 };
//    localRHS *= coeff;
//
//}
//
//int main()
//{
//    const char* reqfilename = "MeshRequest.txt";
//    const char* meshfilename = "Mesh.txt";
//    const char* bcfilename = "BoundaryConditions.txt";
//    const char* ofilename = "Results.txt";
//
//
//    Matrix<int, Dynamic, Dynamic> bcInds;
//    CreateMeshFile(reqfilename, meshfilename, bcInds);
//
//    int NE, NP, NC;
//    point* nodes = nullptr;
//    elem* elems = nullptr;
//    int** CPN = nullptr;
//    int* CPNsizes = nullptr;
//
//    ReadMeshFile(meshfilename, NE, NP, NC, nodes, elems, CPN, CPNsizes);
//    T* lambdaArray = new T[NP];
//    T* qArray = new T[NP];
//    Tabulate(lambda, nodes, NP, lambdaArray);
//    Tabulate(q, nodes, NP, qArray);
//
//    Matrix<double, Dynamic, Dynamic> systemMatrix = MatrixXd::Zero(NP, NP);
//    Vector<double, Dynamic> systemRHS = VectorXd::Zero(NP);
//    for (int i = 0; i < NE; ++i) {
//
//        Matrix<double, 3, 3> localMatrix;
//        Vector<double, 3> localRHS;
//        Calculate(elems[i], nodes, lambdaArray, qArray, localMatrix, localRHS);
//        Add(elems[i], localMatrix, localRHS, systemMatrix, systemRHS);
//
//    }
//
//    int bcSize = CPNsizes[0] + 4;
//    BoundaryCondition* bc = new BoundaryCondition[bcSize];
//    for (int i = 0; i < bcSize; ++i) bc[i].bcType = 0;
//    MakeBoundaryConditions(bc, bcSize, bcInds, nodes, 1, f1, 2);
//    MakeBoundaryConditions(bc, bcSize, bcInds, nodes, 2, f2, 1);
//    MakeBoundaryConditions(bc, bcSize, bcInds, nodes, 3, f3, 2);
//    MakeBoundaryConditions(bc, bcSize, bcInds, nodes, 4, f4, 1);
//    WriteBoundaryFile(bcfilename, bc, bcSize);
//    ReadBoundaryFile(bcfilename, bc, bcSize);
//
//    ApplyBoundaryValues(bc, bcSize, nodes, lambdaArray, systemMatrix, systemRHS);
//
//    //solve
//    Vector<double, Dynamic> result = systemMatrix.colPivHouseholderQr().solve(systemRHS);
//
//    WriteResultFile(ofilename, result);
//
//    cout << "Done!\n";
//    delete[] nodes;
//    for (int i = 0; i < NE; ++i) {
//        delete[] elems[i].nodesNum;
//    }
//    delete[] elems;
//    for (int i = 0; i < NC; ++i) {
//        delete[] CPN[i];
//    }
//    delete[] CPN;
//    delete[] CPNsizes;
//    delete[] bc;
//    return 0;
//}



#include <iostream>
#include <fstream>
#include <Eigen/Sparse>

using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::Triplet;
using Eigen::SparseLU;

double f(double x, double t) {
    return 0;
}

double k(double x, double t) {
    return 1;
}

double u_true(double x, double t) {
    return 7 * exp(-(4 * 3.1415) * (4 * 3.1415) * t) * sin(4 * 3.1415 * x);
}

double initial_condition(double x) {
    if (x >= 0 && x <= 0.3) {
        return x / 0.3;
    }
    return (1 - x) / 0.7;
}

template<typename T1, typename T2>
double scalar_product(const T1& vec1, const T2& vec2, double dx1, double dx2) {
    double res = 0;

    if (vec1.size() != vec2.size()) {
        std::cout << "The dimensions of the vectors do not match" << std::endl;
        return -1;
    }

    for (int i = 0; i < vec1.size(); i++) {
        res += vec1[i] * vec2[i];
    }

    return res * dx1 * dx2;
}


template<typename T>
void print_vector(const T& vec) {
    for (int i = 0; i < vec.size(); i++) {
        std:: cout << vec[i] << " ";
    }
    std::cout << std::endl;
}


void write_data_to_file(std::ofstream& file, double time, const VectorXd& u) {
    file << "Time = " << time << std::endl;

    for (int i = 0; i < u.size(); ++i) {
        file << u(i);
        if (i < u.size() - 1) {
            file << " ";
        }
    }
    file << std::endl;
}

VectorXd solve_heat_equation(double L, double T, int N, int M, double sigma, const VectorXd& v) {
    double dx = L / N;
    double dx2 = dx * dx;
    double dt = T / M;
    double l1 = 0;
    double l2 = 0;
    double alpha, beta, gamma;

    SparseMatrix<double> A(N + 1, N + 1);
    VectorXd b(N + 1), u_prev(N + 1);

    for (int i = 0; i < N + 1; i++) {
        u_prev(i) = v(i);
    }

    for (int m = 0; m < M; m++) {
        std::vector<Triplet<double>> triplets;
        triplets.reserve(3 * (N - 1) + 2);
        triplets.push_back(Triplet<double>(0, 0, 1.)); triplets.push_back(Triplet<double>(N, N, 1.));

        for (int i = 1; i < N; i++) {
            alpha = k(i * dx + dx / 2, (m + 1) * dt) * sigma / dx2;
            gamma = k(i * dx - dx / 2, (m + 1) * dt) * sigma / dx2;
            beta = -(1 / dt + (k(i * dx + dx / 2, (m + 1) * dt) + k(i * dx - dx / 2, (m + 1) * dt)) * sigma / dx2);

            triplets.push_back(Triplet<double>(i, i + 1, alpha));
            triplets.push_back(Triplet<double>(i, i, beta));
            triplets.push_back(Triplet<double>(i, i - 1, gamma));

            b(i) = k(i * dx + dx / 2, m * dt) * u_prev(i + 1) * (sigma - 1) / dx2 \
                + k(i * dx - dx / 2, m * dt) * u_prev(i - 1) * (sigma - 1) / dx2 \
                + f(i * dx, m * dt) \
                - (1 / dt + (k(i * dx + dx / 2, m * dt) + k(i * dx - dx / 2, m * dt)) * (sigma - 1) / dx2) * u_prev(i);
        }

        b(0) = l1;
        b(N) = l2;

        A.setFromTriplets(triplets.begin(), triplets.end());

        SparseLU<SparseMatrix<double>> solver;
        solver.compute(A);
        VectorXd u = solver.solve(b);

        u_prev = u;
    }

    return u_prev;
}

void reverse_solve(double iter_count, double L, double T, int N, int M, double sigma, std::ofstream& output_file) {
    VectorXd v(N + 1), y_N(N + 1), fi(N + 1);
    v.setZero();

    std::cout << "REVERSE SOLVE" << std::endl;

    for (int j = 0; j < N + 1; j++) {
        v(j) += initial_condition(j * L / N);
    }

    fi = solve_heat_equation(L, T, N, M, sigma, v);

    v.setZero();

    for (int i = 0; i < iter_count; i++) {
        std::cout << i << std::endl;

        y_N = solve_heat_equation(L, T, N, M, sigma, v);

        /*for (int j = 0; j < N + 1; j++) {
            v(j) += (fi(j) - y_N(j));
        }*/
        v += (fi - y_N); // ������ ����, ��� �������������� �������� �����������

        write_data_to_file(output_file, i, v);
    }

}

void reverse_solve_variation(double iter_count, double L, double T, int N, int M, double sigma, std::ofstream& output_file) {
    double s;
    VectorXd v(N + 1), y_N(N + 1), fi(N + 1), r_k(N + 1), Ar_k(N + 1);
    v.setZero();

    std::cout << "REVERSE SOLVE VARIATION" << std::endl;

    // ������ ������ ������ ��� �������� ��������
    for (int j = 0; j < N + 1; j++) {
        v(j) += initial_condition(j * L / N);
    }

    fi = solve_heat_equation(L, T, N, M, sigma, v);

    // ��������� �����������
    v.setZero();

    for (int i = 0; i < iter_count; i++) {
        std::cout << i << std::endl;

        y_N = solve_heat_equation(L, T, N, M, sigma, v);

        // ���������� � ������������ ������
        r_k = y_N - fi;

        Ar_k = solve_heat_equation(L, T, N, M, sigma, r_k);

        s = scalar_product(Ar_k, r_k, L / N, L / N) / scalar_product(Ar_k, Ar_k, L / N, L / N);

        v += s * (fi - y_N);

        write_data_to_file(output_file, i, v);
    }
}

int main() {
    double L = 1.0;
    double T = 0.1;
    int N = 10000;
    int M = 100;

    std::ofstream output_revers_var("output_revers_var.txt");
    std::ofstream output_revers("output_revers.txt");
    std::ofstream output_params("output_params.txt");

    output_params << "L = " << L << std::endl;
    output_params << "T = " << T << std::endl;
    output_params << "N = " << N << std::endl;
    output_params << "M = " << M << std::endl;

    reverse_solve(70, L, T, N, M, 0.5, output_revers);
    reverse_solve_variation(10, L, T, N, M, 0.5, output_revers_var);

    output_params.close();
    output_revers.close();

    return 0;
}




//#include <iostream>
//#include <fstream>
//#include <Eigen/Sparse>
//
//using namespace std;
//using namespace Eigen;
//
//class InverseHeatEquationSolver {
//public:
//    InverseHeatEquationSolver(double L, double T, int N, int M, double sigma, int iter_count)
//        : L(L), T(T), N(N), M(M), sigma(sigma), iter_count(iter_count) {}
//
//    void inverse_solve(const string& output_params_file, const string& output_revers_file) {
//        ofstream output_params(output_params_file);
//        ofstream output_revers(output_revers_file);
//
//        output_params << "L = " << L << endl;
//        output_params << "T = " << T << endl;
//        output_params << "N = " << N << endl;
//        output_params << "M = " << M << endl;
//
//        VectorXd v(N + 1), y_N(N + 1), fi(N + 1);
//        v.setZero();
//
//        // ������� ������ ������ ��� ���������� ������� ��
//        for (int j = 0; j < N + 1; j++) {
//            v(j) += initial_condition(j * L / N);
//        }
//        fi = solve_heat_equation(v);
//
//
//       /* write_data_to_file(output_params, 1, v);
//
//        write_data_to_file(output_params, 2, fi);*/
//
//        // ��������� �����������
//        v.setZero();
//
//
//        for (int i = 0; i < iter_count; i++) {
//            y_N = solve_heat_equation(v);
//
//            for (int j = 0; j < N + 1; j++) {
//                v(j) -= (fi(j) - y_N(j));
//            }
//
//            write_data_to_file(output_revers, i, v);
//        }
//
//        output_params.close();
//        output_revers.close();
//    }
//
//private:
//    double L, T, sigma;
//    int N, M;
//    int iter_count = 100; 
//
//    double f(double x, double t) {
//        return 0;
//    }
//
//    double initial_condition(double x) {
//        if (x >= 0 && x <= 0.3) {
//            return x / 0.3;
//        }
//        return (1 - x) / 0.7;
//    }
//
//    /*double final_condition(double x) {
//        if (x >= 0 && x <= 0.3) {
//            return x / 0.3;
//        }
//        return (1 - x) / 0.7;
//    }*/
//
//    void write_data_to_file(ofstream& file, double time, const VectorXd& u) {
//        file << "Time = " << time << endl;
//
//        for (int i = 0; i < u.size(); ++i) {
//            file << u(i);
//            if (i < u.size() - 1) {
//                file << " ";
//            }
//        }
//        file << endl;
//    }
//
//    VectorXd solve_heat_equation(const VectorXd& v) {
//        double dx = L / N;
//        double dx2 = dx * dx;
//        double dt = T / M;
//        double l1 = 0;
//        double l2 = 0;
//
//        SparseMatrix<double> A(N + 1, N + 1);
//        VectorXd b(N + 1), u_prev(N + 1);
//
//        for (int i = 0; i < N + 1; i++) {
//            u_prev(i) = v(i);
//        }
//
//        for (int m = 0; m < M; m++) {
//            A.setZero();
//            b.setZero();
//
//            for (int i = 1; i < N; i++) {
//                double alpha = sigma / dx2;
//                double gamma = sigma / dx2;
//                double beta = -(1 / dt + 2 * sigma / dx2);
//
//                A.coeffRef(i, i) += beta;
//                A.coeffRef(i, i - 1) += gamma;
//                A.coeffRef(i, i + 1) += alpha;
//
//                b(i) = f(i * dx, m * dt)
//                    - (1 / dt + (sigma - 1) / dx2) * u_prev(i);
//            }
//
//            b(0) = l1;
//            b(N) = l2;
//
//            BiCGSTAB<SparseMatrix<double>> solver;
//            solver.compute(A);
//            VectorXd u = solver.solve(b);
//
//            u_prev = u;
//        }
//
//        return u_prev;
//    }
//};
//
//int main() {
//    double L = 1.0;
//    double T = 0.1;
//    int N = 100;
//    int M = 100;
//    double sigma = 0.5;
//
//    InverseHeatEquationSolver solver(L, T, N, M, sigma, 10);
//    solver.inverse_solve("output_params.txt", "output_revers.txt");
//
//    return 0;
//}





//void solve_heat_equation(double L, double T, int N, int M, double sigma, ofstream& output_file) {
//    double dx = L / N;
//    double dx2 = dx * dx;
//    double dt = T / M;
//    double l1 = boundary_condition_0();
//    double l2 = boundary_condition_1();
//    double alpha, beta, gamma;
//
//    SparseMatrix<double> A(N + 1, N + 1);
//    VectorXd b(N + 1), u_prev(N + 1);
//
//    // ��������� �������
//    for (int i = 0; i < N + 1; i++) {
//        u_prev(i) = initial_condition(i * dx);
//    }
//
//    // ������� ��� ������� ���������� ����
//    for (int m = 0; m < M; m++) {
//        // ���������� ������� ������� ��������� � ������� ������ �����
//        vector<Triplet<double>> triplets;
//        triplets.reserve(3 * (N-1) + 2);
//        triplets.push_back(Triplet<double>(0, 0, 1.)); triplets.push_back(Triplet<double>(N, N, 1.));
//
//        for (int i = 1; i < N; i++) {
//            alpha = k(i * dx + dx / 2, (m + 1) * dt) * sigma / dx2;
//            gamma = k(i * dx - dx / 2, (m + 1) * dt) * sigma / dx2;
//            beta = -(1 / dt + (k(i * dx + dx / 2, (m + 1) * dt) + k(i * dx - dx / 2, (m + 1) * dt)) * sigma / dx2);
//
//            triplets.push_back(Triplet<double>(i, i + 1, alpha));
//            triplets.push_back(Triplet<double>(i, i, beta));
//            triplets.push_back(Triplet<double>(i, i - 1, gamma));
//
//            b(i) = k(i * dx + dx / 2, m * dt) * u_prev(i + 1) * (sigma - 1) / dx2 \
//                + k(i * dx - dx / 2, m * dt) * u_prev(i - 1) * (sigma - 1) / dx2 \
//                + f(i * dx, m * dt) \
//                - (1 / dt + (k(i * dx + dx / 2, m * dt) + k(i * dx - dx / 2, m * dt)) * (sigma - 1) / dx2) * u_prev(i);
//        }
//
//        // ���� ��������� �������
//        b(0) = l1;
//        b(N) = l2;
//
//        A.setFromTriplets(triplets.begin(), triplets.end());
//
//        // ������� ������� ���������
//        SparseLU<SparseMatrix<double>> solver;
//        solver.compute(A);
//        VectorXd u = solver.solve(b);
//
//        write_data_to_file(output_file, m * dt, u_prev);
//
//        //// ��������� � �������� ��������
//        // cout << "Time = " << m * dt << ", norm = " << norm(dx, m * dt, u_prev) << endl;
//
//        // ���������� ����������� ������� ��� ���������� ���������� ����
//        u_prev = u;
//    }
//    write_data_to_file(output_file, M * dt, u_prev);
//}
//
//int main() {
//    double Lengt = 1.0; // ������� ������� �� ������������ (0; Lengt)
//    double Time = 1.0; // ������� ������� (0; Time)
//    int N = 100; // ����� �� ������������ 
//    int M = 1000; // ����� �� �������
//
//    ofstream output_params("output_params.txt");
//    ofstream output_results("output_results.txt");
//
//    output_params << "L = " << Lengt << endl;
//    output_params << "T = " << Time << endl;
//    output_params << "N = " << N << endl;
//    output_params << "M = " << M << endl;
//    solve_heat_equation(Lengt, Time, N, M, 0.5, output_results);
//
//    output_results.close();
//    output_params.close();
//
//    return 0;
//}

