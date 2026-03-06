#include "io.h"
#include "expressions.h"
#include <iostream>
int main()
{
    // 1. Read TXT
    BeamData d = readInput("input.txt");

    // 2. Calculate
    // scalar
    double delta = computeDeflection(d, d.I, d.E);
    double fx = computef(d, d.IX, d.EC);
    double fy = computef(d, d.IY, d.EC);
    //  matrices
    auto SA = createSAMatrix(d.nFloors);
    auto Sx  = multiplyscamat(fx, SA);
    auto Sy  = multiplyscamat(fy, SA);
//    auto M = buildMassMatrix(d);
std::vector<std::vector<double>> M(13, std::vector<double>(13, 0.0));

    double diag[13] =
    {
        539.7588,
        526.1904,
        524.1307,
        522.0564,
        522.3585,
        507.0437,
        490.7714,
        492.1757,
        492.0478,
        454.4963,
        380.9595,
        330.399,
        6.271963
    };

    for(int i = 0; i < 13; i++)
        M[i][i] = diag[i];

    auto Dx = multiplymatmat(Sx, M);
    auto Dy = multiplymatmat(Sy, M);
    Matrix phi,LA;
    eigenDecomposition(Dx, phi, LA);
    // 3. Write results
    writeOutput("output.txt", delta, fx, fy, SA, Sx, Sy, M, Dx, Dy, phi, LA);
    return 0;
}
