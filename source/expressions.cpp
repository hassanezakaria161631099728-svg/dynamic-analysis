#include "expressions.h"
#include <cmath>

double computeDeflection(const BeamData& d, double I, double E)

{
        // simple beam formula example
        double delta = (5 * d.w * pow(d.L,4)) /(384 * E * I);
    return delta;
}

double computef(const BeamData& d, double I, double E)

{
        // expression for a factor we multiply on SA to get stiffness matrix
        double f = (pow(d.h,3)) / (6 * E * I );
    return f;
}

std::vector<std::vector<double>>
buildMassMatrix(const BeamData& d)
{
    int n = d.nFloors;

    // initialize all zeros
    std::vector<std::vector<double>> M(n,
        std::vector<double>(n, 0.0));

    // ---- RULE 1 ----
    // for i = 1 → n-2
    if(n >= 3)
    {
        for(int i=1; i<=n-2; i++)
        {
            M[i-1][i-1] = d.M1;
        }

        // ---- RULE 2 ----
        M[n-2][n-2] = d.M2;

        // ---- RULE 3 ----
        M[n-1][n-1] = d.M3;
    }

    return M;
}

std::vector<std::vector<double>>
createSAMatrix(int n)
{
    // create n x n matrix filled with 0
    std::vector<std::vector<double>> A(n,
        std::vector<double>(n, 0.0));

    for(int i=1; i<=n; i++)
    {
        for(int j=1; j<=n; j++)
        {
            if(i < j)
                A[i-1][j-1] = i*i * (3*j - i);

            else if(i == j)
                A[i-1][j-1] = 2*i*i*i;

            else // i > j
                A[i-1][j-1] = j*j * (3*i - j);
        }
    }

    return A;
}

std::vector<std::vector<double>>
multiplymatmat(const std::vector<std::vector<double>>& A,
         const std::vector<std::vector<double>>& B)
{
    int n = A.size();

    std::vector<std::vector<double>> C(n,
        std::vector<double>(n, 0.0));

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            for(int k=0; k<n; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}

std::vector<std::vector<double>>
multiplyscamat(double a,
         const std::vector<std::vector<double>>& B)
{
    int n = B.size();

    std::vector<std::vector<double>> C(n,
        std::vector<double>(n));

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
                C[i][j] = a * B[i][j];
            }
        }

    return C;
}

// ===== Helpers =====

static std::vector<double>
matVec(const Matrix& A,
       const std::vector<double>& x)
{
    int n = A.size();
    std::vector<double> r(n,0.0);

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            r[i] += A[i][j]*x[j];

    return r;
}

static double dot(const std::vector<double>& a,
                  const std::vector<double>& b)
{
    double s=0;
    for(size_t i=0;i<a.size();i++)
        s += a[i]*b[i];
    return s;
}

static void normalize(std::vector<double>& x)
{
    double n = std::sqrt(dot(x,x));
    if(n < 1e-12) return;

    for(double& v : x)
        v /= n;
}

// ===== POWER ITERATION → eigenvector + eigenvalue =====

static void
powerMode(const Matrix& A,
          std::vector<double>& vec,
          double& lambda)
{
    int n = A.size();

    vec.assign(n,1.0);
    normalize(vec);

    for(int k=0;k<200;k++)
    {
        auto y = matVec(A,vec);

        lambda = dot(vec,y);      // Rayleigh

        vec = y;
        normalize(vec);
    }
}

// ===== Deflation =====

static void
deflate(Matrix& A,
        const std::vector<double>& v,
        double lambda)
{
    int n = A.size();

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            A[i][j] -= lambda * v[i] * v[j];
}

// ===== MAIN =====

void eigenDecomposition(
    const Matrix& A,
    Matrix& B,
    Matrix& C)
{
    int n = A.size();

    B.assign(n, std::vector<double>(n,0.0));
    C.assign(n, std::vector<double>(n,0.0));

    Matrix work = A;

    for(int k=0;k<n;k++)
    {
        std::vector<double> vec;
        double lambda;

        powerMode(work, vec, lambda);

        // store eigenvector as column k
        for(int i=0;i<n;i++)
            B[i][k] = vec[i];

        // store eigenvalue on diagonal
        C[k][k] = lambda;

        // remove this mode
        deflate(work, vec, lambda);
    }
}

