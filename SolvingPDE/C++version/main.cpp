#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <set>
#include <tuple>
#include <algorithm>

using namespace std;

double h = 0.5;
int N = round(1/h);


// left part of equation
double f(double x,double y) {
    return 1.1*sin(x) + (3.2*pow(x,2) + 4.4*pow(y,2))*cos(2*x*y);
}

// exact solution
double u(double x, double y) {
    return sin(x) + cos(2*x*y);
}

// set of indices which covers initial area (WITHOUT BOUNDARY)!
set<tuple<int,int>> get_indices() {
    set<tuple<int,int>> idx;

    // we know, that x in [0, 0.5] and y in [0, 0.5] is always in our area
    // STEP 1: Rectangle y in [0,1], x in [0,0.5]
    for (int i=1; i<N/2; i++) {
        for (int j=1; j<N; j++) {
            idx.insert(tuple<int, int>{i, j});
        }
    }

    // STEP 2: Square y in [0, 0.5] x in [0.5, 1]
    for (int i=N/2; i<N; i++) {
        for (int j=1; j<N/2; j++) {
            idx.insert(tuple<int, int>{i, j});
        }
    }

    // STEP 3: Upper right triangle
    for (int i=N/2; i<N; i++) {
        int indicator = i - N/2;
        for (int j=N/2; j<N-indicator; j++) {
            idx.insert(tuple<int, int>{i, j});
        }
    }
    return idx;
}

// Set of indices for boundary
set<tuple<int,int>> get_boundary() {
    set<tuple<int, int>> b;
    // lower and left boundaries
    for(int j=0; j<=N; j++) {
        b.insert(tuple<int, int>(0, j));
        b.insert(tuple<int,int>(j, 0));
    }

    // upper and right boundaries
    for (int j=0; j<=N/2; j++) {
        b.insert(tuple<int, int>(N, j));
        b.insert(tuple<int,int>(j, N));
    }

    // triangle part of boundary
    for (int j=N/2; j<=N; j++) {
        int indicator = j - N/2;
        b.insert(tuple<int, int>(j, N-indicator));
    }
    return b;
}

// create matrix of size x size with val
double** create_matrix(int size, double val) {
    auto** Y = new double*[N+1];
    for (int i=0; i<=N; i++) {
        Y[i] = new double[N+1];
        for (int j=0; j<=N; j++) {
            Y[i][j] = val;
        }
    }
    return Y;
}

// solution matrix
double** solution_matrix(set<tuple<int, int>>* s, set<tuple<int, int>>* b) {
    double** U = create_matrix(N+1, 0);
    set<tuple<int,int>> idx;
    set<tuple<int, int>>::iterator it;
    // we don't need in separated boundary and inner indices, so merge our sets
    merge(s->begin(), s->end(),
              b->begin(), b->end(),
              inserter(idx, idx.begin()));
    for (it = idx.begin(); it != idx.end(); it++) {
        int i = get<0>(*it);
        int j = get<1>(*it);
        U[i][j] = u(i * h, j * h);
    }
    return U;
}

// dot product of two matrices size of N
double dot_product(double** A, double** B) {
    double res = 0;
    for(int i=0; i<=N; i++) {
        for(int j=0; j<=N; j++) {
            res += A[i][j] * B[i][j];
        }
    }
    return res;
}

// norm of matrix size of N
double norm(double** A) {
   return pow(dot_product(A, A), 0.5);
}

// normalize matrix
double** normalize(double** A) {
    double norm_val = norm(A);
    if (norm_val == 0) throw runtime_error("error");
    for(int i=0; i<=N; i++) {
        for(int j=0; j<=N; j++) {
            A[i][j] = A[i][j]/norm_val;
        }
    }
    return A;
}

// remove matrix from memory
void clean_matrix(double** Y) {
    for (int i=0; i<=N; i++)
        delete [] Y[i];
    delete [] Y;
}

// five points operator
double** operator_step(double** U, set<tuple<int, int>>* s) {
    double** Y = create_matrix(N+1, 0);
    set<tuple<int,int>>::iterator iter;
    for(iter = s->begin(); iter != s->end(); iter++) {
        int i = get<0>(*iter);
        int j = get<1>(*iter);
        Y[i][j] = -1.1*(U[i-1][j] - 2*U[i][j] + U[i+1][j])/pow(h,2) - 0.8*(U[i][j-1] - 2*U[i][j] + U[i][j+1])/pow(h,2);
    }
    return Y;
}

// five points operator with f
double** operator_step_f(double** U, set<tuple<int, int>>* s) {
    double** Y = create_matrix(N+1, 0);
    set<tuple<int,int>>::iterator iter;
    for(iter = s->begin(); iter != s->end(); iter++) {
        int i = get<0>(*iter);
        int j = get<1>(*iter);
        Y[i][j] = -1.1*(U[i-1][j] - 2*U[i][j] + U[i+1][j])/pow(h,2) - 0.8*(U[i][j-1] - 2*U[i][j] + U[i][j+1])/pow(h,2) - f(i*h, j*h);
    }
    return Y;
}

double** operator_step_reversed(double** U, double lambda_max, set<tuple<int,int>>* s) {
    double** Y = create_matrix(N+1, 0);
    set<tuple<int,int>>::iterator iter;
    for(iter = s->begin(); iter != s->end(); iter++) {
        int i = get<0>(*iter);
        int j = get<1>(*iter);
        Y[i][j] = lambda_max*U[i][j] + 1.1*(U[i-1][j] - 2*U[i][j] + U[i+1][j])/pow(h,2) + 0.8*(U[i][j-1] - 2*U[i][j] + U[i][j+1])/pow(h,2);
    }
    return Y;
}

// copy matrix
double** copy_matrix(double** Y) {
    auto** Z = new double*[N+1];
    for (int i=0; i<=N; i++) {
        Z[i] = new double[N+1];
        for (int j = 0; j <= N; j++)
            Z[i][j] = Y[i][j];
    }
    return Z;
}

// finding max and min eigenvalues
tuple<double,int,double,int> min_max_eigenvalues(double delta, set<tuple<int,int>>* s) {
    // initial matrix and max eigenvalue
    double** Y = create_matrix(N+1, 0);
    set<tuple<int,int>>::iterator iter;
    for(iter = s->begin(); iter != s->end(); iter++) {
        int i = get<0>(*iter);
        int j = get<1>(*iter);
        Y[i][j] = 1;
    }

    Y = normalize(Y);
    // double ik = Y[0][0]
    // double k = -1.1*(Y[0][1] - 2*Y[1][1] + Y[2][1])/pow(h,2) - 0.8*(Y[1][0] - 2*Y[1][1] + Y[1][2])/pow(h,2);
    double** Z = operator_step(Y, s);

    double lambda_prev = 1;
    double lambda_max = dot_product(Z, Y);
    int max_count = 0;

    while (fabs(1-lambda_max/lambda_prev) > delta) {
        max_count += 1;
        Y = copy_matrix(Z);
        Y = normalize(Y);
        Z = operator_step(Y, s);

        lambda_prev = lambda_max;
        lambda_max = dot_product(Y,Z);
    }

    // Now for lambda min
    Y = create_matrix(N+1, 0);
    for(iter = s->begin(); iter != s->end(); iter++) {
        int i = get<0>(*iter);
        int j = get<1>(*iter);
        Y[i][j] = 1;
    }
    Y = normalize(Y);
    Z = operator_step_reversed(Y, lambda_max, s);

    lambda_prev = 1;
    double lambda_min = dot_product(Z,Y);
    int min_count = 0;

    while (fabs(1-lambda_min/lambda_prev) > delta && norm(Z) != 0) {
        min_count += 1;
        Y = copy_matrix(Z);
        Y = normalize(Y);
        Z = operator_step_reversed(Y, lambda_max, s);

        lambda_prev = lambda_min;
        lambda_min = dot_product(Y,Z);
    }
    lambda_min = lambda_max - lambda_min;
    return {lambda_max, max_count, lambda_min, min_count};
}

double next_step(double** U, int i, int j) {
    return -1.1*(U[i-1][j] - 2*U[i][j] + U[i+1][j])/pow(h,2) - 0.8*(U[i][j-1] - 2*U[i][j] + U[i][j+1])/pow(h,2);
}

double** multiply(double** Y, double val) {
    for (int i=0; i<=N; i++)
        for (int j=0; j<=N; j++)
            Y[i][j] = Y[i][j]*val;
    return Y;
}

double** subtract(double** X1, double** X2) {
    double** Y = create_matrix(N+1, 0);
    for (int i=0; i<=N; i++)
        for (int j=0; j<=N; j++)
            Y[i][j] = X1[i][j]-X2[i][j];
    return Y;
}

// Find solution of system using gradient descent method
double** solve_grad(double delta, set<tuple<int,int>>* s, set<tuple<int,int>>* b) {
    double** U = create_matrix(N+1, 0);
    // fill matrix of solution according to boundary condition
    set<tuple<int, int>>::iterator it;
    for (it=b->begin(); it != b->end(); it++) {
        int i = get<0>(*it);
        int j = get<1>(*it);
        U[i][j] = u(i*h, j*h);
    }
    for (it=s->begin(); it != s->end(); it++) {
        int i = get<0>(*it);
        int j = get<1>(*it);
        U[i][j] = 1;
    }

    // Number of iterations
    int count = 0;

    double** R0 = create_matrix(N+1, 0);
    double** R1 = operator_step_f(U, s);

    while (norm(subtract(R0, R1))*h > delta || count == 0) {
        count += 1;
        R0 = copy_matrix(R1);
        double tau = (dot_product(R0, R0)/dot_product(operator_step(R0, s), R0));
        U = subtract(U, multiply(R0, tau));
        R1 = operator_step_f(U, s);
    }
    clean_matrix(R0);
    clean_matrix(R1);
    return U;
}

// five points operator for jacoby method
double** operator_step_ja(double** U, set<tuple<int, int>>* s, set<tuple<int, int>>* b) {
    double** Y = create_matrix(N+1, 0);
    set<tuple<int,int>>::iterator iter;
    for(iter = s->begin(); iter != s->end(); iter++) {
        int i = get<0>(*iter);
        int j = get<1>(*iter);
        Y[i][j] = (1/(2*1.9)) * (1.1*(U[i-1][j] + U[i+1][j]) + 0.8*(U[i][j-1] + U[i][j+1]) + pow(h,2)*f(i*h, j*h));
    }
    set<tuple<int, int>>::iterator it;
    for (it=b->begin(); it != b->end(); it++) {
        int i = get<0>(*it);
        int j = get<1>(*it);
        Y[i][j] = u(i*h, j*h);
    }
    return Y;
}

double** solve_jacobi(double delta, set<tuple<int,int>>* s, set<tuple<int,int>>* b) {
    double** U = create_matrix(N+1, 0);
    // fill matrix of solution according to boundary condition
    set<tuple<int, int>>::iterator it;
    for (it=b->begin(); it != b->end(); it++) {
        int i = get<0>(*it);
        int j = get<1>(*it);
        U[i][j] = u(i*h, j*h);
    }
    for (it=s->begin(); it != s->end(); it++) {
        int i = get<0>(*it);
        int j = get<1>(*it);
        U[i][j] = 1;
    }

    // Number of iterations
    int count = 0;

    double ** Z = operator_step_ja(U, s, b);
    while (norm(subtract(Z,U))*h > delta || count == 0) {
        count += 1;
        U = copy_matrix(Z);
        Z = operator_step_ja(U,s,b);
    }
    return U;
}

// returns arrays (x,y) for plots, val is in [0,1]
tuple<double*, double*> get_arrays(double **Y, float val) {
    int k = round(val/h);
    auto* x = new double[N+1];
    auto* func_val = new double[N+1];

    for (int i=0; i<=N; i++) {
        x[i] = i * h;
        func_val[i] = Y[k][i];
    }
    return {x, func_val};
}

int main() {
    set<tuple<int, int>> s = get_indices();
    set<tuple<int, int>> b = get_boundary();

    // show the greatest eigenvalue
    auto eigenvalues = min_max_eigenvalues(1e-9, &s);
    cout << "Min eigenvalue: " << get<2>(eigenvalues) << '\n' << "Max eigenvalue: " << get<0>(eigenvalues) << endl;

    // Solution of gradient descent
    double** Sol = solve_grad(1e-5, &s, &b);
    double** U = solution_matrix(&s, &b);
    cout << "Error of gradient descent method: " << norm(subtract(Sol, U))*h << endl;

    // Solution of Jacobi method
    double** Sol_ja = solve_jacobi(1e-2, &s, &b);
    cout << "Error of Jacobi method: " << norm(subtract(Sol_ja, U))*h << endl;

    // data for plots (2d plots, gradient descent)
//    auto approx_sol = get_arrays(Sol, 0.5);
//    auto exact_solution = get_arrays(U, 0.5);
//    double* x = get<0>(approx_sol);
//    double* y_approx = get<1>(approx_sol);
//    double* y_exact = get<1>(exact_solution);
//    ofstream file;
//    file.open("data_for_plots.txt");
//    for (int i=0; i<=N; i++)
//        file << x[i] << ' ';
//    file << endl;
//
//    for (int i=0; i<=N; i++)
//        file << y_approx[i] << ' ';
//    file << endl;
//
//
//    for (int i=0; i<=N; i++)
//        file << y_exact[i] << ' ';
//    file << endl;
//
//    file.close();

    // data for plots (3d plots, gradient descent)
    ofstream file;
    file.open("approximate_solution.txt");
    for (int i=0; i<=N; i++) {
        for (int j = 0; j <= N; j++)
            file << Sol_ja[i][j] << ' ';
        file << endl;
    }
    file.close();
    file.open("exact_solution.txt");
    for (int i=0; i<=N; i++) {
        for (int j = 0; j <= N; j++)
            file << U[i][j] << ' ';
        file << endl;
    }
    file.close();
    return 0;
}
