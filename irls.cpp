#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

float madsigma(VectorXcf residuals, int rank)
{
    float sigma = 0.6745;
    VectorXf residual_magnitues = residuals.array().abs();
    sort(residual_magnitues.data(), residual_magnitues.data()+residual_magnitues.size());
    int middle = (rank - 1 + residual_magnitues.size()) / 2;
    float median = residual_magnitues(middle);
    return median / sigma;
}

VectorXf tukey_biweight(VectorXcf residuals)
{
    VectorXf abs_weights = residuals.array().abs();
    VectorXf weights = residuals.array().abs();
    for (int i = 0; i < weights.size(); i++) {
        if (abs_weights(i) > 1) {
            weights(i) = 0;
        }
    }
    return (1 - weights.array().square()).square();
}

VectorXcf irls(MatrixXcf x, VectorXcf y)
{
    int n = x.rows();
    int p = x.cols();

    VectorXcf prior_weights = VectorXcf::Ones(n);
    VectorXcf curr_weights = VectorXcf::Ones(n);
    VectorXcf new_weights;

    ColPivHouseholderQR<MatrixXcf> qr = x.colPivHouseholderQr();
    VectorXcf b = qr.solve(y);
    VectorXcf b0;

    MatrixXcf r = qr.matrixR().triangularView<Upper>();
    // B/A = (A'\B')'
    MatrixXcf e =  r.transpose().colPivHouseholderQr().
        solve((x * qr.colsPermutation()).transpose()).
        transpose();
    VectorXcf h = (e.array().square()).rowwise().sum();
    VectorXcf adj_factor = 1 / (1 - h.array() / prior_weights.array()).sqrt();
    
    int wxrank = qr.rank();
    float D = 1.4901e-08;

    int iter_limit = 50;
    VectorXcf res, r_adj;
    for (int i = 0; i < iter_limit; i++) {
        if (b0.size()) {
            // test threshhold
        }

        res = y - x * b;
        r_adj = (res.array() * adj_factor.array()) / curr_weights.array();
        float s = madsigma(r_adj, wxrank);

        float tune = 4.685;
        VectorXf weights = tukey_biweight(r_adj.array() / (s * tune));
        b0 = b;
        weights = weights.array().sqrt();
        VectorXcf y_weighted = y.array() * weights.array();
        MatrixXcf x_weighted = x.array() * (weights * VectorXcf::Ones(p).transpose()).array();
        qr = x_weighted.colPivHouseholderQr();
        wxrank = qr.rank();
        b = qr.solve(y_weighted);
    }
    //cout << b << endl;
    return b;
}

MatrixXcf loadMatrix(string file_name, int m, int n)
{
    MatrixXcf a = MatrixXcf::Zero(m, n);
    string line;
    ifstream f2 (file_name);
    int i = 0;
    while ( getline (f2, line) )
    {
        int j = 0;
        stringstream ss(line);
        string token;
        float c;
        while ( getline (ss, token, ',') ) {
            c = ::stof(token);
            a(i, j) = c;
            j++;
        }
        i++;
    }
    f2.close();
    return a;
}

VectorXcf loadVector(string file_name, int m)
{
    VectorXcf b = VectorXcf::Zero(m);
    string line;
    ifstream f2 (file_name);
    int i = 0;
    while ( getline (f2, line) )
    {
        float c = ::stof(line);
        b(i) = c;
        i++;
    }
    f2.close();
    return b;
}

MatrixXcf loadComplexMatrix(string file_name, int m, int n)
{
    MatrixXcf b = MatrixXcf::Zero(m, n);
    string line;
    ifstream f2 (file_name);
    int i = 0;
    while ( getline (f2, line) )
    {
        int j = 0;
        stringstream ss(line);
        string token;
        float c;
        while ( getline (ss, token, ' ') ) {
            c = ::stof(token);
            if (j % 2 == 0) {
                b.real()(i, j/2) = c;
            } else {
                b.imag()(i, j/2) = c;
            }
            j++;
        }
        i++;
    }
    f2.close();
    return b;
}

VectorXcf loadComplexVector(string file_name, int m)
{
    VectorXcf b = VectorXcf::Zero(m);
    string line;
    ifstream f2 (file_name);
    int i = 0;
    while ( getline (f2, line) )
    {
        int j = 0;
        stringstream ss(line);
        string token;
        float c;
        while ( getline (ss, token, ' ') ) {
            c = ::stof(token);
            if (j == 0) {
                b.real()(i) = c;
            } else {
                b.imag()(i) = c;
            }
            j++;
        }
        i++;
    }
    f2.close();
    return b;
}

//int _main()
//{

    //MatrixXcf a = loadMatrix("data/F2_real.csv", 402, 120);
    //VectorXcf y = loadVector("data/y_real.csv", 402);

    //ofstream a_test ("data/a_test");
    //a_test << a;
    //ofstream y_test ("data/y_test");
    //y_test << y;
    
    //std::clock_t start, end;
    //start = std::clock();
    //VectorXcf b = irls(a, y);
    //end = std::clock();
    //ofstream b_test ("data/b_test");
    //b_test << b;

    //cout << "Took: " << (end - start) / (double)(CLOCKS_PER_SEC / 1000) << "ms\n" << endl;

    //return 0;
//}

