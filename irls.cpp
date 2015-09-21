#include <ctime>
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

float madsigma(VectorXf residuals, int rank)
{
    float sigma = 1.4826;
    sort(residuals.data(), residuals.data()+rank);
    float median = residuals(residuals.size() / 2);
    return median * sigma;
}

VectorXf tukey_biweight(VectorXf residuals)
{
    VectorXf abs_weights = residuals.array().abs();
    VectorXf weights = residuals;
    for (int i = 0; i < weights.size(); i++) {
        if (abs_weights(i) > 1) {
            weights(i) = 0;
        }
    }
    return (1 - weights.array().square()).square();
}

void irls(MatrixXf x, VectorXf y)
{
    int n = x.rows();
    int p = x.cols();

    VectorXf prior_weights = VectorXf::Ones(n);
    VectorXf curr_weights = VectorXf::Ones(n);
    VectorXf new_weights;

    ColPivHouseholderQR<MatrixXf> qr = x.colPivHouseholderQr();
    VectorXf b = qr.solve(y);
    VectorXf b0;

    MatrixXf r = qr.matrixR();
    MatrixXf e = x.array() / r.array();
    //float h = (e.array() * e.array()).rowwise().sum().minCoeff();
    VectorXf adj_factor = 100 * VectorXf::Ones(n);
    cout << "E" << endl;
    cout << e << endl << endl;
    cout << "R" << endl;
    cout << r << endl << endl;

    
    int wxrank = qr.rank();
    float D = 1.4901e-08;

    int iter_limit = 50;
    VectorXf res, r_adj;
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
        VectorXf y_weighted = y.array() * weights.array();
        MatrixXf x_weighted = x.array() * (weights * VectorXf::Ones(p).transpose()).array();
        qr = x_weighted.colPivHouseholderQr();
        wxrank = qr.rank();
        b = qr.solve(y_weighted);
    }
    cout << b << endl;
}

int main()
{
    //MatrixXf a = MatrixXf::Random(302, 120);
    //VectorXf b = VectorXf::Random(302);
    
    MatrixXf a(5,3);
    a << 7,     1,    10,
        2,     1,     1,
        8,     9,     5,
        1,     7,     4,
        3,     4,     8;
    VectorXf b(5);
    b << 2,
     3,
     9,
     3,
     9;

    /*
       Looking for something like:
       1.3285
       0.1326
       0.5602
       0.1517
       */

    std::clock_t start, end;

    start = std::clock();
    irls(a, b);
    end = std::clock();


    cout << "Took: " << (end - start) / (double)(CLOCKS_PER_SEC / 1000) << "ms\n" << endl;
}

