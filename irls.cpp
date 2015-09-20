#include <ctime>
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;



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

    //MatrixXf r = qr.matrixQR().triangularView<Upper>();
    //MatrixXf e = x.array() / r.array();
    //float h = (e.array() * e.array()).rowwise().sum().minCoeff();
    VectorXf adj_factor = 100 * VectorXf::Ones(n);
    
    int wxrank = qr.rank();
    float D = 1.4901e-08;

    int iter_limit = 50;
    VectorXf r, r_adj;
    for (int i = 0; i < iter_limit; i++) {
        if (b0.size()) {
            // test threshhold
        }

        r = y - x * b;
        r_adj = (r.array() * adj_factor.array()) / curr_weights.array();
        new_weights = madsigma(r_adj, wxrank);
    }
}

int main()
{
    MatrixXf a = MatrixXf::Random(302,120);
    VectorXf b = VectorXf::Random(302);

    std::clock_t start, end;

    start = std::clock();
    irls(a, b);
    end = std::clock();

    cout << "Took: " << (end - start) / (double)(CLOCKS_PER_SEC / 1000) << "ms\n" << endl;
}

