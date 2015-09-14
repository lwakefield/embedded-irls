#include <ctime>
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
int main()
{
    MatrixXcf a = MatrixXcf::Random(302,120);

    std::clock_t start, end;

    start = std::clock();
    a.householderQr();
    end = std::clock();

    cout << "Took: " << (end - start) / (double)(CLOCKS_PER_SEC / 1000) << "ms\n" << endl;
}
