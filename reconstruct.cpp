#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

#include "irls.cpp"

VectorXcf extractKnown(VectorXcf spectrum, VectorXcf valid_points)
{
    ArrayXf a = valid_points.real().array();
    VectorXcf extracted = VectorXcf::Zero((a > 0).count());

    int j = 0;
    for (int i = 0; i < a.rows(); i++) {
        if (a(i)) {
            extracted(j) = spectrum(j);
            j++;
        }
    }
    return extracted;
}

void reconstruct(VectorXcf spectrum, VectorXcf valid_points)
{
    VectorXcf hamming_320 = loadComplexVector("data/hamming_320.dat", 320);
    VectorXcf invalid_points = VectorXcf::Ones(320).array() - valid_points.array();
    VectorXcf channel_impulse_response_valid_points = VectorXcf::Ones(320);

    float stop_threshold = 0.03;
    int stop_left = 0;
    int stop_right = 0;
    int last_run = 0;
    int left_cursor = 100;
    int right_cursor = 220;
    MatrixXcf idft = loadComplexMatrix("data/idft.dat", 320, 320);
    VectorXcf spectrum_windowed = spectrum.array() * hamming_320.array();
    VectorXcf known = extractKnown(spectrum_windowed, valid_points);

    while (!stop_left && !stop_right) {

    }



    //cout << spectrum_windowed << endl;
}

int main()
{
    cout << "hello world from reconstruct" << endl;
    VectorXcf spectrum = loadComplexVector("data/freq_resp.dat", 320);
    VectorXcf valid_points = loadComplexVector("data/valid_points.dat", 320);

    reconstruct(spectrum, valid_points);

    cout << "done" << endl;
}
