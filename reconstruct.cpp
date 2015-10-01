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

VectorXcf buildZeroedPoints(int left_cursor, int right_cursor, int size)
{
    VectorXcf left = VectorXcf::Ones(left_cursor);
    VectorXcf middle = VectorXcf::Zero(size - left_cursor - right_cursor);
    VectorXcf right = VectorXcf::Ones(right_cursor);
    VectorXcf joined;
    joined << left, middle, right;
    return joined;
}

int getInvalidStart(VectorXcf valid_points)
{
    for (int i = 0; i < valid_points.rows(); i++) {
        if (valid_points.real()(i) == 0) {
            return i;
        }
    }
    return -1;
}

int getInvalidEnd(VectorXcf valid_points)
{
    for (int i = valid_points.rows() - 1; i >= 0; i--) {
        if (valid_points.real()(i) == 0) {
            return i;
        }
    }
    return -1;
}

void reconstruct(VectorXcf spectrum, VectorXcf valid_points)
{
    VectorXcf hamming_320 = loadComplexVector("data/hamming_320.dat", 320);
    VectorXcf invalid_points = VectorXcf::Ones(320).array() - valid_points.array();
    VectorXcf channel_impulse_response_valid_points = VectorXcf::Ones(320);

    int invalid_start = getInvalidStart(valid_points);
    int invalid_end = getInvalidEnd(valid_points);

    float stop_threshold = 0.03;
    bool stop_left = false;
    bool stop_right = false;
    bool last_run = false;
    int left_cursor = 100;
    int right_cursor = 220;
    MatrixXcf idft = loadComplexMatrix("data/idft.dat", 320, 320);
    VectorXcf spectrum_windowed = spectrum.array() * hamming_320.array();
    VectorXcf known = extractKnown(spectrum_windowed, valid_points);

    //while (!stop_left || !stop_right || last_run) {
        //VectorXcf zeroed = buildZeroedPoints(left_cursor, right_cursor, 320);
        //VectorXcf non_zeroed = VectorXcf::Ones(320).array() - zeroed.array();


        MatrixXcf F1_1 = idft.topLeftCorner(left_cursor, invalid_start);
        MatrixXcf F1_2 = idft.topRightCorner(left_cursor, idft.cols() - invalid_end - 1);
        MatrixXcf F1_3 = idft.bottomLeftCorner(idft.rows() - right_cursor, invalid_start);
        MatrixXcf F1_4 = idft.bottomRightCorner(idft.rows() - right_cursor, idft.cols() - invalid_end - 1);
        MatrixXcf F1(F1_1.rows() + F1_3.rows(), F1_1.cols() + F1_2.cols());
        F1 << F1_1, F1_2,
           F1_3, F1_4;

        MatrixXcf F2_1 = idft.block(0, invalid_start, left_cursor, invalid_end - invalid_start);
        MatrixXcf F2_2 = idft.block(right_cursor, invalid_start, idft.rows() - right_cursor, invalid_end - invalid_start);
        MatrixXcf F2(F2_1.rows() + F2_2.rows(), F2_1.cols());
        F2 << F2_1,
           F2_2;

        MatrixXcf F3_1 = idft.block(left_cursor, 0, right_cursor - left_cursor, invalid_start);
        MatrixXcf F3_2 = idft.block(left_cursor, invalid_start, right_cursor - left_cursor, idft.cols() - invalid_end);
        MatrixXcf F3(F3_1.rows(), F3_1.cols() + F3_2.cols());
        F3 << F3_1, F3_2;

        MatrixXcf F4 = idft.block(left_cursor, invalid_start, right_cursor - left_cursor, invalid_end - invalid_start);

        VectorXcf y = -F1 * known;

        VectorXcf unknown = irls(F2, y);


        //printf("F1_1: %d,%d %dx%d\n", 0, 0, (int)F1_1.rows(), (int)F1_1.cols());
        //printf("F1_2: %d,%d %dx%d\n", 0, 0, (int)F1_2.rows(), (int)F1_2.cols());
        //printf("F1_3: %d,%d %dx%d\n", 0, 0, (int)F1_3.rows(), (int)F1_3.cols());
        //printf("F1_4: %d,%d %dx%d\n", 0, 0, (int)F1_4.rows(), (int)F1_4.cols());
    //}



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
