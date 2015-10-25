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

VectorXcf extractKnown(VectorXcf spectrum, VectorXcf valid_points)
{
    ArrayXf a = valid_points.real().array();
    VectorXcf extracted = VectorXcf::Zero((a > 0).count());

    int j = 0;
    for (int i = 0; i < a.rows(); i++) {
        if (a(i) > 0) {
            extracted(j) = spectrum(i);
            j++;
        }
    }
    return extracted;
}

VectorXcf reconstruct(VectorXcf spectrum, VectorXcf valid_points)
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
    VectorXcf unknown;

    while (!stop_left || !stop_right || last_run) {

        printf("[left_cursor, right_cursor]: [%d, %d]\n", left_cursor, right_cursor);

        MatrixXcf F1(left_cursor + idft.rows() - right_cursor, invalid_start + idft.cols() - invalid_end - 1);
        F1.topLeftCorner(left_cursor, invalid_start) = idft.topLeftCorner(left_cursor, invalid_start);
        F1.topRightCorner(left_cursor, idft.cols() - invalid_end - 1) = idft.topRightCorner(left_cursor, idft.cols() - invalid_end - 1);
        F1.bottomLeftCorner(idft.rows() - right_cursor, invalid_start) = idft.bottomLeftCorner(idft.rows() - right_cursor, invalid_start);
        F1.bottomRightCorner(idft.rows() - right_cursor, idft.cols() - invalid_end - 1) = idft.bottomRightCorner(idft.rows() - right_cursor, idft.cols() - invalid_end - 1);

        MatrixXcf F2_1 = idft.block(0, invalid_start, left_cursor, invalid_end - invalid_start + 1);
        MatrixXcf F2_2 = idft.block(right_cursor, invalid_start, idft.rows() - right_cursor, invalid_end - invalid_start + 1);
        MatrixXcf F2(F2_1.rows() + F2_2.rows(), F2_1.cols());
        F2 << F2_1,
           F2_2;

        MatrixXcf F3_1 = idft.block(left_cursor, 0, right_cursor - left_cursor, invalid_start);
        MatrixXcf F3_2 = idft.block(left_cursor, invalid_start, right_cursor - left_cursor, idft.cols() - invalid_end - 1);
        MatrixXcf F3(F3_1.rows(), F3_1.cols() + F3_2.cols());
        F3 << F3_1, F3_2;

        MatrixXcf F4 = idft.block(left_cursor, invalid_start, right_cursor - left_cursor, invalid_end - invalid_start + 1);

        //printf("F1: %dx%d\n", (int)F1.rows(), (int)F1.cols());
        //printf("F2: %dx%d\n", (int)F2.rows(), (int)F2.cols());
        //printf("F3: %dx%d\n", (int)F3.rows(), (int)F3.cols());
        //printf("F4: %dx%d\n", (int)F4.rows(), (int)F4.cols());

        VectorXcf y = -F1 * known;

        unknown = irls(F2, y);

        VectorXcf zeroed_points = (F1 * known) + (F2 * unknown);
        VectorXcf nonzeroed_points = (F3 * known) + (F4 * unknown);

        float zeroed_max = zeroed_points.array().abs().maxCoeff();
        float nonzeroed_max = nonzeroed_points.array().abs().maxCoeff();
        float threshold = stop_threshold * (zeroed_max > nonzeroed_max ? zeroed_max : nonzeroed_max);
        //cout << "threshold: " << threshold << endl;
        //cout << iter_ir_soln.array().abs() << endl;

        VectorXcf iter_ir_soln(zeroed_points.rows() + nonzeroed_points.rows());
        iter_ir_soln << zeroed_points.head(left_cursor),
                  nonzeroed_points,
                  zeroed_points.tail(idft.rows() - right_cursor);
        ArrayXf abs_iter_ir_soln = iter_ir_soln.array().abs();

        last_run = false;
        if (!stop_left) {
            ArrayXf top_rows = zeroed_points.topRows(left_cursor).array().abs();
            if ((top_rows > threshold).any()) {
                stop_left = true;
                left_cursor--;
            } else {
                int quick_search = -1;
                for (int i = 0; i < abs_iter_ir_soln.rows(); i++) {
                    if (abs_iter_ir_soln(i) > threshold) {
                        quick_search = i;
                        break;
                    }
                }
                left_cursor = (quick_search > left_cursor+1) ? quick_search : left_cursor + 1;
            }
            last_run = true;
        }
        if (!stop_right) {
            ArrayXf bottom_rows = zeroed_points.bottomRows(idft.rows() - right_cursor).array().abs();
            if ((bottom_rows > threshold).any()) {
                stop_right = true;
                right_cursor++;
            } else {
                int quick_search = -1;
                for (int i = abs_iter_ir_soln.rows() - 1; i > 0; i--) {
                    if (abs_iter_ir_soln(i) > threshold) {
                        quick_search = i;
                        break;
                    }
                }
                right_cursor = (quick_search != -1 && quick_search < right_cursor) ? quick_search : right_cursor - 1;
            }
            last_run = true;
        }
    }

    VectorXcf spectrum_recovered(known.rows() + unknown.rows());
    spectrum_recovered << known.head(invalid_start),
                       unknown,
                       known.tail(spectrum.rows() - invalid_end - 1);
    VectorXcf spectrum_recovered_unwindowed = spectrum_recovered.array() / hamming_320.array();
    return spectrum_recovered_unwindowed;
}

int main()
{
    cout << "Beginning reconstruction" << endl;
    VectorXcf spectrum = loadComplexVector("data/freq_resp.dat", 320);
    VectorXcf valid_points = loadComplexVector("data/valid_points.dat", 320);

    ofstream test_in("data/test_in.dat");
    test_in << spectrum.array() * valid_points.array();


    std::clock_t start, end;
    start = std::clock();
    VectorXcf recovered = reconstruct(spectrum, valid_points);
    end = std::clock();
    cout << "Took: " << (end - start) / (double)(CLOCKS_PER_SEC / 1000) << "ms\n" << endl;

    ofstream test_out("data/test_out.dat");
    test_out << recovered;

    cout << "done" << endl;
}
