//
// Created by Sjoerd de Jonge on 13/05/2026.
//

#ifndef CANNYEDGEDETECTOR_IMAGE_PROCESSING_H
#define CANNYEDGEDETECTOR_IMAGE_PROCESSING_H
#include "Matrix.h"

/**
 * Image processing is done on the Matrix<double> class for extended precision (the Image class has a vector<uint8_t>).
 */

/**
 * @brief       Calculate the average value of all elements.
 * @param       input: The Matrix<double> to get the average of.
 * @return      The average 'pixel' value.
 */
double getImageAverage(const Matrix<double>& input);

/**
 * @brief       Convolution on an image/matrix with a kernel.
 * @param       input: The input Matrix<double> to perform the convolution on.
 * @param       kernel: The kernel Matrix<double> to perform the convolution with.
 * @return      The resulting Matrix<double> after the convolution.
 * @see         https://en.wikipedia.org/wiki/Kernel_(image_processing)#Convolution
 */
Matrix<double> convolve(const Matrix<double>& input, const Matrix<double>& kernel); // Convolution

/**
 * @brief       Blur an image (represented as Matrix<double>) by convolution with a Gaussian kernel.
 * @param       input: The input Matrix<double> to blur.
 * @param       sigma: The size of the Gaussian kernel.
 * @param       kernel_type: Whether to use a normal Gaussian kernel (0), derivative with respect to X (1),
 *              derivative with respect to Y (2), or derivative with respect to XY (3).
 * @return      The blurred Matrix<double>.
 */
Matrix<double> gaussianBlur(const Matrix<double>& input, double sigma, int kernel_type);

/**
 * @brief       Calculate the gradient magnitude using the X and Y derivatives of an image. The X and Y derivatives are
 *              obtained by convolving an image with the derivative of a Gaussian with respect to X or Y.
 *              The formula for calculating the gradient magnitude is: \n
 *              Gradient magnitude = sqrt(derivX^2 + derivY^2)
 * @param       im_derivX: Image convolved with a derivative of a Gaussian kernel with respect to X.
 * @param       im_derivY: Image convolved with a derivative of a Gaussian kernel with respect to Y.
 * @return      The Matrix<double> representing the gradient magnitude.
 */
Matrix<double> gradientMagnitude(const Matrix<double>& im_derivX, const Matrix<double>& im_derivY);

/**
 * @brief       Compares the value of the current pixel in the gradient magnitude matrix with those in the same direction of the gradient direction and suppresses all values that are not the maximum.
 *
 * The gradient direction for a pixel can be calculated as: \n
 * direction = arctan( derivY / derivX )
 *
 * Direction will be rounded to one of four angles: horizontal (0 or 180), vertical (90 or 270) and both diagonals
 * (45 or 225 && 135 or 315). For each pixel in the gradient magnitude matrix, the two neighbouring pixels in that
 * direction will be compared. If neither of those pixels
 * have a value larger than the current pixel, the value of the current pixel is kept. Otherwise it is thrown away.
 * @param       gradMag: The gradient magnitude image/matrix.
 * @param       im_derivX: Image convolved with a derivative of a Gaussian kernel with respect to X.
 * @param       im_derivY: Image convolved with a derivative of a Gaussian kernel with respect to Y.
 * @return      A Matrix<double> with pixel-thin lines.
 */
Matrix<double> nonMaximumSuppression(const Matrix<double>& gradMag, const Matrix<double>& im_derivX, const Matrix<double>& im_derivY);

/**
 * @brief       Set two thresholds, a high threshold and a low threshold.
 * 1. If a gradient pixel is above the high threshold, mark it as an edge pixel.
 * 2. If a gradient pixel is below the low threshold, mark it as a non-edge pixel.
 * 3. If a gradient pixel is between the low and high threshold, mark it as an edge pixel if it is connected
 * to an edge pixel directly, or indirectly via other low-high pixels.
 * @param       gradMag: Gradient magnitude matrix/image.
 * @param       high_threshold: High threshold.
 * @param       low_threshold: Low threshold.
 */
void hysteresisThresholding(Matrix<double> &gradMag, double high_threshold, double low_threshold);

#endif //CANNYEDGEDETECTOR_IMAGE_PROCESSING_H
