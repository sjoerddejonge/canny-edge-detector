//
// Created by Sjoerd de Jonge on 13/05/2026.
//

#include "image_processing.h"
#include "gaussian.h"
#include "constants.h"

// TODO: Decide whether input and output should be Matrix<double>, as it does represent an image (which is uint8_t).
// TODO: Make templated Matrix<T> for input, kernel, and output if needed.
/// Convolution on a Matrix with a kernel Matrix. The input Matrix represents a picture.
Matrix<double> convolve(const Matrix<double>& input, const Matrix<double>& kernel) {
    auto output = Matrix<double>(input.getWidth(), input.getHeight(), input.getLayers());

    // Kernel should have an odd size (width AND height):
    if (kernel.getWidth() % 2 == 0 || kernel.getHeight() % 2 == 0){
        throw std::invalid_argument(
            std::format("convolve: Kernel width/height must be an odd number. Kernel width: {}. Kernel height: {}",
                kernel.getWidth(), kernel.getHeight())
        );
    }
    // Kernel cannot be a 3D Matrix:
    if (kernel.getLayers() > 1){
        throw std::invalid_argument("convolve: Unable to convolve image, kernel cannot be a 3D Matrix");
    }
    // Kernel cannot be larger than input:
    if (kernel.getWidth() > input.getWidth() || kernel.getHeight() > input.getHeight()) {
        throw std::invalid_argument("convolve: Kernel should not be larger than input image");
    }

    const bool grayscale = input.getLayers() == 1;

    // Creating constant copies of the image parameters for shorter naming:
    const int in_w = input.getWidth();   // Input Matrix width
    const int in_h = input.getHeight();  // Input Matrix height
    const int in_l = input.getLayers();  // Input Matrix layers. For color images either 3 or 4 channels (RGB or RGBA).
    const std::vector<double>& im_data = input.getData();
    // Creating constant copies of the kernel parameters for shorter naming:
    const int k_w = kernel.getWidth();  // Kernel Matrix width
    const int k_h = kernel.getHeight(); // Kernel Matrix height
    const std::vector<double>& k_data = kernel.getData();

    double acc; // Accumulator for the new pixel value. If image is in color, this represents the blue channel.
    double acc_g; // Accumulator for the green pixel value if image is in color.
    double acc_r; // Accumulator for the red pixel value if image is in color.

    // Looping through the image:
    for (int y = 0; y < in_h; ++y){
        for (int x = 0; x < in_w; ++x){
            acc = 0;
            acc_g = 0;
            acc_r = 0;

            // For each pixel in the image, this loops through the kernel elements:
            for (int r = 0; r < k_h; ++r){   // For each kernel row
                for (int c = 0; c < k_w; ++c){   // For each kernel column
                    // The center (origin) of the kernel is placed on the current image pixel. The kernel will overlap
                    // the neighbouring pixels. The X and Y coordinates of these pixels are calculated first.
                    // In c++ positive integers are always truncated towards 0, so using floor() is not required here.
                    int posX = x-k_w/2+c;   // posX = 'position X', the X-coordinate on the image which overlaps with the current point (c,r) on the kernel.
                    int posY = y-k_h/2+r;  // posY = 'position Y', Y-coordinate on the image that is overlapped with (c,r) of the kernel
                    // Here it is verified whether the coordinates lie within the image borders, if not they are
                    // mirrored into the image:
                    if ( (0 > posX) || (posX >= in_w) || (0 > posY) || (posY >= in_h) ){
                        // If negative value, mirror the coordinates into the image:
                        posX = abs(posX);
                        posY = abs(posY);
                        // If coordinates are larger than image width, mirror coordinates into image:
                        posX = (in_w-1) - abs(posX - (in_w-1));
                        posY = (in_h-1) - abs(posY - (in_h-1));
                        //... and no extra if-statements are needed! :)
                        // Best way to visualize this is to calculate an example, say im_width = 240, im_height = 240
                        // and posX = 245 (outside the image): (240-1) - abs(245 - (240-1)) = 239 - 6 = 233, which is
                        // as close to 239 as 245 is.
                        // But if posX = 233 (inside the image), we get: (240-1) - abs(233 - (240-1)) = 239 - 6 = 233,
                        // so nothing changes.
                        // Btw, 'im_width-1' (etc) is used, because, for example, a width of 240 pixels means the maximum
                        // image coordinate is 239, because coordinates starts at 0 (but width/height starts at 1).
                    }
                    // The value of the overlapping kernel and image are multiplied and added to the accumulator.
                    acc += im_data[in_l * (posY * in_w + posX) + 0] * k_data[r*k_w + c];
                    // If image is not grayscale, green and red channels need to be calculated too:
                    if (!grayscale) {
                        acc_g += im_data[in_l * (posY * in_w + posX) + 1] * k_data[r * k_w + c];
                        acc_r += im_data[in_l * (posY * in_w + posX) + 2] * k_data[r * k_w + c];
                    }
                }
            }
            // Write the accumulated value to the BGR channels of the pixels of the output image:
            //output.data[channels * (y * im_width + x) + 0] = acc; <- this was the old 'public-variable' version
            output.setData(in_l * (y * in_w + x) + 0, acc);
            //std::cout << acc << std::endl;
            if(grayscale) {
                output.setData(in_l * (y * in_w + x) + 1, acc);
                output.setData(in_l * (y * in_w + x) + 2, acc);
            } else {
                output.setData(in_l * (y * in_w + x) + 1, acc_g);
                output.setData(in_l * (y * in_w + x) + 2, acc_r);
            }

        }
    }

    return output;
}

/**
 * Matrix<double> gaussianBlur( ... )
 * Constructs a Gaussian kernel and uses the convolve() function for convolution using that kernel.
 * Kernel types are:
 * 0: Normal Gaussian blur
 * 1: Gaussian blur with derivative with respect to X
 * 2: Gaussian blur with derivative with respect to Y
 * 3: Gaussian blur with derivative with respect to XY
 */
Matrix<double> gaussianBlur(const Matrix<double>& input, double sigma, int kernel_type) {
    if (kernel_type < 0 || kernel_type > 3){
        throw std::invalid_argument("gaussianBlur: invalid kernel_type");
    }
    int size = 2* static_cast<int>(ceil(sigma * 4))+1;   // Calculate the kernel size based on the sigma
    Matrix<double> g_hor(size, 1, 1);      // Empty horizontal kernel matrix of 1 x size
    Matrix<double> g_ver(1, size, 1);      // Empty vertical  kernel matrix of size x 1
    // Constructing Gaussian kernels:
    switch (kernel_type){
        case 0:     // Normal Gaussian blur
            g_hor = constructGaussianKernel(g_hor.getWidth(), g_hor.getHeight(), sigma);
            g_ver = constructGaussianKernel(g_ver.getWidth(), g_ver.getHeight(), sigma);
            break;
        case 1:     // Derivative of X
            g_hor = constructGaussianKernelDerivative(g_hor.getWidth(), g_hor.getHeight(), sigma);
            g_ver = constructGaussianKernel((g_ver.getWidth()), g_ver.getHeight(), sigma);
            break;
        case 2:     // Derivative of Y
            g_hor = constructGaussianKernel(g_hor.getWidth(), g_hor.getHeight(), sigma);
            g_ver = constructGaussianKernelDerivative((g_ver.getWidth()), g_ver.getHeight(), sigma);
            break;
        case 3:     // Derivative of XY
            g_hor = constructGaussianKernelDerivative(g_hor.getWidth(), g_hor.getHeight(), sigma);
            g_ver = constructGaussianKernelDerivative((g_ver.getWidth()), g_ver.getHeight(), sigma);
            break;
        default:    // Normal Gaussian blur
            g_hor = constructGaussianKernel(g_hor.getWidth(), g_hor.getHeight(), sigma);
            g_ver = constructGaussianKernel(g_ver.getWidth(), g_ver.getHeight(), sigma);
            break;
    }
    // Convolve our image with the two Gaussian matrices, resulting in a new matrix:
    Matrix<double> output = convolve(input, g_hor);
    output = convolve(output, g_ver);
    return output;
}

/**
 * Matrix<double> gradientMagnitude( ... )
 * Calculate the gradient magnitude using the X and Y derivatives of an image. The X and Y derivatives are
 * obtained by convolving an image with the derivative of a Gaussian with respect to X or Y.
 * The formula for calculating the gradient magnitude is:
 * Gradient magnitude = sqrt(derivX^2 + derivY^2)
 */
Matrix<double> gradientMagnitude(const Matrix<double>& im_derivX, const Matrix<double>& im_derivY) {
    const int width = im_derivX.getWidth();
    const int height = im_derivX.getHeight();
    const int channels = im_derivX.getLayers();
    Matrix<double> output = im_derivX;   // Create a copy of im_derivX that will be the output

    // Check if both input matrices are equal in dimensions:
    if (width != im_derivY.getWidth() || height != im_derivY.getHeight()) {
        throw std::invalid_argument("gradientMagnitude: to calculate gradient magnitude both input images should have similar size");
    }

    // Loop through all the pixels of the output image/matrix:
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Index for the 'blue' value of the image, of the current pixel:
            int b_index = (channels*(y * width+ x)+0);  // image is grayscale so green and red values are the same

            // Gradient magnitude value for that pixel, using the pixels of im_derivX and im_derivY:
            // The formula is: Magnitude = sqrt( derivativeX^2 + derivativeY^2 )
            const double grad_value = sqrt(pow(im_derivX.getData()[b_index],2) + pow(im_derivY.getData()[b_index],2));

            // Set the blue pixel value for the output matrix as the calculated gradient:
            output.setData(b_index,grad_value);

            // Green and red channel are same as blue channel because image is in grayscale:
            output.setData(b_index+1,grad_value); // green channel = same as blue channel
            output.setData(b_index+2,grad_value); // red channel = same as blue channel

        }
    }
    return output;
}

/**
 * Matrix<double> nonMaximumSuppression( ... )
 * Compares the value of the current pixel in the gradient magnitude matrix with those in the same direction of the
 * gradient direction and suppresses all values that are not the maximum. The gradient direction for a pixel can be
 * calculated as:
 * direction = arctan( derivY / derivX )
 * Direction will be rounded to one of four angles: horizontal (0 or 180), vertical (90 or 270) and both diagonals
 * (45 or 225 && 135 or 315). For each pixel in the gradient magnitude matrix, the two neighbouring pixels in that
 * direction will be compared. If neither of those pixels have a value larger than the current pixel, the value of
 * the current pixel is kept. Otherwise it is thrown away.
 */
Matrix<double> nonMaximumSuppression(const Matrix<double>& gradMag, const Matrix<double>& im_derivX, const Matrix<double>& im_derivY){
    Matrix<double> output = gradMag;
    const int channels = gradMag.getLayers();
    const int height = gradMag.getHeight();
    const int width = gradMag.getWidth();

    for (int y = 0; y < height; ++y){
        for (int x = 0; x < width; ++x){
            const int index = (channels*(y * width+ x)+0);
            // Formula for direction: direction = arctan(derivY/derivX)
            double direction = std::atan(im_derivY.getData()[index]/im_derivX.getData()[index]);
            direction = direction * 180/PI; // Convert direction to degrees
            if (direction < 0){ // if direction is negative degrees
                direction = 360 + direction; // find positive equivalent
            }

            // Rounding the direction to be 0, 45, 90, 135, 180, 225, 270 or 315:
            // Creating an array with absolute difference between the direction and the given angles:
            double difference[8] { abs(direction - 0), abs(direction - 45), abs(direction - 90),
                                   abs(direction - 135), abs(direction - 180), abs(direction - 225),
                                   abs(direction - 270), abs(direction - 315) };
            // The smallest difference then gives which angle is closest to the direction, which will be used onwards.
            // Here is a list of which array indices correspond with which angle:
            // index 0  =   angle 0
            // index 1  =   angle 45
            // index 2  =   angle 90
            // index 3  =   angle 135
            // index 4  =   angle 180
            // index 5  =   angle 225
            // index 6  =   angle 270
            // index 7  =   angle 315
            // Pre-determined angles:
            constexpr int angles[8] { 0, 45, 90, 135, 180, 225, 270, 315};
            // Finding the index of the smallest difference:
            const int min_element_index = static_cast<int>(std::distance(&difference[0], std::min_element(difference,difference+8)));
            // Set the direction as one of the pre-determined angles:
            direction = angles[min_element_index];

            // Horizontal direction:
            if (direction == 0 || direction == 180){
                // Look at pixels to the left and right:
                if (x-1 >= 0){
                    // If the pixel to the left has a larger value, suppress the value in the current pixel:
                    if (gradMag.getData()[(channels*(y * width+ x)+0)] < gradMag.getData()[(channels*(y * width+ x-1)+0)]){
                        output.setData((channels*(y * width+ x)+0), 0); // B value
                        output.setData((channels*(y * width+ x)+1), 0); // G value
                        output.setData((channels*(y * width+ x)+2), 0); // R value
                    }
                }
                if (x+1 < width){
                    // If the pixel to the right has a larger value, suppress the value in the current pixel:
                    if (gradMag.getData()[(channels*(y * width+ x)+0)] < gradMag.getData()[(channels*(y * width+ x+1)+0)]){
                        output.setData((channels*(y * width+ x)+0), 0); // B value
                        output.setData((channels*(y * width+ x)+1), 0); // G value
                        output.setData((channels*(y * width+ x)+2), 0); // R value
                    }
                }
            }

            // Diagonal direction towards up-left
            if (direction == 45 || direction == 225){
                // If the pixel diagonally up-left exists:
                if (x-1 >=0 && y-1 >= 0){
                    // If the pixel diagonally up-left has a larger value, suppress the value in the current pixel:
                    if (gradMag.getData()[(channels*(y * width+ x)+0)] < gradMag.getData()[(channels*((y-1) * width+ x-1)+0)]){
                        output.setData((channels*(y * width+ x)+0), 0); // B value
                        output.setData((channels*(y * width+ x)+1), 0); // G value
                        output.setData((channels*(y * width+ x)+2), 0); // R value
                    }
                }
                // If the pixel diagonally down-right exists:
                if (x+1 < width && y+1 < height ){
                    // If the pixel diagonally down-right has a larger value, suppress the value in the current pixel:
                    if (gradMag.getData()[(channels*(y * width+ x)+0)] < gradMag.getData()[(channels*((y+1) * width+ x+1)+0)]){
                        output.setData((channels*(y * width+ x)+0), 0); // B value
                        output.setData((channels*(y * width+ x)+1), 0); // G value
                        output.setData((channels*(y * width+ x)+2), 0); // R value
                    }
                }
            }

            // Vertical direction:
            if (direction == 90 || direction == 270){
                // If the pixel above exists:
                if (y-1 >= 0){
                    // If the pixel above has a larger value, suppress the value in the current pixel:
                    if (gradMag.getData()[(channels*(y * width+ x)+0)] < gradMag.getData()[(channels*((y-1) * width+ x)+0)]){
                        output.setData((channels*(y * width+ x)+0), 0); // B value
                        output.setData((channels*(y * width+ x)+1), 0); // G value
                        output.setData((channels*(y * width+ x)+2), 0); // R value
                    }
                }
                // If the pixel below exists:
                if (y+1 < height){
                    // If the pixel below has a larger value, suppress the value in the current pixel:
                    if (gradMag.getData()[(channels*(y * width+ x)+0)] < gradMag.getData()[(channels*((y+1) * width+ x)+0)]){
                        output.setData((channels*(y * width+ x)+0), 0); // B value
                        output.setData((channels*(y * width+ x)+1), 0); // G value
                        output.setData((channels*(y * width+ x)+2), 0); // R value
                    }
                }
            }

            // Diagonal direction towards up-right:
            if (direction == 135 || direction == 315){
                // If the pixel diagonally up-right exists:
                if (x+1 < width && y-1 >= 0){
                    // If the pixel diagonally up-right has a larger value, suppress the value in the current pixel:
                    if (gradMag.getData()[(channels*(y * width+ x)+0)] < gradMag.getData()[(channels*((y-1) * width+ x+1)+0)]){
                        output.setData((channels*(y * width+ x)+0), 0); // B value
                        output.setData((channels*(y * width+ x)+1), 0); // G value
                        output.setData((channels*(y * width+ x)+2), 0); // R value
                    }
                }
                // If the pixel diagonally down-left exists:
                if (x-1 >= 0 && y+1 < height ){
                    // If the pixel diagonally down-left has a larger value, suppress the value in the current pixel:
                    if (gradMag.getData()[(channels*(y * width+ x)+0)] < gradMag.getData()[(channels*((y+1) * width+ x-1)+0)]){
                        output.setData((channels*(y * width+ x)+0), 0); // B value
                        output.setData((channels*(y * width+ x)+1), 0); // G value
                        output.setData((channels*(y * width+ x)+2), 0); // R value
                    }
                }
            }
        }
    }
    return output;

}

// TODO: Change Matrix<double>& input to Matrix<T>?
double getImageAverage(const Matrix<double>& input) {
    double average = 0;
    // Loop through the image to compute the sum of all elements (stored inside the double 'average'):
    for (const double i : input.getData()){
        average += i;
    }
    average /= static_cast<double>(input.getData().size()); // Calculate the average by dividing the sum by the amount of elements in the image
    return average;
}

/**
 * void hysteresisThresholding( ... )
 * gradMag = Gradient magnitude matrix/image
 * high_threshold = High threshold
 * low_threshold = Low threshold
 * Set two thresholds, a high threshold and a low threshold.
 *      If a gradient pixel is above the high threshold, mark it as an edge pixel.
 *      If a gradient pixel is below the low threshold, mark it as a non-edge pixel.
 *      If a gradient pixel is between the low and high threshold, mark it as an edge pixel if it is connected
 *      to a edge pixel directly, or indirectly via other low-high pixels.
 */
// TODO: Should return a Matrix<double> like the other functions.
void hysteresisThresholding(Matrix<double> &gradMag, const double high_threshold, const double low_threshold){
    const int channels = gradMag.getLayers();
    const int height = gradMag.getHeight();
    const int width = gradMag.getWidth();
    int count_edge = 0; // Counts the primary edge pixels (above high threshold), useful info for debugging
    int count_edge_sec = 0; // Counts the secondary edge pixels (between low and high threshold)

    // First loop, gather all edge pixels:
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // If pixel value is above HIGH threshold, set its value to max (edge pixel)
            if (gradMag.getData()[channels*(y * width+x)+0] >= high_threshold){
                gradMag.setData(channels*(y * width+x)+0, 255);
                gradMag.setData(channels*(y * width+x)+1, 255);
                gradMag.setData(channels*(y * width+x)+2, 255);
                count_edge++;

                // Then look at neighbour pixels, if those are above LOW threshold, set their values to max as well.
                // The neighbours of pixel 'X' are:
                //   ___________
                //  | 1 | 2 | 3 |
                //  | 4 | X | 5 |
                //  | 6 | 7 | 8 |
                //   ‾‾‾‾‾‾‾‾‾‾‾
                // The neighbouring pixels are checked in this order.
                // If 1, 2, 3 or 4 is added to the edge pixels, the loop will go back to that pixel before advancing.
                // This is to make sure no additional edge pixels are missed.

                // 1. Pixel above-left:
                if (x-1 >= 0){
                    if (y-1 >= 0){
                        if ((gradMag.getData()[channels*((y-1) * width+x-1)+0] >= low_threshold) &&
                            (gradMag.getData()[channels*((y-1) * width+x-1)+0] < high_threshold) ){
                            gradMag.setData(channels*((y-1) * width+x-1)+0, 255);
                            gradMag.setData(channels*((y-1) * width+x-1)+1, 255);
                            gradMag.setData(channels*((y-1) * width+x-1)+2, 255);
                            count_edge_sec++;
                            // The next pixel that will be looked at in the loop is (x-1, y-1):
                            x = (x-1)-1;    // x is decreased by 2 because as soons as the loop ends, it will increase by 1
                            y = y-1;
                            continue;
                        }
                    }
                }

                // 2. Pixel above:
                if (y-1 >= 0){
                    if ((gradMag.getData()[channels*((y-1) * width+x)+0] >= low_threshold) &&
                        (gradMag.getData()[channels*((y-1) * width+x)+0] < high_threshold)){
                        gradMag.setData(channels*((y-1) * width+x)+0, 255);
                        gradMag.setData(channels*((y-1) * width+x)+1, 255);
                        gradMag.setData(channels*((y-1) * width+x)+2, 255);
                        count_edge_sec++;
                        // The next pixel that will be looked at in the loop is (x, y-1):
                        x = (x)-1;
                        y = y-1;
                        continue;
                    }
                }

                // 3. Pixel above-right:
                if (x+1 < width){
                    if (y-1 >= 0){
                        if ((gradMag.getData()[channels*((y-1) * width+x+1)+0] >= low_threshold) &&
                            (gradMag.getData()[channels*((y-1) * width+x+1)+0] < high_threshold)){
                            gradMag.setData(channels*((y-1) * width+x+1)+0, 255);
                            gradMag.setData(channels*((y-1) * width+x+1)+1, 255);
                            gradMag.setData(channels*((y-1) * width+x+1)+2, 255);
                            count_edge_sec++;
                            // The next pixel that will be looked at in the loop is (x+1, y-1):
                            x = (x+1)-1;
                            y = y-1;
                            continue;
                        }
                    }
                }

                // 4. Pixel to the left:
                if (x-1 >= 0){
                    if ((gradMag.getData()[channels*(y * width+x-1)+0] >= low_threshold) &&
                        (gradMag.getData()[channels*(y * width+x-1)+0] < high_threshold)){
                        gradMag.setData(channels*(y * width+x-1)+0, 255);
                        gradMag.setData(channels*(y * width+x-1)+1, 255);
                        gradMag.setData(channels*(y * width+x-1)+2, 255);
                        count_edge_sec++;
                        // The next pixel that will be looked at in the loop is (x-1, y):
                        x = (x-1)-1;
                        // y = y;
                        continue;
                    }
                }

                // 5. Pixel to the right:
                if (x+1 < width){
                    if (gradMag.getData()[channels*(y * width+x+1)+0] >= low_threshold){
                        gradMag.setData(channels*(y * width+x+1)+0, 255);
                        gradMag.setData(channels*(y * width+x+1)+1, 255);
                        gradMag.setData(channels*(y * width+x+1)+2, 255);
                        count_edge_sec++;
                    }
                }

                // 6. Pixel below-left:
                if (x-1 >= 0){
                    if (y+1 < height){
                        if (gradMag.getData()[channels*((y+1) * width+x-1)+0] >= low_threshold){
                            gradMag.setData(channels*((y+1) * width+x-1)+0, 255);
                            gradMag.setData(channels*((y+1) * width+x-1)+1, 255);
                            gradMag.setData(channels*((y+1) * width+x-1)+2, 255);
                            count_edge_sec++;
                        }
                    }
                }

                // 7. Pixel below:
                if (y+1 < height){
                    if (gradMag.getData()[channels*((y+1) * width+x)+0] >= low_threshold){
                        gradMag.setData(channels*((y+1) * width+x)+0, 255);
                        gradMag.setData(channels*((y+1) * width+x)+1, 255);
                        gradMag.setData(channels*((y+1) * width+x)+2, 255);
                        count_edge_sec++;
                    }
                }

                // 8. Pixel below-right:
                if (x+1 < width){
                    if (y+1 < height){
                        if (gradMag.getData()[channels*((y+1) * width+x+1)+0] >= low_threshold){
                            gradMag.setData(channels*((y+1) * width+x+1)+0, 255);
                            gradMag.setData(channels*((y+1) * width+x+1)+1, 255);
                            gradMag.setData(channels*((y+1) * width+x+1)+2, 255);
                            count_edge_sec++;
                        }
                    }
                }
            }
            // If pixel value is below LOW threshold, discard it
            else if (gradMag.getData()[channels*(y * width+ x)+0] < low_threshold){
                gradMag.setData(channels*(y * width+ x)+0, 0);
                gradMag.setData(channels*(y * width+ x)+1, 0);
                gradMag.setData(channels*(y * width+ x)+2, 0);
            }
        }
    }
    // Second loop, remove all pixels between low and high that are not connected:
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (gradMag.getData()[channels*(y * width+x)+0] < high_threshold){
                gradMag.setData(channels*(y * width+x)+0, 0);
                gradMag.setData(channels*(y * width+x)+1, 0);
                gradMag.setData(channels*(y * width+x)+2, 0);
            }
        }
    }
    //std::cout << "Primary edge pixels: " << count_edge << std::endl;
    //std:: cout << "Secondary edge pixels: " << count_edge_sec << std::endl;
    // Function is of type void so no return.
}

