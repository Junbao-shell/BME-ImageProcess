///////////////////////////////////////////////////////////
/// @copyright copyright description
/// 
/// @brief CT water model uniformity test
/// 
/// @file water_model_uniformity.cpp
/// 
/// @author GaoJunbao(junbaogao@foxmail.com)
/// 
/// @date 2022-06-29
///////////////////////////////////////////////////////////

// Current Cpp header
// System header
// C/C++ standard library header
#include <iostream>
#include <cmath>
// External library header
// Current module header
// Root directory header

constexpr float voxel_size = 0.33; // unit mm; 
constexpr int boundary_offset_distance = 10; // unit mm

struct Parameters
{
    float VoxelSizeX;
    float VoxelSizeY;
    float VoxelSizeZ;
};

struct Point
{
    int x;
    int y;
};

/**
 * @brief statistic the ROI value;
 * @note this funciton assumed that the project in the FOV Center
 * 
 * @param image image data
 * @param fov image matrix size(one side length)
 * @param project_diameter project size, e.g. 200 mm water model, 300 mm water model
 */
void CalculateROI(const float *image, const int fov, const int project_diameter)
{
    const int project_voxel_number_in_diameter = project_diameter / voxel_size;  // 水模直径
    const int roi_radius = std::ceil(project_voxel_number_in_diameter / 10 / 2.0); // ROI radius, unit voxel count
    const int boundary_offset_voxel = boundary_offset_distance / voxel_size;     // roi to project distance, unit voxel count
    
    const int center_to_circle_size = (project_voxel_number_in_diameter / 2) - boundary_offset_voxel - roi_radius;

    // center point
    Point center {std::ceil(fov / 2.0), std::ceil(fov / 2.0)};
    Point north_point {center.x, center.y - center_to_circle_size};
    Point east_point {center.x + center_to_circle_size, center.y};
    Point south_point {center.x, center.y + center_to_circle_size};
    Point west_point {center.x - center_to_circle_size, center.y};

}

float GetRoiMeans(float *image, const int image_size, const Point point, const int radius)
{
    for (int i = 0; i < image_size; ++i)
    {
        for (int j = 0; j < image_size; ++j)
        {
            if ((std::pow((i - point.x), 2) + std::pow((j - point.y), 2)) < radius * radius)
            {
                image[j * image_size + i] = 2.0f;
            }
        }
    }
}

int main(int argc, char **argv)
{
    const int voxel_num_x = 512;
    const int image_size = voxel_num_x * voxel_num_x;
    float *image = new float[image_size]();

    Point center {256, 256};
    Point roi {128, 256};

    GetRoiMeans(image, voxel_num_x, roi, 20);

    delete[] image;
    image = nullptr;

    return 0;
}