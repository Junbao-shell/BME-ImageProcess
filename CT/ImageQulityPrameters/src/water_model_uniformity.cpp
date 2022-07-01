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
#include <fstream>
#include <iomanip>
#include <experimental/filesystem>
#include <cstring>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <utility>
// External library header
// Current module header
// Root directory header

constexpr int ratio_roi_project = 10;
constexpr int boundary_offset_distance = 10; // unit mm
constexpr float voxel_size = 0.33; // unit mm; 

struct Point
{
    int x;
    int y;
};

void PrintROI(const char *position, const Point point, const float means, const float variance)
{
    std::cout << std::setw(16) << (std::string(position) + std::string(" roi")) 
        << "(" << point.x << ", " << point.y <<  ")" << ": " 
        << std::fixed << std::setprecision(2) 
        << "means: "<< std::setw(8) << means << ", variance: " << std::setw(8) << variance << std::endl;
}

float SetRoiValue(float *image, const int image_size, const Point point, const int radius, const float value)
{
    for (int i = 0; i < image_size; ++i)
    {
        for (int j = 0; j < image_size; ++j)
        {
            if ((std::pow((i - point.x), 2) + std::pow((j - point.y), 2)) < radius * radius)
            {
                image[j * image_size + i] = value;
            }
        }
    }
}

/**
 * @brief Get the Roi Statistic object
 * 
 * @param image target image to statistic
 * @param image_size image size one side
 * @param point point corrdinate in 2D
 * @param radius roi radius
 * @return auto return pair struct 
 */
auto GetRoiStatistic(const float *image, const int image_size, const Point point, const int radius)
{
    std::vector<float> roi {};
    for (int i = 0; i < image_size; ++i)
    {
        for (int j = 0; j < image_size; ++j)
        {
            if ((std::pow((i - point.x), 2) + std::pow((j - point.y), 2)) < radius * radius)
            {
                roi.emplace_back(image[j * image_size + i]);
            }
        }
    }

    // caculate mean 
    auto roi_sum = std::accumulate(std::begin(roi), std::end(roi), 0.0);
    float roi_mean = roi_sum / roi.size();

    // caculate variance
    float roi_variance = 0.f;
    for (auto i : roi)
    {
        roi_variance += std::pow(i - roi_mean, 2);
    }
    roi_variance /= roi.size();
    auto roi_std_variance = std::sqrt(roi_variance);

    return std::make_pair(roi_mean, roi_std_variance);
}

std::vector<std::string> GetFileList(const std::string &path)
{
    namespace fs = std::experimental::filesystem;
    if (!fs::exists(path))
    {
        std::cout << "image not exist" << std::endl;
        return std::vector<std::string>{};
    }

    std::vector<std::string> file_list{};

    fs::directory_iterator end_iterator;
    for (fs::directory_iterator iter(path); iter != end_iterator; ++iter)
    {
        if (!fs::is_regular_file(iter->status()) || iter->path().filename().extension() != ".raw")
        {
            continue;
        }
        file_list.emplace_back(iter->path());
    }
    std::sort(file_list.begin(), file_list.end());

    return file_list;
}

std::vector<std::string> GetDirList(const std::string &path)
{
    namespace fs = std::experimental::filesystem;
    if (!fs::exists(path))
    {
        std::cout << "image not exist" << std::endl;
        return std::vector<std::string>{};
    }

    std::vector<std::string> file_list{};

    fs::directory_iterator end_iterator;
    for (fs::directory_iterator iter(path); iter != end_iterator; ++iter)
    {
        if (!fs::is_directory(iter->status()))
        {
            continue;
        }
        file_list.emplace_back(iter->path());
    }
    std::sort(file_list.begin(), file_list.end());

    return file_list;
}

auto GetCoef(const std::string &str)
{
    int pos1 = str.find('_');
    std::string substr = str.substr(pos1 + 1);
    int pos2 = substr.find('_');
    const float bhw = std::stof(substr.substr(0, pos2));
    
    substr = substr.substr(pos2 + 1);
    pos1 = substr.find('_');
    substr = substr.substr(pos1 + 1);
    pos2 = substr.find('_');
    const float asca = std::stof(substr.substr(0, pos2));

    substr = substr.substr(pos2 + 1);
    pos1 = substr.find('_');
    substr = substr.substr(pos1 + 1);
    pos2 = substr.find('_');
    const float gamma = std::stof(substr.substr(0, pos2));

    substr = substr.substr(pos2 + 1);
    pos1 = substr.find('_');
    substr = substr.substr(pos1 + 1);
    pos2 = substr.find('_');
    const float lambda = std::stof(substr.substr(0, pos2));

    return std::make_tuple(bhw, asca, gamma, lambda);
}

/**
 * @brief statistic the ROI value;
 * @note this funciton assumed that the project in the FOV Center
 * 
 * @param image image data
 * @param fov image matrix size(one side length)
 * @param project_diameter project size, e.g. 200 mm water model, 300 mm water model
 */
auto CalculateROI(float *image, const int fov, const int project_diameter)
{
    const int project_voxel_number_in_diameter = project_diameter / voxel_size;  // water model diameter
    const int roi_radius = std::ceil(project_voxel_number_in_diameter / ratio_roi_project / 2); // ROI radius, unit voxel count

    const int boundary_offset_voxel = boundary_offset_distance / voxel_size; // roi to project distance, unit voxel count
    const int center_to_circle_size = std::ceil(project_voxel_number_in_diameter / 2) - boundary_offset_voxel - roi_radius;

    // center point
    Point center {std::ceil(fov / 2), std::ceil(fov / 2)};
    Point north_point {center.x, center.y - center_to_circle_size};
    Point east_point {center.x + center_to_circle_size, center.y};
    Point south_point {center.x, center.y + center_to_circle_size};
    Point west_point {center.x - center_to_circle_size, center.y};

    std::vector<std::pair<float, float>> res {};
    auto [center_mean, center_vari] = GetRoiStatistic(image, fov, center, roi_radius);
    PrintROI("center", center, center_mean, center_vari);
    res.emplace_back(std::make_pair(center_mean, center_vari));
    
    auto [north_mean, north_vari] = GetRoiStatistic(image, fov, north_point, roi_radius);
    PrintROI("north", north_point, north_mean, north_vari);
    res.emplace_back(std::make_pair(north_mean, north_vari));

    auto [east_mean, east_vari] = GetRoiStatistic(image, fov, east_point, roi_radius);
    PrintROI("east", east_point, east_mean, east_vari);
    res.emplace_back(std::make_pair(east_mean, east_vari));
    
    auto [south_mean, south_vari] = GetRoiStatistic(image, fov, south_point, roi_radius);
    PrintROI("south", south_point, south_mean, south_vari);
    res.emplace_back(std::make_pair(south_mean, south_vari));
    
    auto [west_mean, west_vari] = GetRoiStatistic(image, fov, west_point, roi_radius);
    PrintROI("west", west_point, west_mean, west_vari);
    res.emplace_back(std::make_pair(west_mean, west_vari));

    // horizontal line
    const int horizontal_num = std::floor((center.x - west_point.x - 2 * roi_radius) / (2 * roi_radius));
    const int horizontal_gap = std::floor((center.x - west_point.x - 2 * roi_radius - horizontal_num * roi_radius * 2) / (horizontal_num + 1));

    for (int i = 0; i < horizontal_num * 2; ++i)
    {
        Point horizontal_point{0, center.y};
        if (i < horizontal_num)
        {
            horizontal_point.x = center.x - (horizontal_num - i) * roi_radius * 2 - (horizontal_num - i) * horizontal_gap;
        }
        else
        {
            horizontal_point.x = center.x + (i - horizontal_num + 1) * roi_radius * 2 + (i - horizontal_num + 1) * horizontal_gap;
        }
    
        auto [horizontal_mean, horizontal_vari] = GetRoiStatistic(image, fov, horizontal_point, roi_radius);
        auto point_name = "horizontal " + std::to_string(i + 1);
        PrintROI(point_name.c_str(), horizontal_point, horizontal_mean, horizontal_vari);
        res.emplace_back(std::make_pair(horizontal_mean, horizontal_vari));
    }

    return res;
}

void Binning(const float **image, const int binning_number, const int image_size, float *mean_image)
{
    for (int i = 0; i < binning_number; ++i)
    {
        for (int v = 0; v < image_size; ++v)
        {
            mean_image[v] += image[i][v];
        }
    }

    for (int v = 0; v < image_size; ++v)
    {
        mean_image[v] /= binning_number;
    }
}

void CT_ImageUniform(const std::string path)
{
    // image parameters
    const int voxel_num_x = 1024;
    const int image_size = voxel_num_x * voxel_num_x;
    float *image = new float[image_size]();
    float *image_bin = new float[image_size]();

    FILE *fp;
    std::fstream out;
    out.open("./ct_water_300_uniform_log002.xls", std::ios::out | std::ios::app);
    out << "index" << "\t" << "center" << "\t" << "north" << "\t" << "sorth" << "\t" << "east" << "\t" << "west" << "\t" 
        << "point1" << "\t" << "point2" << "\t" << "point3" << "\t" << "point4" << "\t" << "point5" << "\t" << "point6" << "\t" 
        << "max means" << "\t" << "max variance" << "\t\t" << "bhw" << "\t" << "asca" << "\t" << "gamma" << "\t" << "lamma" 
        << std::endl;

    int image_index = 0;
    auto dir_list = GetDirList(path);
    for (auto sub_dir : dir_list)
    {
        auto file_list = GetFileList(sub_dir);
        for (int i = 0; i < file_list.size(); ++i)
        {
            if (i > 66 && i < 79)
            {
                fp = fopen(file_list.at(i).c_str(), "rb");
                fread(image, sizeof(float), image_size, fp);
                fclose(fp);

                for (int v = 0; v < image_size; ++v)
                {
                    image_bin[v] += image[v];
                }
            }
        }

        for (int v = 0; v < image_size; ++v)
        {
            image_bin[v] /= 12;
        }

        // sava the binned memory
        std::experimental::filesystem::path curr_file(file_list.at(0));
        auto parent_path = curr_file.parent_path().string();
        auto pos = std::string(curr_file.parent_path()).rfind('/');
        const std::string bin_image_file_name = parent_path.substr(pos);
        const std::string save_path = "./out001/" + bin_image_file_name + ".raw";

        // according the bin_image_file_name calculate the coef, bhw, asca, gamma, lambda
        auto [bhw_coef, asca_coef, gamma_coef, lambda_coef] = GetCoef(bin_image_file_name);

        std::cout << "bhw: " << bhw_coef<< ", asca:" << asca_coef << ", gamma: " << gamma_coef << ", lambda: " << lambda_coef << std::endl;

        fp = fopen(save_path.c_str(), "wb");
        fwrite(image_bin, sizeof(float), image_size, fp);
        fclose(fp);
        
        image_index++;
        auto res = CalculateROI(image_bin, voxel_num_x, 300);

        // max distance ceter roi with the round 4 roi
        float max_means = std::abs(res.at(0).first - res.at(1).first);
        for (int m = 2; m < 5; ++m)
        {
            const float diff = std::abs(res.at(0).first - res.at(m).first);
            if (max_means < diff)
            {
                max_means = diff;
            }
        }

        // calculate the rois mean of mean in horizontal 
        float means_rois = 0;
        for (int m = 0; m < res.size(); ++m)
        {
            if (m == 1 || m == 3)
            {
                continue;
            }
            means_rois += res.at(m).first;
        }
        means_rois /= 9;

        float means_variance = 0;
        for (int m = 0; m < res.size(); ++m)
        {
            if (m == 1 || m == 3)
            {
                continue;
            }
            
            means_variance += std::pow((res.at(m).first - means_rois), 2);
        }
        means_variance = std::sqrt(means_variance / 9);
        

        std::cout << "max means: " << max_means << ", means variance: " << means_variance << std::endl;

        out << image_index << "\t" 
            << res.at(0).first << "\t" 
            << res.at(1).first << "\t" 
            << res.at(2).first << "\t" 
            << res.at(3).first << "\t" 
            << res.at(4).first << "\t" 
            << res.at(5).first << "\t" 
            << res.at(6).first << "\t" 
            << res.at(7).first << "\t" 
            << res.at(8).first << "\t" 
            << res.at(9).first << "\t" 
            << res.at(10).first << "\t"
            << max_means << "\t" << means_variance << "\t\t"
            << bhw_coef << "\t" << asca_coef << "\t" << gamma_coef << "\t" << lambda_coef << "\t"
            << std::endl;
    }
    out.close();

    delete[] image;
    image = nullptr;
    delete[] image_bin;
    image_bin = nullptr;
}

int main(int argc, char **argv)
{
    const std::string path = "/home/nv/gaojunbao/data/20220630/result/";
    
    CT_ImageUniform(path);

    return 0;
}