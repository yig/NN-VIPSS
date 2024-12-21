#pragma once

#include <iostream>
#include <array>
#include <cmath>

struct RGBColor {
    double r, g, b;
    
    RGBColor(double red, double green, double blue) : r(red), g(green), b(blue) {}

    // Linear interpolation between two colors
    static RGBColor interpolate(const RGBColor& c1, const RGBColor& c2, double t) {
        return RGBColor(
            c1.r + t * (c2.r - c1.r),
            c1.g + t * (c2.g - c1.g),
            c1.b + t * (c2.b - c1.b)
        );
    }
};

// Blend function
RGBColor ErrorColorBlend(double u) {
    if (u < 0.0 || u > 1.0) {
        throw std::out_of_range("u must be in the range [0, 1]");
    }
    // Define the keypoints and corresponding colors
    std::array<std::pair<double, RGBColor>, 7> points = {{
        {0.0,   RGBColor(0.0, 0.0, 9.0 / 16.0)}, // RGBColor[0, 0, 9/16]
        {1.0/9.0, RGBColor(0.0, 0.0, 1.0)},      // Blue
        {23.0/63.0, RGBColor(0.0, 1.0, 1.0)},    // Cyan
        {13.0/21.0, RGBColor(1.0, 1.0, 0.0)},    // Yellow
        {47.0/63.0, RGBColor(1.0, 0.647, 0.0)},  // Orange
        {55.0/63.0, RGBColor(1.0, 0.0, 0.0)},    // Red
        {1.0,   RGBColor(0.5, 0.0, 0.0)}         // RGBColor[1/2, 0, 0]
    }};

    // Find the interval for u
    for (size_t i = 0; i < points.size() - 1; ++i) {
        const auto& [t1, c1] = points[i];
        const auto& [t2, c2] = points[i + 1];
        if (u >= t1 && u <= t2) {
            double t = (u - t1) / (t2 - t1); // Normalize u in the interval [t1, t2]
            return RGBColor::interpolate(c1, c2, t);
        }
    }

    // Default (should never be reached if u is within [0, 1])
    return points.back().second;
}

// int main() {
//     try {
//         // Test the Blend function
//         for (double u = 0.0; u <= 1.0; u += 0.1) {
//             RGBColor color = Blend(u);
//             std::cout << "u: " << u
//                       << " -> RGB(" << color.r << ", " << color.g << ", " << color.b << ")\n";
//         }
//     } catch (const std::exception& e) {
//         std::cerr << e.what() << std::endl;
//     }

//     return 0;
// }
