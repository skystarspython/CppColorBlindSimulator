#include <iostream>
#include <vector>
#include <cmath>
#include <map>

double linearRGB_from_sRGB(int v) {
    double fv = v / 255.0;
    if (fv < 0.04045) return fv / 12.92;
    return std::pow((fv + 0.055) / 1.055, 2.4);
}

int sRGB_from_linearRGB(double v) {
    if (v <= 0.) return 0;
    if (v >= 1.) return 255;
    if (v < 0.0031308) return 0.5 + (v * 12.92 * 255);
    return 0 + 255 * (std::pow(v, 1.0 / 2.4) * 1.055 - 0.055);
}

struct BrettelParams {
    std::vector<double> rgbCvdFromRgb_1;
    std::vector<double> rgbCvdFromRgb_2;
    std::vector<double> separationPlaneNormal;
};

std::map<std::string, BrettelParams> brettel_params = {
    {"protan", {
                    {0.14510, 1.20165, -0.34675,
                     0.10447, 0.85316, 0.04237,
                     0.00429, -0.00603, 1.00174},
                    {0.14115, 1.16782, -0.30897,
                     0.10495, 0.85730, 0.03776,
                     0.00431, -0.00586, 1.00155},
                    {0.00048, 0.00416, -0.00464}
                }},
    {"deutan", {
                    {0.36198, 0.86755, -0.22953,
                     0.26099, 0.64512, 0.09389,
                     -0.01975, 0.02686, 0.99289},
                    {0.37009, 0.88540, -0.25549,
                     0.25767, 0.63782, 0.10451,
                     -0.01950, 0.02741, 0.99209},
                    {-0.00293, -0.00645, 0.00938}
                }},
    {"tritan", {
                    {1.01354, 0.14268, -0.15622,
                     -0.01181, 0.87561, 0.13619,
                     0.07707, 0.81208, 0.11085},
                    {0.93337, 0.19999, -0.13336,
                     0.05809, 0.82565, 0.11626,
                     -0.37923, 1.13825, 0.24098},
                    {0.03960, -0.02831, -0.01129}
                }}
};

std::vector<int> brettel(const std::vector<int>& srgb, const std::string& t, double severity) {
    std::vector<double> rgb(3);
    rgb[0] = linearRGB_from_sRGB(srgb[0]);
    rgb[1] = linearRGB_from_sRGB(srgb[1]);
    rgb[2] = linearRGB_from_sRGB(srgb[2]);

    BrettelParams params = brettel_params[t];
    std::vector<double> separationPlaneNormal = params.separationPlaneNormal;
    std::vector<double> rgbCvdFromRgb_1 = params.rgbCvdFromRgb_1;
    std::vector<double> rgbCvdFromRgb_2 = params.rgbCvdFromRgb_2;

    double dotWithSepPlane = rgb[0] * separationPlaneNormal[0] + rgb[1] * separationPlaneNormal[1] + rgb[2] * separationPlaneNormal[2];
    std::vector<double> rgbCvdFromRgb = (dotWithSepPlane >= 0 ? rgbCvdFromRgb_1 : rgbCvdFromRgb_2);

    std::vector<double> rgb_cvd(3);
    rgb_cvd[0] = rgbCvdFromRgb[0] * rgb[0] + rgbCvdFromRgb[1] * rgb[1] + rgbCvdFromRgb[2] * rgb[2];
    rgb_cvd[1] = rgbCvdFromRgb[3] * rgb[0] + rgbCvdFromRgb[4] * rgb[1] + rgbCvdFromRgb[5] * rgb[2];
    rgb_cvd[2] = rgbCvdFromRgb[6] * rgb[0] + rgbCvdFromRgb[7] * rgb[1] + rgbCvdFromRgb[8] * rgb[2];

    rgb_cvd[0] = rgb_cvd[0] * severity + rgb[0] * (1.0 - severity);
    rgb_cvd[1] = rgb_cvd[1] * severity + rgb[1] * (1.0 - severity);
    rgb_cvd[2] = rgb_cvd[2] * severity + rgb[2] * (1.0 - severity);

    return { sRGB_from_linearRGB(rgb_cvd[0]), sRGB_from_linearRGB(rgb_cvd[1]), sRGB_from_linearRGB(rgb_cvd[2]) };
}

std::vector<int> monochrome_with_severity(const std::vector<int>& srgb, double severity) {
    int z = std::round(srgb[0] * 0.299 + srgb[1] * 0.587 + srgb[2] * 0.114);
    int r = z * severity + (1.0 - severity) * srgb[0];
    int g = z * severity + (1.0 - severity) * srgb[1];
    int b = z * severity + (1.0 - severity) * srgb[2];
    return { r, g, b };
}

int main() {
    std::vector<int> inputRGB = { 0, 128, 0 };
    std::string colorBlindType = "deutan";
    double severity = 1.0;

    std::vector<int> outputRGB = brettel(inputRGB, colorBlindType, severity);

    std::cout << "Input RGB: (" << inputRGB[0] << ", " << inputRGB[1] << ", " << inputRGB[2] << ")\n";
    std::cout << "Output RGB for " << colorBlindType << ": (" << outputRGB[0] << ", " << outputRGB[1] << ", " << outputRGB[2] << ")\n";

    return 0;
}
