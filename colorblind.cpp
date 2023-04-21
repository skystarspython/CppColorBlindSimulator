#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

vector<double> powGammaLookup(256);

struct ColorMatrix {
    vector<double> R;
    vector<double> G;
    vector<double> B;
};

struct RBlind {
    double cpu, cpv, am, ayi;
};

map<string, ColorMatrix> colorMatrixMatrixes = {
    {"Normal", {{100, 0, 0}, {0, 100, 0}, {0, 0, 100}}},
    {"Protanopia", {{56.667, 43.333, 0}, {55.833, 44.167, 0}, {0, 24.167, 75.833}}},
    {"Protanomaly", {{81.667, 18.333, 0}, {33.333, 66.667, 0}, {0, 12.5, 87.5}}},
    {"Deuteranopia", {{62.5, 37.5, 0}, {70, 30, 0}, {0, 30, 70}}},
    {"Deuteranomaly", {{80, 20, 0}, {25.833, 74.167, 0}, {0, 14.167, 85.833}}},
    {"Tritanopia", {{95, 5, 0}, {0, 43.333, 56.667}, {0, 47.5, 52.5}}},
    {"Tritanomaly", {{96.667, 3.333, 0}, {0, 73.333, 26.667}, {0, 18.333, 81.667}}},
    {"Achromatopsia", {{29.9, 58.7, 11.4}, {29.9, 58.7, 11.4}, {29.9, 58.7, 11.4}}},
    {"Achromatomaly", {{61.8, 32, 6.2}, {16.3, 77.5, 6.2}, {16.3, 32, 51.6}}},
};

map<string, RBlind> rBlind = {
    {"protan", {0.735, 0.265, 1.273463, -0.073894}},
    {"deutan", {1.14, -0.14, 0.968437, 0.003331}},
    {"tritan", {0.171, -0.003, 0.062921, 0.292119}},
};

void initPowGammaLookup() {
    for (int i = 0; i < 256; i++) {
        powGammaLookup[i] = pow(i / 255.0, 2.2);
    }
}

vector<double> matrixFunction(const ColorMatrix& a, const vector<double>& b) {
    vector<double> result(3);
    result[0] = (b[0] * a.R[0] + b[1] * a.R[1] + b[2] * a.R[2]) / 100;
    result[1] = (b[0] * a.G[0] + b[1] * a.G[1] + b[2] * a.G[2]) / 100;
    result[2] = (b[0] * a.B[0] + b[1] * a.B[1] + b[2] * a.B[2]) / 100;
    return result;
}

double inversePow(double a) {
    return 255 * (a <= 0 ? 0 : a >= 1 ? 1 : pow(a, 1 / 2.2));
}

vector<double> blindMK(const vector<double>& a, const string& str) {
    const RBlind& rb = rBlind.at(str);
    double r = a[0], g = a[1], b = a[2];
    double j = powGammaLookup[r], k = powGammaLookup[g], l = powGammaLookup[b];
    double m = 0.430574 * j + 0.34155 * k + 0.178325 * l;
    double n = 0.222015 * j + 0.706655 * k + 0.07133 * l;
    double o = 0.020183 * j + 0.129553 * k + 0.93918 * l;
    double p = m + n + o;
    double q = 0, s = 0;
    if (p != 0) {
        q = m / p;
        s = n / p;
    }
    double u, t = 0;
    u = q < rb.cpu ? (rb.cpv - s) / (rb.cpu - q) : (s - rb.cpv) / (q - rb.cpu);
    double v = s - q * u;
    double x = (rb.ayi - v) / (u - rb.am);
    double y = u * x + v;
    double z = x * n / y;
    double A = n;
    double B = (1 - (x + y)) * n / y;
    double C = 3.063218 * z - 1.393325 * A - 0.475802 * B;
    double D = -0.969243 * z + 1.875966 * A + 0.041555 * B;
    double E = 0.067871 * z - 0.228834 * A + 1.069251 * B;
    double F = 0, G = 0, H = 0;
    F = s - z;
    G = t - B;
    double dr = 3.063218 * F - 1.393325 * t - 0.475802 * G;
    double dg = -0.969243 * F + 1.875966 * t + 0.041555 * G;
    double db = 0.067871 * F - 0.228834 * t + 1.069251 * G;
    double I = dr ? ((C < 0 ? 0 : 1) - C) / dr : 0;
    double J = dg ? ((D < 0 ? 0 : 1) - D) / dg : 0;
    double K = db ? ((E < 0 ? 0 : 1) - E) / db : 0;
    double L = max(max(I > 1 || I < 0 ? 0 : I, J > 1 || J < 0 ? 0 : J), K > 1 || K < 0 ? 0 : K);
    C += L * dr;
    D += L * dg;
    E += L * db;
    return { inversePow(C), inversePow(D), inversePow(E) };
}

vector<double> anomylize(const vector<double>& a, const vector<double>& b) {
    double c = 1.75;
    double d = 1 * c + 1;
    vector<double> result(3);
    result[0] = (c * b[0] + 1 * a[0]) / d;
    result[1] = (c * b[1] + 1 * a[1]) / d;
    result[2] = (c * b[2] + 1 * a[2]) / d;
    return result;
}

vector<double> monochrome(const vector<double>& a) {
    double b = round(0.299 * a[0] + 0.587 * a[1] + 0.114 * a[2]);
    return { b, b, b };
}

vector<double> applyColorBlindFunction(const vector<double>& a, const string& colorBlindType) {
    if (colorBlindType == "Normal") {
        return a;
    }
    else if (colorBlindType == "Protanopia") {
        return blindMK(a, "protan");
    }
    else if (colorBlindType == "Protanomaly") {
        return anomylize(a, blindMK(a, "protan"));
    }
    else if (colorBlindType == "Deuteranopia") {
        return blindMK(a, "deutan");
    }
    else if (colorBlindType == "Deuteranomaly") {
        return anomylize(a, blindMK(a, "deutan"));
    }
    else if (colorBlindType == "Tritanopia") {
        return blindMK(a, "tritan");
    }
    else if (colorBlindType == "Tritanomaly") {
        return anomylize(a, blindMK(a, "tritan"));
    }
    else if (colorBlindType == "Achromatopsia") {
        return monochrome(a);
    }
    else if (colorBlindType == "Achromatomaly") {
        return anomylize(a, monochrome(a));
    }
    else {
        throw invalid_argument("Unknown color blind type.");
    }
}

vector<int> colorBlindConvert(const vector<int>& rgb, const string& colorBlindType) {
    vector<double> input(rgb.begin(), rgb.end());
    vector<double> output = applyColorBlindFunction(input, colorBlindType);
    vector<int> result(output.size());
    transform(output.begin(), output.end(), result.begin(), [](double d) { return static_cast<int>(round(d)); });
    return result;
}

int main() {
    initPowGammaLookup();
    vector<int> inputRGB = { 255, 0, 0 };
    string colorBlindType = "Protanopia";
    vector<int> outputRGB = colorBlindConvert(inputRGB, colorBlindType);
    cout << "Original RGB: (" << inputRGB[0] << ", " << inputRGB[1] << ", " << inputRGB[2] << ")\n";
    cout << "Color blind type: " << colorBlindType << "\n";
    cout << "Output RGB: (" << outputRGB[0] << ", " << outputRGB[1] << ", " << outputRGB[2] << ")\n";
    return 0;
}
