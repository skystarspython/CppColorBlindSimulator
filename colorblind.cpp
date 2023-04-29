#include <iostream>
#include <vector>
#include <cmath>
#include <map>

struct RGB {
    double r;
    double g;
    double b;
};

struct BlindType {
    double cpu;
    double cpv;
    double am;
    double ayi;
};

std::map<std::string, BlindType> rBlind = {
    {"protan", {0.735, 0.265, 1.273463, -0.073894}},
    {"deutan", {1.14, -0.14, 0.968437, 0.003331}},
    {"tritan", {0.171, -0.003, 0.062921, 0.292119}}
};

std::vector<double> powGammaLookup(256);

void initializePowGammaLookup() {
    for (int i = 0; i < 256; ++i) {
        powGammaLookup[i] = std::pow(i / 255.0, 2.2);
    }
}

double inversePow(double v) {
    return (255 * (v <= 0 ? 0 : v >= 1 ? 1 : std::pow(v, 1 / 2.2)));
}

RGB blindMK(RGB rgb, const std::string& t) {
    double gamma = 2.2;
    double wx = 0.312713;
    double wy = 0.329016;
    double wz = 0.358271;

    double cr = powGammaLookup[static_cast<int>(rgb.r)];
    double cg = powGammaLookup[static_cast<int>(rgb.g)];
    double cb = powGammaLookup[static_cast<int>(rgb.b)];

    double cx = (0.430574 * cr + 0.341550 * cg + 0.178325 * cb);
    double cy = (0.222015 * cr + 0.706655 * cg + 0.071330 * cb);
    double cz = (0.020183 * cr + 0.129553 * cg + 0.939180 * cb);

    double sum_xyz = cx + cy + cz;
    double cu = 0;
    double cv = 0;

    if (sum_xyz != 0) {
        cu = cx / sum_xyz;
        cv = cy / sum_xyz;
    }

    double nx = wx * cy / wy;
    double nz = wz * cy / wy;
    double clm;
    double dy = 0;

    if (cu < rBlind[t].cpu) {
        clm = (rBlind[t].cpv - cv) / (rBlind[t].cpu - cu);
    }
    else {
        clm = (cv - rBlind[t].cpv) / (cu - rBlind[t].cpu);
    }

    double clyi = cv - cu * clm;
    double du = (rBlind[t].ayi - clyi) / (clm - rBlind[t].am);
    double dv = (clm * du) + clyi;

    double sx = du * cy / dv;
    double sy = cy;
    double sz = (1 - (du + dv)) * cy / dv;

    double sr = (3.063218 * sx - 1.393325 * sy - 0.475802 * sz);
    double sg = (-0.969243 * sx + 1.875966 * sy + 0.041555 * sz);
    double sb = (0.067871 * sx - 0.228834 * sy + 1.069251 * sz);

    double dx = nx - sx;
    double dz = nz - sz;

    double dr = (3.063218 * dx - 1.393325 * dy - 0.475802 * dz);
    double dg = (-0.969243 * dx + 1.875966 * dy + 0.041555 * dz);
    double db = (0.067871 * dx - 0.228834 * dy + 1.069251 * dz);

    double adjr = dr ? ((sr < 0 ? 0 : 1) - sr) / dr : 0;
    double adjg = dg ? ((sg < 0 ? 0 : 1) - sg) / dg : 0;
    double adjb = db ? ((sb < 0 ? 0 : 1) - sb) / db : 0;

    double adjust = std::max(
        ((adjr > 1 || adjr < 0) ? 0 : adjr),
        std::max(((adjg > 1 || adjg < 0) ? 0 : adjg),
        ((adjb > 1 || adjb < 0) ? 0 : adjb))
    );

    sr = sr + (adjust * dr);
    sg = sg + (adjust * dg);
    sb = sb + (adjust * db);

    return { inversePow(sr), inversePow(sg), inversePow(sb) };
}

RGB anomylize(RGB a, RGB b) {
    double v = 1.75, d = v * 1 + 1;

    return {
        (v * b.r + a.r * 1) / d,
        (v * b.g + a.g * 1) / d,
        (v * b.b + a.b * 1) / d
    };
}

RGB monochrome(RGB r) {
    double z = std::round(r.r * 0.299 + r.g * 0.587 + r.b * 0.114);
    return { z, z, z };
}

int main() {
    initializePowGammaLookup();

    RGB input_rgb = { 0, 128, 0 }; // 输入RGB颜色
    std::string color_blindness_type = "protan"; // 色盲类型（如：protan，deutan，tritan）

    RGB output_rgb = blindMK(input_rgb, color_blindness_type);

    std::cout << "Output RGB: (" << output_rgb.r << ", " << output_rgb.g << ", " << output_rgb.b << ")" << std::endl;

    return 0;
}
