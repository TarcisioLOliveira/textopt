/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of textopt.
 *
 *   textopt is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   textopt is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with textopt.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <algorithm>
#include <iostream>
#include "newton.hpp"
#include <iostream>
#include "texture_shallow.hpp"
#include "param.hpp"
#include "smooth.hpp"

namespace texture_shallow{

void map_exact(std::vector<double>& map_z, const std::vector<double>& orig_z, double ap, double vc){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    // Calculate maximum width
    const double z1 = -std::sqrt(r*r - param::y1*param::y1) + r;
    const double z2 = -std::sqrt(r*r - param::y2*param::y2) + r;
    const double z12min = std::min(z1, z2);
    const double z12max = std::max(z1, z2);
    const double h_max = ap + Az_uet;
    if(h_max < z12min){
        // Depth fits entirely within tool radius
        w_max = 2*std::sqrt(r*r - (h_max - r)*(h_max - r));
    } else if(h_max < z12max){
        // Depth fits partially within tool radius
        // Can only happen for alpha1 != alpha2
        w_max = std::sqrt(r*r - (h_max - r)*(h_max - r));
        if(alpha1 <= alpha2){
            w_max += (h_max - b1off)/tan1;
        } else {
            w_max += (h_max - b2off)/tan2;
        }
    } else {
        // Depth encompasses tool's straight edges
        w_max = (h_max - b1off)/tan1 + (h_max - b2off)/tan2;
    }

    f = w_max;

    const double delta_uet = vc/f_uet;

    max_z = 0;
    min_z = -h_max;

    const double line_root1_const = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double line_root2_const = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    const double perimeter = 2*M_PI*cylinder_radius;

    #pragma omp parallel for
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        const double mult = smooth::floor(y/w_max + 0.5);

        // Distance travelled by tool
        const double xoffset_uet = mult*perimeter;

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            // Consider distance travelled by tool
            const double xcirc = x + xoffset_uet;

            // Tool path
            const double t0 = (std::floor(xcirc / delta_uet) + 0.5)/(f_uet);
            const double t = newton_t(xcirc, t0, vc);
            double z_uet = zt(t);

            // Clearance angle
            const double x_uet_c = xcirc + 0.5*delta_uet - std::floor(xcirc / delta_uet + 0.5)*delta_uet;
            const double z_clear = tanc*(delta_uet - x_uet_c);

            // Final UET height
            const double z_uet_min = std::min(z_uet, z_clear);

            const double line_root1 = line_root1_const + z_uet_min;
            const double line_root2 = line_root2_const + z_uet_min;

            // Tool shape
            double shape_z;
            if(y <= y1 + mult*f){
                shape_z = -tan1*(y - mult*f) + line_root1;
            } else if(y <= y2 + mult*f){
                const double yy = y - mult*f;
                shape_z = -std::sqrt(r*r - yy*yy) + r - ap + z_uet_min;
            } else {
                shape_z = tan2*(y - mult*f) + line_root2;
            }

            // Get pixels
            double& z = map_z[X + Y*tex_width];
            const double& oz = orig_z[X + Y*tex_width];
            // Write results
            z = std::min(oz, shape_z);
        }
    }

    max_z = *std::max_element(map_z.begin(), map_z.end());
}

void map(std::vector<double>& map_z, const std::vector<double>& orig_z, double ap, double vc){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    // Calculate maximum width
    const double z1 = -std::sqrt(r*r - param::y1*param::y1) + r;
    const double z2 = -std::sqrt(r*r - param::y2*param::y2) + r;
    const double z12min = std::min(z1, z2);
    const double z12max = std::max(z1, z2);
    const double h_max = ap + Az_uet;
    if(h_max < z12min){
        // Depth fits entirely within tool radius
        w_max = 2*std::sqrt(r*r - (h_max - r)*(h_max - r));
    } else if(h_max < z12max){
        // Depth fits partially within tool radius
        // Can only happen for alpha1 != alpha2
        w_max = std::sqrt(r*r - (h_max - r)*(h_max - r));
        if(alpha1 <= alpha2){
            w_max += (h_max - b1off)/tan1;
        } else {
            w_max += (h_max - b2off)/tan2;
        }
    } else {
        // Depth encompasses tool's straight edges
        w_max = (h_max - b1off)/tan1 + (h_max - b2off)/tan2;
    }

    f = w_max;

    const double delta_uet = vc/f_uet;

    max_z = 0;
    min_z = -h_max;

    const double line_root1_const = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double line_root2_const = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    const double perimeter = 2*M_PI*cylinder_radius;

    #pragma omp parallel for
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        const double mult = smooth::floor(y/w_max + 0.5);

        // Distance travelled by tool
        const double xoffset_uet = mult*perimeter;

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            // Consider distance travelled by tool
            const double xcirc = x + xoffset_uet;

            // Tool path
            const double t0 = (std::floor(xcirc / delta_uet) + 0.5)/(f_uet); // Initial estimate of t for Newton
            const double t = newton_t(xcirc, t0, vc);
            double z_uet = zt(t);

            // Clearance angle
            const double x_uet_c = xcirc + 0.5*delta_uet - smooth::floor(xcirc / delta_uet + 0.5)*delta_uet;
            const double z_clear = tanc*(delta_uet - x_uet_c);

            // Final UET height
            const double z_uet_min = smooth::min({z_uet, z_clear});

            const double line_root1 = line_root1_const + z_uet_min;
            const double line_root2 = line_root2_const + z_uet_min;

            // Tool shape
            double shape_z;
            if(y <= y1 + mult*f){
                shape_z = -tan1*(y - mult*f) + line_root1;
            } else if(y <= y2 + mult*f){
                const double yy = y - mult*f;
                shape_z = -std::sqrt(r*r - yy*yy) + r - ap + z_uet_min;
            } else {
                shape_z = tan2*(y - mult*f) + line_root2;
            }

            // Get pixels
            double& z = map_z[X + Y*tex_width];
            const double& oz = orig_z[X + Y*tex_width];
            // Write results
            z = smooth::min({oz, shape_z});
        }
    }

    max_z = *std::max_element(map_z.begin(), map_z.end());
}

void dzdvc(const std::vector<double>& orig_z, double ap, double vc, std::vector<double>& dzdvc){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    // Calculate maximum width
    const double z1 = -std::sqrt(r*r - param::y1*param::y1) + r;
    const double z2 = -std::sqrt(r*r - param::y2*param::y2) + r;
    const double z12min = std::min(z1, z2);
    const double z12max = std::max(z1, z2);
    const double h_max = ap + Az_uet;
    if(h_max < z12min){
        // Depth fits entirely within tool radius
        w_max = 2*std::sqrt(r*r - (h_max - r)*(h_max - r));
    } else if(h_max < z12max){
        // Depth fits partially within tool radius
        // Can only happen for alpha1 != alpha2
        w_max = std::sqrt(r*r - (h_max - r)*(h_max - r));
        if(alpha1 <= alpha2){
            w_max += (h_max - b1off)/tan1;
        } else {
            w_max += (h_max - b2off)/tan2;
        }
    } else {
        // Depth encompasses tool's straight edges
        w_max = (h_max - b1off)/tan1 + (h_max - b2off)/tan2;
    }

    f = w_max;

    const double delta_uet = vc/f_uet;
    const double ddelta_uet = 1.0/f_uet;

    const double line_root1_const = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double line_root2_const = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    const double perimeter = 2*M_PI*cylinder_radius;

    #pragma omp parallel for
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        const double mult = smooth::floor(y/w_max + 0.5);

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            // Distance travelled by tool
            const double xoffset_uet = mult*perimeter;

            // Consider distance travelled by tool
            const double xcirc = x + xoffset_uet;

            // Tool path
            const double t0 = (std::floor(xcirc / delta_uet) + 0.5)/(f_uet); // Initial estimate of t for Newton
            const double t = newton_t(xcirc, t0, vc);
            double z_uet = zt(t);
            double dz_uet = dzdt(t)*dtdvc(t, vc);

            // Clearance angle
            const double x_uet_c = xcirc + 0.5*delta_uet - smooth::floor(xcirc / delta_uet + 0.5)*delta_uet;
            const double dx_uet_c = 0.5*ddelta_uet - (smooth::floor(xcirc / delta_uet + 0.5)*ddelta_uet + smooth::floor_deriv(xcirc / delta_uet + 0.5)*delta_uet*(-xcirc/(delta_uet*delta_uet)));

            const double z_clear = tanc*(delta_uet - x_uet_c);
            const double dz_clear = tanc*(ddelta_uet - dx_uet_c);

            // Final UET height
            const double z_uet_min = smooth::min({z_uet, z_clear});
            const double dz_uet_min = smooth::min_deriv({z_uet, z_clear},{dz_uet, dz_clear});

            const double line_root1 = line_root1_const + z_uet_min;
            const double line_root2 = line_root2_const + z_uet_min;
            const double dline_root1 = dz_uet_min;
            const double dline_root2 = dz_uet_min;

            // Tool shape
            double shape_z;
            double dshape_z;
            if(y <= y1 + mult*f){
                shape_z = -tan1*(y - mult*f) + line_root1;
                dshape_z = dline_root1;
            } else if(y <= y2 + mult*f){
                const double yy = y - mult*f;
                shape_z = -std::sqrt(r*r - yy*yy) + r - ap + z_uet_min;
                dshape_z = dz_uet_min;
            } else {
                shape_z = tan2*(y - mult*f) + line_root2;
                dshape_z = dline_root2;
            }

            // Get pixels
            double& dz = dzdvc[X + Y*tex_width];
            const double& oz = orig_z[X + Y*tex_width];
            // Write results
            dz = smooth::min_deriv({oz, shape_z}, 1)*dshape_z;
        }
    }
}

void dzdap(const std::vector<double>& orig_z, double ap, double vc, std::vector<double>& dzdap){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    // Calculate maximum width
    const double z1 = -std::sqrt(r*r - param::y1*param::y1) + r;
    const double z2 = -std::sqrt(r*r - param::y2*param::y2) + r;
    const double z12min = std::min(z1, z2);
    const double z12max = std::max(z1, z2);
    const double h_max = ap + Az_uet;

    double dw_max;
    if(h_max < z12min){
        // Depth fits entirely within tool radius
        w_max = 2*std::sqrt(r*r - (h_max - r)*(h_max - r));
        dw_max = 2*(-(h_max - r))/std::sqrt(r*r - (h_max - r)*(h_max - r));
    } else if(h_max < z12max){
        // Depth fits partially within tool radius
        // Can only happen for alpha1 != alpha2
        w_max = std::sqrt(r*r - (h_max - r)*(h_max - r));
        dw_max = 2*(-(h_max - r))/std::sqrt(r*r - (h_max - r)*(h_max - r));
        if(alpha1 <= alpha2){
            w_max += (h_max - b1off)/tan1;
            dw_max += 1.0/tan1;
        } else {
            w_max += (h_max - b2off)/tan2;
            dw_max += 1.0/tan2;
        }
    } else {
        // Depth encompasses tool's straight edges
        w_max = (h_max - b1off)/tan1 + (h_max - b2off)/tan2;
        dw_max = 1.0/tan1 + 1.0/tan2;
    }

    f = w_max;
    const double df = dw_max;

    const double delta_uet = vc/f_uet;

    min_z = -ap;

    const double line_root1_const = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double line_root2_const = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;
    const double dline_root1_const = -1;
    const double dline_root2_const = -1;

    const double perimeter = 2*M_PI*cylinder_radius;

    #pragma omp parallel for
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        const double mult = smooth::floor(y/w_max + 0.5);
        const double dmult = smooth::floor_deriv(y/w_max + 0.5)*(-y*dw_max/(w_max*w_max));

        // Distance travelled by tool
        const double xoffset_uet = mult*perimeter;
        const double dxoffset_uet = dmult*perimeter;

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            // Consider distance travelled by tool
            const double xcirc = x + xoffset_uet;
            const double dxcirc = dxoffset_uet;

            // Tool path
            const double t0 = (std::floor(xcirc / delta_uet) + 0.5)/(f_uet); // Initial estimate of t for Newton
            const double t = newton_t(xcirc, t0, vc);
            const double dtdf = dxcirc/dxdt(t, vc);
            double z_uet = zt(t);
            double dz_uet = dzdt(t)*dtdf;

            // Clearance angle
            const double x_uet_c = xcirc + 0.5*delta_uet - smooth::floor(xcirc / delta_uet + 0.5)*delta_uet;
            const double dx_uet_c = dxcirc  - smooth::floor_deriv(xcirc / delta_uet + 0.5)*dxcirc;

            const double z_clear = tanc*(delta_uet - x_uet_c);
            const double dz_clear = tanc*(-dx_uet_c);

            // Final UET height
            const double z_uet_min = smooth::min({z_uet, z_clear});
            const double dz_uet_min = smooth::min_deriv({z_uet, z_clear},{dz_uet, dz_clear});

            const double line_root1 = line_root1_const + z_uet_min;
            const double line_root2 = line_root2_const + z_uet_min;
            const double dline_root1 = dline_root1_const + dz_uet_min;
            const double dline_root2 = dline_root2_const + dz_uet_min;

            // Tool shape
            double shape_z;
            double dshape_z;
            if(y <= y1 + mult*f){
                shape_z = -tan1*(y - mult*f) + line_root1;
                dshape_z = tan1*(dmult*f + mult*df) + dline_root1;
            } else if(y <= y2 + mult*f){
                const double yy = y - mult*f;
                const double dyy = -(dmult*f + mult*df);
                shape_z = -std::sqrt(r*r - yy*yy) + r - ap + z_uet_min;
                dshape_z = yy*dyy/std::sqrt(r*r - yy*yy) - 1.0 + dz_uet_min;
            } else {
                shape_z = tan2*(y - mult*f) + line_root2;
                dshape_z = -tan2*(dmult*f + mult*df) + dline_root2;
            }

            // Get pixels
            double& dz = dzdap[X + Y*tex_width];
            const double& oz = orig_z[X + Y*tex_width];
            // Write results
            dz = smooth::min_deriv({oz, shape_z}, 1)*dshape_z;
        }
    }
}

}
