#include <span>
#include <queue>
#include <optional>
#include <CLI/CLI.hpp>

#include "adgrid/timer.h"
#include "adgrid/csg.h"
#include "adgrid/grid_mesh.h"
#include "adgrid/grid_refine.h"
#include "adgrid/external/implicit_functions/implicit_functions.h"
#include "local_vipss.hpp"
#include "adgrid/external/implicit_functions/Hermite_RBF.h"


class HRBFDistanceFunction : public ImplicitFunction<double>
{
public:
    HRBFDistanceFunction()
    {
        
    }

    double evaluate(double x, double y, double z) const override
    {
        // std::cout << "eval pt : " << x << " " << y << " " << z << std::endl;
        R3Pt newPt(x, y, z);
        return LocalVipss::NNDistFunction(newPt);
    }

    double evaluate_gradient(double x, double y, double z, double &gx, double &gy, double &gz) const override
    {
        double delt = 1e-4;
        // std::cout << "pt : " << x << " " << y << " " << z << std::endl;
        R3Pt newPt(x, y, z);
        double dist_val = LocalVipss::NNDistFunction(newPt);
        // std::cout << "dist_val : " << dist_val << std::endl;
        R3Pt newPt_x(x + delt, y, z);
        double dist_x = LocalVipss::NNDistFunction(newPt_x);
        R3Pt newPt_y(x, y + delt, z);
        double dist_y = LocalVipss::NNDistFunction(newPt_y);
        R3Pt newPt_z(x, y, z + delt);
        double dist_z = LocalVipss::NNDistFunction(newPt_z);

        gx = (dist_x - dist_val) / delt;
        gy = (dist_y - dist_val) / delt;
        gz = (dist_z - dist_val) / delt;

        // std::cout << "grad : " << gx << " " << gy << " " << gz << std::endl;
        return dist_val;
    }

private:
    double radius_;
};


void GenerateAdaptiveGridOut(const std::array<size_t, 3>& resolution, 
                            const std::array<double, 3>& bbox_min,
                            const std::array<double, 3>& bbox_max,
                            const std::string& outdir,
                            const std::string& fillname);

