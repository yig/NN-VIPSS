#include "adgrid.h"
#include <iostream>

// struct
//     {
//         std::string grid_file;
//         std::string function_file;
//         double threshold;
//         double alpha = std::numeric_limits<double>::infinity();
//         int max_elements = -1;
//         double smallest_edge_length = 0;
//         std::string method = "IA";
//         std::string csg_file;
//         bool bfs = false;
//         bool dfs = false;
//         bool curve_network = false;
//         bool discretize_later = false;
//         //bool analysis_mode = false;
//     } args;

void GenerateAdaptiveGridOut(const std::array<size_t, 3>& resolution, 
                             const std::array<double, 3>& bbox_min, 
                             const std::array<double, 3>& bbox_max,
                             const std::string& outdir,
                             const std::string& filename)
{
    std::cout << "start to call  GenerateAdaptiveGridOut" << std::endl;
    std::cout << "bbox min " << bbox_min[0] << " " << bbox_min[1] << " " << bbox_min[2] << std::endl;
    std::cout << "bbox max " << bbox_max[0] << " " << bbox_max[1] << " " << bbox_max[2] << std::endl;

    double expand_scale = 0.2;
    double dx = bbox_max[0] - bbox_min[0];
    double dy = bbox_max[1] - bbox_min[1];
    double dz = bbox_max[2] - bbox_min[2];

    std::array<double, 3> expand_bbox_min = {bbox_min[0] - expand_scale * dx, 
                                            bbox_min[1] - expand_scale * dy,
                                            bbox_min[2] - expand_scale * dz};
    std::array<double, 3> expand_bbox_max = {bbox_max[0] + expand_scale * dx, 
                                            bbox_max[1] + expand_scale * dy,
                                            bbox_max[2] + expand_scale * dz};


    std::array<double, 3> new_bbox_min;

    mtet::MTetMesh grid = generate_tet_mesh(resolution, expand_bbox_min, expand_bbox_max, grid_mesh::TET5);

    int max_elements = std::numeric_limits<int>::max();
    double threshold = 0.0005;
    double alpha = 1.0;
    
    llvm_vecsmall::SmallVector<csg_unit, 20> csg_tree = {};
    /// Read implicit function
    std::vector<std::shared_ptr<HRBFDistanceFunction>> functions;
    std::shared_ptr<HRBFDistanceFunction> hrbf_func = std::make_shared<HRBFDistanceFunction>();
    functions.push_back(hrbf_func);

    const size_t funcNum = functions.size();
    ///
    /// the lambda function for function evaluations
    ///  @param[in] data            The 3D coordinate
    ///  @param[in] funcNum         The number of functions
    ///
    ///  @return        A vector of `Eigen::RowVector4d`.The vector size is the function number. Each eigen vector represents the value at 0th index and gradients at {1, 2, 3} index.
    auto implicit_func = [&](std::span<const Scalar, 3> data, size_t funcNum){
        llvm_vecsmall::SmallVector<Eigen::RowVector4d, 20> vertex_eval(funcNum);
        for(size_t funcIter = 0; funcIter < funcNum; funcIter++){
            auto &func = functions[funcIter];
            Eigen::Vector4d eval;
            eval[0] = func->evaluate_gradient(data[0], data[1], data[2], eval[1], eval[2], eval[3]);
            vertex_eval[funcIter] = eval;
        }
        return vertex_eval;
    };
    
    ///
    /// the lambda function for csg tree iteration/evaluation.
    /// @param[in] funcInt          Given an input of value range std::array<double, 2> for an arbitrary number of functions
    /// @return   A value range of this CSG operation in a form of `std::array<double, 2>` and a list of active function in a form of    `llvm_vecsmall::SmallVector<int, 20>>`
    ///
    auto csg_func = [&](llvm_vecsmall::SmallVector<std::array<double, 2>, 20> funcInt){
        if (1){
            throw std::runtime_error("ERROR: no csg file provided");
            std::pair<std::array<double, 2>, llvm_vecsmall::SmallVector<int, 20>> null_csg = {{},{}};
            return null_csg;
        }else{
            return iterTree(csg_tree, 1, funcInt);
        }
    };

    int mode = 0;
    //perform main grid refinement algorithm:
    tet_metric metric_list;
    //an array of 10 timings: {total time getting the multiple indices, total time,time spent on single function, time spent on double functions, time spent on triple functions time spent on double functions' zero crossing test, time spent on three functions' zero crossing test, total subdivision time, total evaluation time,total splitting time}
    std::array<double, timer_amount> profileTimer = {0,0,0,0,0,0,0,0,0,0};
    std::cout << "start to call gridRefine .... "<< std::endl;
    if (!gridRefine(mode, false, threshold, alpha, max_elements, funcNum, implicit_func, csg_func, grid, metric_list, profileTimer))
    {
        throw std::runtime_error("ERROR: unsuccessful grid refinement");
    }
    std::cout << " finished "<< std::endl;
    // save timing records
    save_timings("timings.json",time_label, profileTimer);
    //profiled time(see details in time.h) and profiled number of calls to zero
    for (int i = 0; i < profileTimer.size(); i++){
        timeProfileName time_type = static_cast<timeProfileName>(i);
        std::cout << time_label[i] << ": " << profileTimer[i] << std::endl;
    }
    // save tet metrics
    save_metrics("stats.json", tet_metric_labels, metric_list);
    
    if (1){
        std::string outfile = outdir + "/" + filename;
        /// save the grid output for discretization tool
        save_mesh_json(outfile + "grid.json", grid);
        std::cout << "saved grid json file .... "<< std::endl;
        /// save the grid output for isosurfacing tool
        save_function_json(outfile + "function_value.json", grid, metric_list.vertex_func_grad_map, funcNum);
        /// write grid and active tets
        mtet::save_mesh(outfile + "tet_grid.msh",  grid);
        mtet::save_mesh(outfile + "active_tets.msh", grid, std::span<mtet::TetId>(metric_list.activeTetId));
    }
}
