#include <iostream>
//#include <unistd.h>
#include "src/vipss_unit.hpp"
// #include <unistd.h>
// #include <gflags/gflags.h>
#include <CLI/CLI.hpp>

typedef std::chrono::high_resolution_clock Clock;
using namespace std;

// DEFINE_string(input, "", "The input point file name.");
// DEFINE_string(output, "", "The output file name");
// DEFINE_double(lambda, 0.0, "The smooth term for vipss");
// DEFINE_int32(volumeDim, 100, "The volume dimension for surfacing marching cubes");
// DEFINE_int32(surfacing, 1, "Whether generate surface using Nature Neighbor HRBF");
// DEFINE_int32(initPV, 1, "Whether to use partial vipss to init point normals, this method runs faster");
// DEFINE_int32(hard_constraints, 1, "Whether to use normal length as hard constraints");

void SplitPath(const std::string& fullfilename,std::string &filepath);
void SplitFileName (const std::string& fullfilename,std::string &filepath,std::string &filename,std::string &extname);

int main(int argc, char** argv)
{
    VIPSSUnit vipss_unit;
    std::string out_normal_path;
    std::string out_mesh_path;
    struct 
    {
        std::string input;
        std::string output;
        double lambda = 0.0;
        int volumeDim = 100;
        bool surfacing = true;
        bool initPV = true;
        bool hardConstraints = true;
        bool rbf_base = false;
        double opt_threshold = 1e-7;
        double adgrid_threshold = 0.001;
        bool use_adgrid = true;
        bool octree_sample = false;
    }args;
    
    // gflags::ParseCommandLineFlags(&argc, &argv, true);
    CLI::App app{"local vipss Command Line"};
    app.add_option("-i, --input ", args.input, "input file")->required();
    app.add_option("-o, --output ", args.output, "output file, if not given , output will be saved in same folder as input");
    app.add_option("-l, --lambda ", args.lambda, "lambda value");
    app.add_option("-s, --volume_size ", args.volumeDim, "surface volumeDim");
    app.add_option("-P, --initPV", args.initPV, "enable or disable partial vipss for normal initialization");
    app.add_option("-S, --surfacing", args.surfacing, "reconstruct surface or not");
    app.add_option("-H, --hardConstraints", args.hardConstraints, "use hard constraints for energy optimization");
    app.add_option("-R, --use_rbfBase",args.rbf_base, "use simplified rbf base");
    app.add_option("-t, --opt_threshold",args.opt_threshold, "use simplified rbf base");
    app.add_option("-a, --adgrid_threshold",args.adgrid_threshold, "adptive gird generation threshold");
    app.add_option("-A, --use_adgrid",args.use_adgrid, "use adptive gird to generate mesh");
    app.add_option("-O, --octree_sample",args.octree_sample, "use adptive gird to generate mesh");

    
    CLI11_PARSE(app, argc, argv);
    LocalVipss::use_octree_sample_ = args.octree_sample;

    // std::cout << "LocalVipss::feature_preserve_sample_ " <<  LocalVipss::feature_preserve_sample_ << std::endl;

    vipss_unit.hard_constraints_ = args.hardConstraints;
    if(args.input.size() > 0)
    {
        std::string in_path, in_filename, in_extname;
        SplitFileName(args.input, in_path, in_filename, in_extname);
        vipss_unit.input_data_path_ = args.input;
        vipss_unit.input_data_ext_  = in_extname;
        std::cout << "input path : " << in_path << std::endl;
        std::cout << "in_filename : " << in_filename << std::endl;
        std::cout << "in_extname : " << in_extname << std::endl;
        LocalVipss::InitWithPartialVipss = args.initPV;
        LocalVipss::use_rbf_base_ = args.rbf_base;
        vipss_unit.is_surfacing_ = args.surfacing;
        vipss_unit.adgrid_threshold_ = args.adgrid_threshold;
        vipss_unit.use_adgrid_ = args.use_adgrid;
        

        
        if(args.output.size() > 0)
        {
            std::string out_path, out_filename, out_extname;
            SplitFileName(args.output, out_path, out_filename, out_extname);
            vipss_unit.out_dir_ = out_path;
            vipss_unit.file_name_ = out_filename;
            vipss_unit.local_vipss_.out_dir_ = vipss_unit.out_dir_;

            // std::cout << "output path : " << out_path << std::endl;
            // std::cout << "out_filename : " << out_filename << std::endl;
            // std::cout << "out_extname : " << out_extname << std::endl;
            vipss_unit.out_normal_path_ = out_path + "/" + out_filename + "_out_normal";
            vipss_unit.out_surface_path_ = args.output;
            vipss_unit.out_debug_path_ = out_path + "/" + out_filename + "_debug.txt";
        } else {
            vipss_unit.out_dir_ = in_path;
            vipss_unit.file_name_ = in_filename;
            vipss_unit.out_normal_path_ = in_path + "/" + in_filename + "_out_normal";
            vipss_unit.out_surface_path_ = in_path + "/" + in_filename + "_out_surface";
            vipss_unit.out_debug_path_ = in_path + "/" + in_filename + "_debug.txt";
        }
        vipss_unit.local_vipss_.out_dir_ = vipss_unit.out_dir_;
        vipss_unit.local_vipss_.filename_ = vipss_unit.file_name_;
        vipss_unit.user_lambda_ = args.lambda;
        vipss_unit.use_hrbf_surface_ = true;
        vipss_unit.volume_dim_ = args.volumeDim; 
        vipss_unit.axi_plane_ = AXI_PlANE::XYZ;
        vipss_unit.hrbf_type_ = HRBF_SURFACE_TYPE::LOCAL_HRBF_NN;
        // vu.hrbf_type_ = HRBF_SURFACE_TYPE::GLOBAL_HRBF;
        vipss_unit.Run();
    } else {
        std::cout << "There is no valid input point path !" << std::endl;
    }
    std::cout << "Success !" << std::endl;
    return 0;
}

void SplitFileName (const std::string& fullfilename,std::string &filepath,std::string &filename,std::string &extname) {
    size_t pos;
    pos = fullfilename.find_last_of('.');
    filepath = fullfilename.substr(0,pos);
    extname = fullfilename.substr(pos);
    pos = filepath.find_last_of("\\/");
    filename = filepath.substr(pos+1);
    pos = fullfilename.find_last_of("\\/");
    filepath = fullfilename.substr(0,pos+1);
    //cout<<modname<<' '<<extname<<' '<<filepath<<endl;
}

void SplitPath(const std::string& fullfilename,std::string &filepath){
    std::string filename;
    std::string extname;
    SplitFileName(fullfilename, filepath, filename,extname);
}
