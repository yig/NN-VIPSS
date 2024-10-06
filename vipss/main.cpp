#include <iostream>
//#include <unistd.h>
#include "src/vipss_unit.hpp"
// #include <unistd.h>
#include <gflags/gflags.h>


typedef std::chrono::high_resolution_clock Clock;
using namespace std;

DEFINE_string(input, "", "The input point file name.");
DEFINE_string(output, "", "The output file name");
DEFINE_double(lambda, 0.0, "The smooth term for vipss");
DEFINE_int32(volumeDim, 100, "The volume dimension for surfacing marching cubes");
DEFINE_int32(surfacing, 1, "Whether generate surface using Nature Neighbor HRBF");
DEFINE_int32(initPV, 1, "Whether to use partial vipss to init point normals, this method runs faster");


void SplitPath(const std::string& fullfilename,std::string &filepath);
void SplitFileName (const std::string& fullfilename,std::string &filepath,std::string &filename,std::string &extname);

int main(int argc, char** argv)
{
    VIPSSUnit vipss_unit;
    std::string out_normal_path;
    std::string out_mesh_path;
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    if(FLAGS_input.size() > 0)
    {
        std::string in_path, in_filename, in_extname;
        SplitFileName(FLAGS_input, in_path, in_filename, in_extname);
        vipss_unit.input_data_path_ = FLAGS_input;
        vipss_unit.input_data_ext_  = in_extname;
        std::cout << "input path : " << in_path << std::endl;
        std::cout << "in_filename : " << in_filename << std::endl;
        std::cout << "in_extname : " << in_extname << std::endl;

        LocalVipss::InitWithPartialVipss = FLAGS_initPV;

        if(FLAGS_output.size() > 0)
        {
            vipss_unit.is_surfacing_ = FLAGS_surfacing;
            std::string out_path, out_filename, out_extname;
            SplitFileName(FLAGS_output, out_path, out_filename, out_extname);
            // std::cout << "output path : " << out_path << std::endl;
            // std::cout << "out_filename : " << out_filename << std::endl;
            // std::cout << "out_extname : " << out_extname << std::endl;
            vipss_unit.out_normal_path_ = out_path + "/" + out_filename + "_out_normal";
            vipss_unit.out_surface_path_ = FLAGS_output;
            vipss_unit.out_debug_path_ = out_path + "/" + out_filename + "_debug.txt";
        } else {
            vipss_unit.out_normal_path_ = in_path + "/" + in_filename + "_out_normal";
            vipss_unit.out_surface_path_ = in_path + "/" + in_filename + "_out_surface";
            vipss_unit.out_debug_path_ = in_path + "/" + in_filename + "_debug.txt";
        }
        vipss_unit.user_lambda_ = FLAGS_lambda;
        vipss_unit.use_hrbf_surface_ = true;
        vipss_unit.volume_dim_ = FLAGS_volumeDim; 
        vipss_unit.axi_plane_ = AXI_PlANE::XYZ;
        vipss_unit.hrbf_type_ = HRBF_SURFACE_TYPE::LOCAL_HRBF_NN;
        // vu.hrbf_type_ = HRBF_SURFACE_TYPE::GLOBAL_HRBF;
        vipss_unit.Run();
    } else {
        std::cout << "There is no valid input point path !" << std::endl;
    }
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
