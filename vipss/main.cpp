#include <iostream>
//#include <unistd.h>
#include "src/rbfcore.h"
#include "src/readers.h"
#include "src/voronoi_gen.h"
#include "src/local_vipss.hpp"
#include "src/vipss_unit.hpp"
#include "src/sample.h"
#include "src/kdtree.hpp"
#include "src/incremental_vipss.h"

typedef std::chrono::high_resolution_clock Clock;
using namespace std;

void test_sphere_creation()
{
    double cx = 0, cy = 0, cz = 0;
    double radius = 1.0;
    int pt_num = 100;
    auto spts = CreateSpherePoints(cx, cy, cz, radius, pt_num);
    std::string out_dir = "./spheres";
    writePLYFile(out_dir, spts);
    // writeXYZ(out_dir, spts);
}

void test_voronoi()
{
    std::string path = "../data/doghead200/doghead200";
    // std::string path = "../../data/doghead_wire/doghead_wire";
    VoronoiGen v_gen;
    v_gen.filename_ = "doghead200";
    v_gen.out_dir_ = "../data/doghead200/";
    v_gen.loadData(path);
    v_gen.Run();
}

void test_local_vipss()
{
    std::string data_dir =  "../../data/";
    
    LocalVipss l_vp;
    l_vp.filename_ = "arma";
    // l_vp.filename_ = "planck";
    l_vp.out_dir_ = data_dir + l_vp.filename_ + "/";
    std::string path = data_dir + l_vp.filename_ + "/" + l_vp.filename_ + ".ply";

    l_vp.angle_threshold_ = 25;
    l_vp.user_lambda_ = 0.0;
    l_vp.max_iter_ = 30;
    l_vp.use_hrbf_surface_ = false;
    l_vp.volume_dim_ = 100;
    l_vp.Init(path);
    l_vp.Run();
}

void test_vipss_unit()
{
    VIPSSUnit vu;
    vu.data_dir_ = "../../data/";
    // vu.file_name_ = "arma_100k";
    vu.file_name_ = "walrus";
    vu.user_lambda_ = 0.005;
    // vu.init_lambda_ = 0.0001;
    vu.init_with_cluster_merge_ = false;
    vu.merge_angle_ = 45;
    // vu.init_lambda_ = 0;
    vu.use_hrbf_surface_ = true;
    vu.volume_dim_ = 100;
    vu.axi_plane_ = AXI_PlANE::XYZ;
    vu.hrbf_type_ = HRBF_SURFACE_TYPE::LOCAL_HRBF_NN;
    // vu.hrbf_type_ = HRBF_SURFACE_TYPE::GLOBAL_HRBF;
    vu.Run(); 
    return;
}

void test_incre_vipss()
{
    IncreVipss incre_vipss;
    incre_vipss.lambda_ = 0.001;
    incre_vipss.input_path_ = "/home/jjxia/Documents/projects/VIPSS/data/wireframes/spring/input.xyz";
    // incre_vipss.input_path_ = "/home/jjxia/Documents/projects/VIPSS/data/fertility/input.xyz";
    // incre_vipss.input_path_ = "/home/jjxia/Documents/projects/VIPSS/data/surfaces_500/hand/input.xyz";
    // incre_vipss.input_path_ = "/home/jjxia/Documents/projects/VIPSS/data/wireframes/doghead/input.xyz";;

    incre_vipss.out_dir_ = "../out_dir/";
    incre_vipss.residual_beta_ = 1.0;
    incre_vipss.init_sample_num_ = 50;
    incre_vipss.Init();
    incre_vipss.Run();
}

arma::mat generate_mat(size_t n)
{
    arma::mat h_mat;
    h_mat.randu(n, n);

    return h_mat.t(); 
}

void test_arma_mat_inv()
{
    std::vector<int> mat_sizes{200, 400, 600, 800};

    for(size_t i = 0; i < mat_sizes.size(); ++i)
    {
        int cur_size = mat_sizes[i];

        arma::mat new_mat = generate_mat(cur_size);
        auto t0 = Clock::now();
        for(size_t j =0; j < 100; ++j)
        {
            arma::mat inv_mat = arma::inv(new_mat);
        }
        auto t1 = Clock::now();
        double inv_time = std::chrono::nanoseconds(t1 - t0).count() / 1e9;
        printf("mat size %d , running time %f \n",cur_size, inv_time/ 100.0);

    }
}

void SplitPath(const std::string& fullfilename,std::string &filepath);
void SplitFileName (const std::string& fullfilename,std::string &filepath,std::string &filename,std::string &extname);
RBF_Paras Set_RBF_PARA();
int main(int argc, char** argv)
{
    // test_voronoi();
    // test_pt_vipss();
    // test_local_vipss();
    // test_sphere_creation();
    test_vipss_unit();
    // test_incre_vipss();
    // test_arma_mat_inv();
    return 0;
    cout << argc << endl;
    string infilename;
    string outpath, pcname, ext, inpath;

    int n_voxel_line = 100;
    double user_lambda = 0;

    bool is_surfacing = false;
    bool is_outputtime = false;
    
    return 0;
}


RBF_Paras Set_RBF_PARA(){

    RBF_Paras para;
    RBF_InitMethod initmethod = Lamnbda_Search;

    RBF_Kernal Kernal = XCube;
    int polyDeg = 1;
    double sigma = 0.9;
    double rangevalue = 0.001;

    para.Kernal = Kernal;para.polyDeg = polyDeg;para.sigma = sigma;para.rangevalue = rangevalue;
    para.Hermite_weight_smoothness = 0.0;
    para.Hermite_ls_weight = 0;
    para.Hermite_designcurve_weight = 00.0;
    para.Method = RBF_METHOD::Hermite_UnitNormal;
    para.InitMethod = initmethod;
    para.user_lamnbda = 0;
    para.isusesparse = false;
    return para;
}


inline void SplitFileName (const std::string& fullfilename,std::string &filepath,std::string &filename,std::string &extname) {
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
