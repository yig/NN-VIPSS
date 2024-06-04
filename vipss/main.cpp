#include <iostream>
//#include <unistd.h>
#include "src/rbfcore.h"
#include "src/readers.h"
#include "src/voronoi_gen.h"
#include "src/pt_vipss.h"
#include "src/local_vipss.hpp"

using namespace std;

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
    l_vp.filename_ = "walrus";
    // l_vp.filename_ = "planck";
    l_vp.out_dir_ = data_dir + l_vp.filename_ + "/";
    std::string path = data_dir + l_vp.filename_ + "/" + l_vp.filename_ + ".ply";

    l_vp.angle_threshold_ = 25;
    l_vp.user_lambda_ = 0.003;
    l_vp.max_iter_ = 30;
    l_vp.use_hrbf_surface_ = true;
    l_vp.Init(path);
    l_vp.Run();
}

void test_pt_vipss()
{
    std::string path = "../data/doghead200/doghead200";
    // std::string path = "../../data/doghead_wire/doghead_wire";
    // VoronoiGen v_gen;
    // v_gen.filename_ = "doghead200";
    // v_gen.out_dir_ = "../data/doghead200/";
    // v_gen.loadData(path);
    // v_gen.Run();

    PtVipss p_vipss;
    printf("start to init p_vipss! \n");
    p_vipss.Init(path);
    p_vipss.out_name_ = "doghead200";
    p_vipss.out_dir_ = "../data/doghead200/";

    p_vipss.Run();
}


void SplitPath(const std::string& fullfilename,std::string &filepath);
void SplitFileName (const std::string& fullfilename,std::string &filepath,std::string &filename,std::string &extname);
RBF_Paras Set_RBF_PARA();
int main(int argc, char** argv)
{

    // test_voronoi();
    // test_pt_vipss();
    test_local_vipss();
    return 0;
    cout << argc << endl;



    string infilename;
    string outpath, pcname, ext, inpath;

    int n_voxel_line = 100;

    double user_lambda = 0;

    bool is_surfacing = false;
    bool is_outputtime = false;

    //int c;
    
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
