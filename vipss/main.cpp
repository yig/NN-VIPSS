#include <iostream>
//#include <unistd.h>
#include "src/rbfcore.h"
#include "src/readers.h"
#include "src/voronoi_gen.h"

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



void SplitPath(const std::string& fullfilename,std::string &filepath);
void SplitFileName (const std::string& fullfilename,std::string &filepath,std::string &filename,std::string &extname);
RBF_Paras Set_RBF_PARA();
int main(int argc, char** argv)
{

    test_voronoi();
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
