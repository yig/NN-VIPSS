#include "test_timing.h"
#include <cblas.h>

void test_vipss()
{
    int num_threads = openblas_get_num_threads();
    openblas_set_num_threads(1);
    
    printf("OpenBLAS is using %d threads.\n", num_threads);
    for(int i = 1; i <= 10; ++i)
    {

        std::string path = "../../data/kitten/kitten_" + std::to_string(i) + "0k.xyz";
        std::cout << " path " << path << std::endl;
        test_vipss_timing::test_vipss_timing(path);
        
    }
    
}


int main()
{
    test_vipss();
    return 0;
}