#include "sp_mat.h"


SpiMat ConvertToEigenSp(const arma::sp_imat& inmat)
{
    int cols = int(inmat.n_cols);
    int rows = int(inmat.n_rows);
    SpiMat new_mat(rows, cols);
    std::vector<TripletInt> eles;
    arma::sp_imat::const_iterator it     = inmat.begin();
    arma::sp_imat::const_iterator it_end = inmat.end();
    for(; it != it_end; ++it)
    {
        eles.push_back(TripletInt(it.row(), it.col(), *it));
    }
    new_mat.setFromTriplets(eles.begin(), eles.end());
    return new_mat;
} 

SpiMat ConvertToEigenSp(const arma::sp_umat& inmat)
{
    int cols = inmat.n_cols;
    int rows = inmat.n_rows;
    SpiMat new_mat(rows, cols);
    std::vector<TripletInt> eles;
    arma::sp_umat::const_iterator it     = inmat.begin();
    arma::sp_umat::const_iterator it_end = inmat.end();
    for(; it != it_end; ++it)
    {
        eles.push_back(TripletInt(it.row(), it.col(), *it));
    }
    new_mat.setFromTriplets(eles.begin(), eles.end());
    return new_mat;
} 

// template<typename InputIterator, typename DupFunctor>
// void set_from_triplets(const InputIterator& begin, const InputIterator& end, SparseMatrixType& mat, DupFunctor dup_func)
void set_from_triplets(std::vector<Triplet>& in_tris, SpMat& mat)
{

  Eigen::internal::scalar_sum_op<double, double> dup_func;

  enum { IsRowMajor = SpMat::IsRowMajor };
  typedef typename SpMat::Scalar Scalar;
  typedef typename SpMat::StorageIndex StorageIndex;
  std::cout << " rows : " << mat.rows()  << " cols : " << mat.cols() << std::endl;
//   SpMat trMat(mat.rows(),mat.cols());
  Eigen::SparseMatrix<Scalar,IsRowMajor?Eigen::ColMajor:Eigen::RowMajor,StorageIndex> trMat(mat.rows(),mat.cols());

  std::cout << " trMat.outerSize() " << trMat.outerSize() << std::endl;

//   if(begin!=end)
  {
    // pass 1: count the nnz per inner-vector
    typename SpMat::IndexVector wi(trMat.outerSize());
    wi.setZero();
    // for(InputIterator it(begin); it!=end; ++it)
    for(auto it = in_tris.begin(); it != in_tris.end(); ++it)
    {
      eigen_assert(it->row()>=0 && it->row()<mat.rows() && it->col()>=0 && it->col()<mat.cols());
      wi(IsRowMajor ? it->col() : it->row())++;
    }

    // pass 2: insert all the elements into trMat
    trMat.reserve(wi);
    // for(InputIterator it(begin); it!=end; ++it)
    for(auto it = in_tris.begin(); it != in_tris.end(); ++it)
    {
      trMat.insertBackUncompressed(it->row(),it->col()) = it->value();
    }
    // pass 3:
    trMat.collapseDuplicates(dup_func);
  }
  in_tris.clear();
  // pass 4: transposed copy -> implicit sorting
  mat = trMat;

  // pass 4: transposed copy -> implicit sorting
//   mat = trMat;
}



void removeDuplicatePairs(std::vector<MyTriplet>& vec, std::unordered_set<MyTriplet, PairHash>& seen) 
{    
    auto it = vec.begin();
    while (it != vec.end()) {
        MyTriplet pairValue = (*it); // Access the inner pair
        auto target = seen.find(pairValue);
        if (target != seen.end()) {
            const_cast<MyTriplet*>(&(*target))->add(it->value());
            // it = vec.erase(it); // Remove the element and update the iterator
        } else {
            seen.insert(pairValue);
            ++it;
        }
    }
}


void addSparseMatrixParallel(std::vector<SpMat>&diag_h_vec, std::vector<SpMat>&top_h_vec)
{
  int thread_num = 16;
  while(diag_h_vec.size() < 16)
  {
    diag_h_vec.push_back(SpMat(diag_h_vec[0].rows(), diag_h_vec[0].cols()));
    top_h_vec.push_back(SpMat(diag_h_vec[0].rows(), diag_h_vec[0].cols()));
  }
  // int base = 1;
  int level = 1;
  for(int base =1; base < thread_num; base *= 2)
  {
    int cur_th_num = thread_num / base;
    #pragma omp parallel for
    for(int tid = 0; tid < cur_th_num; tid ++)
    {
      if(tid % 2 == 0)
      {
        int cur_id = tid * base;
        diag_h_vec[cur_id] += diag_h_vec[cur_id + base];
        diag_h_vec[cur_id + base].resize(0,0);
        diag_h_vec[cur_id + base].data().squeeze();
      } else {
        int cur_id = (tid-1) * base;
        top_h_vec[cur_id] += top_h_vec[cur_id + base];
         top_h_vec[cur_id + base].resize(0,0);
          top_h_vec[cur_id + base].data().squeeze();
      }
    }
  }

  
      // #pragma omp parallel for
      // for(int tid = 0; tid < 16; tid ++)
      // {
      //     // int cur_id = tid/2;
      //     if(tid % 2 == 0)
      //     {
      //         diag_h_vec[tid] += diag_h_vec[tid + 1];
      //     } else {
      //         top_h_vec[tid-1] += top_h_vec[tid];
      //     }
      // }
      // #pragma omp parallel for
      // for(int tid = 0; tid < 8; tid ++)
      // {
      //     if(tid % 2 == 0)
      //     {
      //         int cur_id = tid * 2;
      //         diag_h_vec[cur_id] += diag_h_vec[cur_id + 2];
      //     } else {
      //         int cur_id = (tid-1) * 2;
      //         top_h_vec[cur_id] += top_h_vec[cur_id + 2];
      //     }
      // }
      // #pragma omp parallel for
      // for(int tid = 0; tid < 4; tid ++)
      // {
      //     if(tid % 2 == 0)
      //     {
      //         int cur_id = tid * 4;
      //         diag_h_vec[cur_id] += diag_h_vec[cur_id + 4];
      //     } else {
      //         int cur_id = (tid-1) * 4;
      //         top_h_vec[cur_id] += top_h_vec[cur_id + 4];
      //     }
      // }
      // #pragma omp parallel for
      // for(int tid = 0; tid < 2; tid ++)
      // {
      //     if(tid % 2 == 0)
      //     {
      //         int cur_id = tid * 8;
      //         diag_h_vec[cur_id] += diag_h_vec[cur_id + 8];
      //     } else {
      //         int cur_id = (tid-1) * 8;
      //         top_h_vec[cur_id] += top_h_vec[cur_id + 8];
      //     }
      // }
}