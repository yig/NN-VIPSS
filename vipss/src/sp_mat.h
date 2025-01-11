#pragma once
#include <Eigen/Sparse>
#include <armadillo>
#include <unordered_set>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix<int> SpiMat;
typedef Eigen::Triplet<double> Triplet;
typedef Eigen::Triplet<int> TripletInt;
typedef Eigen::SparseVector<int> SpiVec;

enum { IsRowMajor = SpMat::IsRowMajor };
// bool IsRowMajor = SpMatD::IsRowMajor;
typedef typename SpMat::Scalar Scalar;
typedef typename SpMat::StorageIndex StorageIndex;
// SpMatD:

// typedef Eigen::SparseMatrix<double, Eigen::ColMajor,int> SpMat;
// typedef Eigen::SparseMatrix<Scalar, Eigen::ColMajor,long> SpMat;
typedef Eigen::SparseMatrix<Scalar,IsRowMajor?Eigen::ColMajor:Eigen::RowMajor, long> SpMatL;


class MyTriplet
{
public:
  MyTriplet() : m_row(0), m_col(0), m_value(0) {}

  MyTriplet(const StorageIndex& i, const StorageIndex& j, const Scalar& v = Scalar(0))
    : m_row(i), m_col(j), m_value(v)
  {}
  /** \returns the row index of the element */
  const StorageIndex& row() const { return m_row; }

  /** \returns the column index of the element */
  const StorageIndex& col() const { return m_col; }

  /** \returns the value of the element */
  const Scalar& value() const { return m_value; }
  
  void add(const Scalar v) { m_value += v;} 
  bool operator==(const MyTriplet& other) const {
        return m_row == other.m_row && m_col == other.m_col;
    }
//   void add( Scalar v) { m_value += v;} 
protected:
  StorageIndex m_row, m_col;
  Scalar m_value;
};

struct PairHash {
    // Hash function for std::pair<int, int>
    size_t operator()(const MyTriplet p) const {
        return std::hash<int>()(p.row()) ^ std::hash<int>()(p.col());
    }
};

void removeDuplicatePairs(std::vector<MyTriplet>& vec, std::unordered_set<MyTriplet, PairHash>& seen); 



SpiMat ConvertToEigenSp(const arma::sp_imat& inmat);
SpiMat ConvertToEigenSp(const arma::sp_umat& inmat);
void set_from_triplets(const std::vector<Triplet>& in_tris, SpMat& trMat);
void addSparseMatrixParallel(std::vector<SpMat>&diag_h_vec, std::vector<SpMat>&top_h_vec);