#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rank_matrix_cpp(NumericMatrix mat,
                              bool descending = true,
                              bool use_abs = false) {
  
  int nrows = mat.nrow();
  int ncols = mat.ncol();
  
  NumericMatrix ranks(nrows, ncols);
  
  auto get_val = [use_abs](double x) {
    return use_abs ? std::fabs(x) : x;
  };
  
  auto comp = [descending](double va, double vb) {
    if(descending) return va > vb;
    return va < vb;
  };
  
  for(int i = 0; i < nrows; i++) {
    
    std::vector<double> row_vec(ncols);
    
    for(int j = 0; j < ncols; j++) {
      row_vec[j] = mat(i, j);
    }
    
    std::vector<size_t> idx(ncols);
    std::iota(idx.begin(), idx.end(), 0);
    
    std::sort(
      idx.begin(),
      idx.end(),
      [&get_val, &comp, &row_vec](size_t a, size_t b) {
        double va = get_val(row_vec[a]);
        double vb = get_val(row_vec[b]);
        
        if(va == vb) return a < b;
        
        return comp(va, vb);
      }
    );
    
    size_t k = 0;
    
    while(k < static_cast<size_t>(ncols)) {
      
      size_t start = k;
      double val = get_val(row_vec[idx[k]]);
      
      while(k < static_cast<size_t>(ncols) &&
            get_val(row_vec[idx[k]]) == val) {
        ++k;
      }
      
      size_t end = k - 1;
      
      double avg_rank =
        (static_cast<double>(start + 1) + static_cast<double>(end + 1)) / 2.0;
      
      for(size_t j = start; j <= end; ++j) {
        size_t id = idx[j];
        ranks(i, id) = avg_rank;
      }
    }
  }
  
  return ranks;
}


// [[Rcpp::export]]
List flatten_cor_matrix_cpp(NumericMatrix cormat,
                            Nullable<NumericMatrix> mrmat = R_NilValue,
                            Nullable<NumericMatrix> pmat = R_NilValue,
                            Nullable<NumericMatrix> padjmat = R_NilValue,
                            CharacterVector row_names = CharacterVector::create()) {
  
  int m = cormat.nrow();
  
  bool has_mr = mrmat.isNotNull();
  bool has_p = pmat.isNotNull();
  bool has_pa = padjmat.isNotNull();
  bool has_names = row_names.size() == m;
  
  size_t num = static_cast<size_t>(m) * static_cast<size_t>(m - 1LL) / 2;
  
  CharacterVector rows;
  CharacterVector cols;
  
  if(has_names) {
    rows = CharacterVector(num);
    cols = CharacterVector(num);
  }
  
  NumericVector cors(num);
  NumericVector mrs;
  NumericVector ps;
  NumericVector pas;
  
  if(has_mr) mrs = NumericVector(num);
  if(has_p) ps = NumericVector(num);
  if(has_pa) pas = NumericVector(num);
  
  NumericMatrix mr = has_mr ? NumericMatrix(mrmat) : NumericMatrix();
  NumericMatrix p = has_p ? NumericMatrix(pmat) : NumericMatrix();
  NumericMatrix pa = has_pa ? NumericMatrix(padjmat) : NumericMatrix();
  
  size_t cnt = 0;
  
  for(int i = 0; i < m - 1; i++) {
    for(int j = i + 1; j < m; j++) {
      
      if(has_names) {
        rows[cnt] = row_names[i];
        cols[cnt] = row_names[j];
      }
      
      cors[cnt] = cormat(i, j);
      
      if(has_mr) mrs[cnt] = mr(i, j);
      if(has_p) ps[cnt] = p(i, j);
      if(has_pa) pas[cnt] = pa(i, j);
      
      cnt++;
    }
  }
  
  List res;
  
  if(has_names) {
    res["row"] = rows;
    res["column"] = cols;
  }
  
  res["cor"] = cors;
  
  if(has_mr) res["mr"] = mrs;
  if(has_p) res["p"] = ps;
  if(has_pa) res["p.adj"] = pas;
  
  return res;
}


// [[Rcpp::export]]
List flatten_rect_cor_matrix_cpp(NumericMatrix cormat,
                                 Nullable<NumericMatrix> mrmat = R_NilValue,
                                 Nullable<NumericMatrix> pmat = R_NilValue,
                                 Nullable<NumericMatrix> padjmat = R_NilValue,
                                 CharacterVector row_names = CharacterVector::create(),
                                 CharacterVector col_names = CharacterVector::create()) {
  
  int nr = cormat.nrow();
  int nc = cormat.ncol();
  
  bool has_mr = mrmat.isNotNull();
  bool has_p = pmat.isNotNull();
  bool has_pa = padjmat.isNotNull();
  bool has_row_names = row_names.size() == nr;
  bool has_col_names = col_names.size() == nc;
  
  size_t num = static_cast<size_t>(nr) * static_cast<size_t>(nc);
  
  CharacterVector rows;
  CharacterVector cols;
  
  if(has_row_names && has_col_names) {
    rows = CharacterVector(num);
    cols = CharacterVector(num);
  }
  
  NumericVector cors(num);
  NumericVector mrs;
  NumericVector ps;
  NumericVector pas;
  
  if(has_mr) mrs = NumericVector(num);
  if(has_p) ps = NumericVector(num);
  if(has_pa) pas = NumericVector(num);
  
  NumericMatrix mr = has_mr ? NumericMatrix(mrmat) : NumericMatrix();
  NumericMatrix p = has_p ? NumericMatrix(pmat) : NumericMatrix();
  NumericMatrix pa = has_pa ? NumericMatrix(padjmat) : NumericMatrix();
  
  size_t cnt = 0;
  
  for(int i = 0; i < nr; i++) {
    for(int j = 0; j < nc; j++) {
      
      if(has_row_names && has_col_names) {
        rows[cnt] = row_names[i];
        cols[cnt] = col_names[j];
      }
      
      cors[cnt] = cormat(i, j);
      
      if(has_mr) mrs[cnt] = mr(i, j);
      if(has_p) ps[cnt] = p(i, j);
      if(has_pa) pas[cnt] = pa(i, j);
      
      cnt++;
    }
  }
  
  List res;
  
  if(has_row_names && has_col_names) {
    res["row"] = rows;
    res["column"] = cols;
  }
  
  res["cor"] = cors;
  
  if(has_mr) res["mr"] = mrs;
  if(has_p) res["p"] = ps;
  if(has_pa) res["p.adj"] = pas;
  
  return res;
}