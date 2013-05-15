#-------------------------------------------------------------------------------
# Common functions useful for any simulation study with HLM residuals.
# 
# March 2013
# Adam Loy
#-------------------------------------------------------------------------------

### Matrix multiplication
prodCpp <- '
using Eigen::Map;
using Eigen::MatrixXd;
const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
const Map<MatrixXd> B(as<Map<MatrixXd> >(BB));
return wrap( A * B );
'

cxxprod <- cxxfunction(signature(AA = "matrix", BB = "matrix"),
                       prodCpp, plugin = "RcppEigen")

### Crossproduct of a singular matrix
crossprodCpp <- '
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
const int m(A.rows()), n(A.cols());
MatrixXd AtA(MatrixXd(n, n).setZero().
selfadjointView<Lower>().rankUpdate(A.adjoint()));
return wrap(AtA);
'

cxxcrossprod <- cxxfunction(signature(AA = "matrix"),
                            crossprodCpp, plugin = "RcppEigen")


### Matrix addition/subtraction
subCpp <- '
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
const MapMatd B(as<MapMatd>(BB));
const MapMatd C(as<MapMatd>(CC));
return wrap(B - C);
'

cxxmatsub <- cxxfunction(signature(BB = "matrix", CC = "matrix"),
                         subCpp, plugin = "RcppEigen")
