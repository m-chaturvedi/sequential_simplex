#ifndef SRC_HELPERS_H_
#define SRC_HELPERS_H_

#include <armadillo>
#include <string>
#include <climits>
#include <boost/assert.hpp>

#define assert_msg BOOST_ASSERT_MSG
#define ABS(x) fabs(x)

typedef arma::mat MAT;
typedef arma::vec VEC;
typedef arma::uvec UVEC;
typedef arma::rowvec ROWVEC;

typedef bool BOOL;
typedef double FL;
typedef int INTEGER;
// uword is a type of index
typedef arma::uword IND_TYPE;

const int INTEGER_MAX=INT_MAX;
const FL EPS = 1e-16;

// CONFIG
const FL MACF_TEST_PERCENT = 0.95;

const std::string IP_DIR = "res/input/";
const std::string OP_DIR = "res/output/";

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

const BOOL WITH_LINEAR_REG = false;
extern INTEGER MAX_ITERS;
const BOOL EBAR_ZERO = false;

#endif
