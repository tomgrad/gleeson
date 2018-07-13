#ifndef QVOTER_H
#define QVOTER_H

#include <cmath>
#include <vector>

using std::pow;

typedef std::vector<double> container_type;

class QVoter {
public:
  QVoter(const size_t k_, const unsigned q_) : q(q_), kmax(k_) {
    FB.resize(kmax);
    RB.resize(kmax);
    mFB.resize(kmax);
    mRB.resize(kmax);
  };

  void operator()(const container_type &x, container_type &dxdt,
                  const double /* t */) {
    coefs(x);
    for (size_t i = 0; i < kmax; ++i) {
      auto const &rk = x[3 * i];
      auto const &pk = x[3 * i + 1];
      auto const &qk = x[3 * i + 2];
      auto const k = i + 1;

      dxdt[3 * i] = x[3 * i + 1];
      dxdt[3 * i + 1] = -x[3 * i];
      dxdt[3 * i + 2] = 0;
    }
  }
  void set_p(double p_) { p = p_; }

private:
  size_t kmax;
  unsigned q;
  double p;
  double bs, bi, gs, gi;
  std::vector<double> FB, RB, mFB, mRB;

  void coefs(const container_type &x) {
    for (size_t i = 0; i < kmax; ++i) {

      auto const &pk = x[3 * i + 1];
      auto const &qk = x[3 * i + 2];
      auto const &k = i + 1;

      FB[i] = (1 - p) * pow(pk, q) + p / 2;
      RB[i] = (1 - p) * pow((1 - qk), q) + p / 2;
      mFB[i] = (1 - p) * (pk * (k - q) + q) * pow(pk, q) + p * k * pk / 2;
      mRB[i] = (1 - p) * qk * (k - q) * pow(1 - qk, q) + p * k * qk / 2;
    }
  }
};

struct Observer {
  size_t count = 0;
  std::vector<double> x;
  //  template <class State> void operator()(const State &x, double t) {
  void operator()(const container_type &x_, double) {
    ++count;
    double m = 0;
    for (size_t i = 0; i < x_.size(); i += 3) {
      m += x_[i];
    }
    x.push_back(m);
  }
};

#endif // QVOTER_H