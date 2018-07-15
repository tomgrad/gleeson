#ifndef QVOTER_H
#define QVOTER_H

#include <cmath>
#include <vector>
#include <iostream>
using std::pow;

typedef std::vector<double> double_v;
class QVoter
{
public:
  QVoter(const size_t k_, const unsigned q_) : q(q_), kmax(k_)
  {
    FB.resize(kmax);
    RB.resize(kmax);
    mFB.resize(kmax);
    mRB.resize(kmax);
  };

  void operator()(const double_v &x, double_v &dxdt,
                  const double /* t */)
  {
    coefs(x);
    for (size_t i = 0; i < kmax; ++i)
    {
      auto const rk = x[3 * i];
      auto const pk = x[3 * i + 1];
      auto const qk = x[3 * i + 2];
      double k = 1.0 + i;

      dxdt[3 * i] = -rk * RB[i] + (1 - rk) * FB[i];
      // auto tmp1=pk * FB[i];
      // auto tmp2=- pk * rk / (1 - rk) * RB[i];
      // auto tmp3=- mFB[i] / k;
      // auto tmp4=rk / (1 - rk) / k * mRB[i];
      
      // dxdt[3 * i + 1] = tmp1+tmp2+tmp3+tmp4 + bs * (1 - pk) - gs * pk;
      dxdt[3 * i + 1] = pk * FB[i] - pk * rk / (1 - rk) * RB[i] - mFB[i] / k + rk / (1 - rk) / k * mRB[i] + bs * (1 - pk) - gs * pk;

      dxdt[3 * i + 2] = qk * RB[i] - qk / rk * (1 - rk) * FB[i] - mRB[i] / k + (1 - rk) / rk / k * mFB[i] + bi * (1 - qk) - gi * qk;
    }
  }

  double P(const unsigned k)
  {
    // random regular
    if (k == kmax)
      return 1;
    else
      return 0.0;

    // return 1.0/kmax;
  }
  void set_p(double p_) { p = p_; }

private:
  size_t kmax;
  unsigned q;
  double p=0.1;
  double bs, bi, gs, gi;
  std::vector<double> FB, RB, mFB, mRB;

  void coefs(const double_v &x)
  {
//    double sFB = 0, sRB = 0, smFB = 0, smRB = 0;

    bs = bi = gs = gi = 0;
    double mbs = 0, mbi = 0, mgs = 0, mgi = 0;
    for (size_t i = 0; i < kmax; ++i)
    {
      auto const rk = x[3 * i];
      auto const pk = x[3 * i + 1];
      auto const qk = x[3 * i + 2];
      const double k = i + 1.0;

      FB[i] = (1 - p) * pow(pk, q) + p / 2;
      RB[i] = (1 - p) * pow(1 - qk, q) + p / 2;
      mFB[i] = (1 - p) * (pk * (k - q) + q) * pow(pk, q) + p * k * pk / 2;
      mRB[i] = (1 - p) * qk * (k - q) * pow(1 - qk, q) + p * k * qk / 2;

      bs += P(k) * (1 - rk) * (k * FB[i] - mFB[i]);
      bi += P(k) * (1 - rk) * mFB[i];
      gs += P(k) * rk * (k * RB[i] - mRB[i]);
      gi += P(k) * rk * mRB[i];

      mbs += P(k) * (1 - rk) * k * (1 - pk);
      mbi += P(k) * (1 - rk) * k * pk;
      mgs += P(k) * rk * k * (1 - qk);
      mgi += P(k) * rk * k * qk;
    }
    bs /= mbs;
    bi /= mbi;
    gs /= mgs;
    gi /= mgi;
    // std::cout << mbi << ' ' << mgi << ' ' << mbs << ' ' << mgs << std::endl;
  }
};

#endif // QVOTER_H