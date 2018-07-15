#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <iostream>

#include "qvoter.h"

using namespace boost::numeric::odeint;
using namespace std;

int main(int ac, char *av[])
{
  const size_t kmax = 50;
  const unsigned Q = 8;
  QVoter ode(kmax, Q);
  ode.set_p(0.07);

  double_v x(3 * kmax);
  generate(begin(x), end(x), []() { return 0.1; });
// x[3*(kmax-1)]=0.5;
// x[3*(kmax-1)+1]=0.5;
// x[3*(kmax-1)+2]=0.5;

  double dt = 0.1, T = 100;

  vector<double> m;
  m.reserve(T / dt);

  auto obs = [&](const double_v &x_, double) {
    double r = 0;
    for (size_t i = 0; i < kmax; ++i)
      r += x_[3 * i] * ode.P(i + 1);
    m.push_back(r);
  };

  integrate_const(runge_kutta4<double_v>(), ref(ode), x, 0.0, T, dt, obs);

  for (auto y : m)
    cout << y << '\n';
  cout << endl;
}