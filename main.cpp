#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <cli_parser.h>
#include <cmath>
#include "qvoter.h"

using namespace boost::numeric::odeint;
using namespace std;

int main(int ac, char *av[])
{
  cli_parser p(ac, av);

  const auto stable = p.get_flag("stable");
  const auto unstable = p.get_flag("unstable");

  const auto kmax = p.get<unsigned>("kmax", 1000);
  const auto kmin = p.get<unsigned>("kmin", 20);
  const auto gamma = p.get<double>("gamma", 3);
  const auto Q = p.get<unsigned>("q", 4);
  auto P = p.get<double>("p", 0);
  auto T = p.get<double>("T", 10);
  auto init = p.get<double>("init", 0.5);
  QVoter ode(kmin, kmax, Q);
  ode.set_p(P);
  ode.set_gamma(gamma);

  double_v x(3 * kmax);
  double dt = 0.1;

  auto set_ic = [&](double x_) {
    generate(begin(x), end(x), [&]() { return x_; });
  };

  auto rho = [&](const double_v &x_) {
    double r = 0;
    for (size_t i = kmin - 1; i < kmax; ++i)
      r += x_[3 * i] * ode.P(i + 1);
    return r;
  };

  if (stable)
  {
    double delta = 0.001;
    double px = 0;
    T = 100;

    set_ic(0.999);
    while (fabs(px - rho(x)) > delta)
    {
      px = rho(x);
      integrate_const(runge_kutta4<double_v>(), ref(ode), x, 0.0, T, dt);
    }
    cout << P << "\t" << px << endl;

    set_ic(0.49);
    px = 0;
    while (fabs(px - rho(x)) > delta)
    {
      px = rho(x);
      integrate_const(runge_kutta4<double_v>(), ref(ode), x, 0.0, T, dt);
    }
    if (px > 0.49)
      cout << P << "\t" << px << endl;
  }

  else if (unstable)
  {
  }

  else
  {
    generate(begin(x), end(x), [&]() { return init; });

    vector<double> m;

    m.reserve(T / dt);

    auto obs = [&](const double_v &x_, double) {
      m.push_back(rho(x_));
    };

    integrate_const(runge_kutta4<double_v>(), ref(ode), x, 0.0, T, dt, obs);

    for (auto y : m)
      cout << y << '\n';
    cout << endl;
  }
}