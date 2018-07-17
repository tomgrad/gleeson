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

  const auto steady = p.get_flag("steady");
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

  if (steady)
  {
    double delta = 0.0001;
    double px = 0;
    T = 100;

    double xs1 = 0.999;
    set_ic(xs1);
    while (fabs(px - xs1) > delta)
    {
      px = xs1;
      integrate_const(runge_kutta4<double_v>(), ref(ode), x, 0.0, T, dt);
      xs1 = rho(x);
    }

    cout << P << "\t" << xs1 << "\t1" << endl;

    double xs2 = 0.499;
    set_ic(xs2);
    px = 0;
    while (fabs(px - xs2) > delta)
    {
      px = xs2;
      integrate_const(runge_kutta4<double_v>(), ref(ode), x, 0.0, T, dt);
      xs2 = rho(x);
    }
    if (xs2 < 0.499)
    {
      xs2 = xs1;
      cout << P << "\t"
           << "0.5"
           << "\t0" << endl;
    }
    else
    {
      double ux;
      cout << P << "\t" << xs2 << "\t1" << endl;

      if (xs1 - xs2 > delta)
      {

        while (xs1 - xs2 > delta)
        {
          ux = (xs1 + xs2) / 2;
          set_ic(ux);
          integrate_const(runge_kutta4<double_v>(), ref(ode), x, 0.0, T, dt);
          if (rho(x) < ux)
            xs2 = ux;
          else
            xs1 = ux;
        }

        cout << P << "\t" << xs2 << "\t0" << endl;
      }
    }
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