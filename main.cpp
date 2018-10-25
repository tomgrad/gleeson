#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include "../modules/cli_parser.h"
#include <cmath>
#include "qvoter.h"

using namespace boost::numeric::odeint;
using namespace std;

int main(int ac, char *av[])
{
  cli_parser p(ac, av);

struct
{
  bool steady;
  unsigned kmax;
  unsigned kmin;
  double gamma;
  unsigned Q;
  double P;
  double T;
  double init;
} par;

par.steady = p.get_flag("steady");
par.kmin = p.get<unsigned>("kmin");
par.kmax = p.get<unsigned>("kmax");
par.gamma = p.get<double>("gamma");
par.Q = p.get<unsigned>("q");
par.P = p.get<double>("p");
par.T = p.get<double>("T", 10);
par.init = p.get<double>("init", 0.5);


  // const auto steady = p.get_flag("steady");
  // const auto kmax = p.get<unsigned>("kmax", 1000);
  // const auto kmin = p.get<unsigned>("kmin", 20);
  // const auto gamma = p.get<double>("gamma", 3);
  // const auto Q = p.get<unsigned>("q", 4);
  // auto P = p.get<double>("p", 0);
  // auto T = p.get<double>("T", 10);
  // auto init = p.get<double>("init", 0.5);


  QVoter ode(par.kmin, par.kmax, par.Q);
  ode.set_p(par.P);
  ode.set_gamma(par.gamma);

  double_v x(3 * par.kmax);
  double dt = 0.1;

  auto set_ic = [&](double x_) {
    generate(begin(x), end(x), [&]() { return x_; });
  };

  auto rho = [&](const double_v &x_) {
    double r = 0;
    for (size_t i = par.kmin - 1; i < par.kmax; ++i)
      r += x_[3 * i] * ode.P(i + 1);
    return r;
  };

  if (par.steady)
  {
    double delta = 0.0001;
    double px = 0;
    double T = 100;

    double xs1 = 0.999;
    set_ic(xs1);
    while (fabs(px - xs1) > delta)
    {
      px = xs1;
      integrate_const(runge_kutta4<double_v>(), ref(ode), x, 0.0, T, dt);
      xs1 = rho(x);
    }

    cout << par.P << "\t" << xs1 << "\t1" << endl;

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
      cout << par.P << "\t"
           << "0.5"
           << "\t0" << endl;
    }
    else
    {
      double ux;
      cout << par.P << "\t" << xs2 << "\t1" << endl;

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

        cout << par.P << "\t" << xs2 << "\t0" << endl;
      }
    }
  }

  else
  {
    generate(begin(x), end(x), [&]() { return par.init; });

    vector<double> m;

    m.reserve(par.T / dt);

    auto obs = [&](const double_v &x_, double) {
      m.push_back(rho(x_));
    };

    integrate_const(runge_kutta4<double_v>(), ref(ode), x, 0.0, par.T, dt, obs);

    for (auto y : m)
      cout << y << '\n';
    cout << endl;
  }
}