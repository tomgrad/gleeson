#include <iostream>
#include <algorithm>
#include <boost/numeric/odeint.hpp>


#include "qvoter.h"

using namespace boost::numeric::odeint;
using namespace std;

int main(int ac, char *av[]) {
const size_t N=3;
  QVoter ode(4);
  ode.set_p(0.1);
  Observer obs;
  container_type x(3*N);

generate(begin(x), end(x), [](){return 0.1;});

  double dt = 0.01, T = 10;

  // integrate_const(runge_kutta4<container_type>(), ref(k), x, 0.0, Ttrans,
  // dt); // transient
  integrate_const(runge_kutta4<container_type>(), ref(ode), x, 0.0, T, dt,
                  ref(obs));


for (auto y:obs.x)
  cout << y << '\n';
cout << endl;

}