#ifndef QVOTER_H
#define QVOTER_H

#include <vector>

typedef std::vector<double> container_type;

class QVoter {
public:
  QVoter(unsigned q_):q(q_){};

  void operator()(const container_type &x, container_type &dxdt,
                  const double /* t */) {

    for (size_t i = 0; i < x.size(); i += 3) {
      dxdt[i] = x[i + 1];
      dxdt[i + 1] = -x[i];
      dxdt[i + 2] = 0;
    }
  }
  void set_p(double p_){p=p_;}

private:
  unsigned q;
  double p;
  double bs, bi, gs, gi;
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