#ifndef KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_

class MotherCuts {

 public:
  MotherCuts() = default;
  MotherCuts(const MotherCuts&) = default;
  MotherCuts(MotherCuts&&) = default;
  MotherCuts& operator=(MotherCuts&&) = default;
  MotherCuts& operator=(const MotherCuts&) = default;
  ~MotherCuts() = default;

 protected:
  float distance_{-1.f};
  float chi2_geo_{-1.f};
  float ldl_{-1.f};
  float chi2_topo_{-1.f};
};

#endif //KFPARTICLESIMPLE_KFSIMPLE_MOTHERCUTS_HPP_
