#ifndef KFPARTICLESIMPLE_SRC_NONLINEARCUTBASE_HPP_
#define KFPARTICLESIMPLE_SRC_NONLINEARCUTBASE_HPP_

#include "Constants.hpp"

class NonLinearCutBase{

 public:

  virtual bool ApplyCut(const SelectionValues& values) = 0;

 protected:

};

#endif //KFPARTICLESIMPLE_SRC_NONLINEARCUTBASE_HPP_
