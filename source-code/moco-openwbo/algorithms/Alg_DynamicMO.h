#ifndef ALG_DYNAMIC_H
#define ALG_DYNAMIC_H

class DynamicMO{
public:
  virtual void increment() = 0;
  virtual void checkSols() = 0;
  virtual void checkLower() = 0;
};

#endif
