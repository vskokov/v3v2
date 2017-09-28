#ifndef DataStack_h
#define DataStack_h

#include <cstddef>
#include <cassert>
#include <cmath>


class DataStack
{  
 public:
  DataStack();
  ~DataStack();

  void push(double number, double weight);

  double mean() {   // weighted mean
    if (count_w>0.)
      return sum/count_w;
    else
      return 0.;
  }

  double variance() {
    if (count_w>0.)
      return (sum2-sum*sum/count_w)/count_w;
    else
      return 0.;
  }

  double meanERR() {  // relative error of the weighted mean
    if (count_w >0.)
      return 1./sqrt(count_w);
    else
      return 0.;
  }
  
  double count_eff() {
    if (count_w>0.)
      return count_w*count_w/count_w2;
    else
      return 0.;
  }
  
 protected:
  struct node
  {
    double entry;
    double weight;
    struct node *link;
  };
  node *last;  // pointer to last node
  double sum;  // weighted sum of all entries
  double sum2; // weighted sum of all entries^2
  unsigned long int count;  // # of entries
  double count_w;  // sum of weights
  double count_w2;  // sum of weights^2
};

#endif
