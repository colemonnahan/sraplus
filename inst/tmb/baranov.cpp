// baranov catches in TMB
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_SCALAR(m);

  DATA_VECTOR(sel);

  DATA_VECTOR(b_a);

  DATA_SCALAR(catches);

  PARAMETER(log_f);

  vector<Type> f_at_a = exp(log_f) * sel;

  vector<Type> c_at_a;

  c_at_a = (f_at_a / (f_at_a + m) * b_a * (Type(1.0) - exp(-(f_at_a + m))));

  Type catch_hat = sum(c_at_a);

  Type ss;

  ss = pow(catches- catch_hat,2);

  return ss;

}
