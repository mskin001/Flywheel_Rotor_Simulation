function [s] = power_law_epoxy(mstiff)
  
  global t
  
  jf0 = 0.0003841;
  A = 1.181e-5;
  n = 0.2049;
  
  j = jf0 + A*t^n;
  
  s(1) = mstiff(1);
  s(2) = j;
  s(3) = s(2);
  
