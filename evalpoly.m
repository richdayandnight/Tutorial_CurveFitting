function [y] = evalpoly(poly,x)
  y = [];
  for i=1:length(x)
    y(i) = horner(poly,x(i));
  endfor
end