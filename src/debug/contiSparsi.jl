
using Symbolics

@variables x1 x2 y1 y2 h dist
d = ((x1-x2)^2 + (y1-y2)^2)^(1/2)
d2 = (x1-x2)^2 + (y1-y2)^2
dx = x1-x2
dy = y1-y2
dxy = [dx,  dy]
# K = 1/h^2 * 1/(d/h+1)^2
W = (1/h^2) * 1/((d/h)+1)^2 
W = (1/h^2) * 1/((d/h)^2+1)

Dx1 = Differential(x1)
Dx1K = simplify(expand_derivatives( Dx1(W) ))
Dx1K = simplify(substitute(Dx1K, Dict([x1-x2 => dx, y1-y2 => dy, ((x1-x2)^2 + (y1-y2)^2)^0.5 => d, (x1-x2)^2 + (y1-y2)^2 => d2 ])))
Dx1K = (-2(x1 - x2)) / ((h + ((x1 - x2)^2 + (y1 - y2)^2)^0.5)^3*(((x1 - x2)^2 + (y1 - y2)^2)^0.5))
Dx1K = -2*dx / ((h + d)^3*d)


(-2(x1 - x2)) / ((h^2 + x1^2 - 2x1*x2 + x2^2 + y1^2 - 2y1*y2 + y2^2)^2)
(-2dx) / ((h^2 + d2)^2)
