function [ ent ] = gamEnt( shape, scale )

ent = shape - log(scale) + gammaln(shape) + (1-shape).*psi(shape);

end