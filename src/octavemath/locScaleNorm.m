function [D] = locScaleNorm( D, n )

D = scaleNorm( locationNorm( D, n ), n );

end