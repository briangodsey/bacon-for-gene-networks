function [D] = scaleNorm( D, n )

if ( size(D,1)/n ~= round(size(D,1)/n) )
   disp('Warning: The number of rows is not divisble by the n parameter value'); 
end

if ( n == 0 )
    n = size(D,1);
end

for i = 1:floor( size(D,1)/n )
    Dsub = D( ((i-1)*n+1):(i*n), : );
    Dmeans = repmat( mean(Dsub,2), 1, size(D,2) );
    residuals = Dsub - Dmeans;
    noise = repmat( std(residuals,0,2), 1, size(D,2) );
    D( ((i-1)*n+1):(i*n), : ) = residuals ./ noise + Dmeans;
end

end