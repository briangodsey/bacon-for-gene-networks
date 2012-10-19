function [ ent ] = gaussEnt( prec )

%ent = log( (2*pi*exp(1)./prec).^0.5 );
ent = 0.5*( 1 + log(2*pi) -log(prec) );


end