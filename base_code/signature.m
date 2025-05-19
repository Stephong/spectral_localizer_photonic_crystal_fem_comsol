function sig = signature(X)

% Uses  Sylvester's law of inertia
% and the sparse LDL^adjoint factorization
% that comes with Matlab.
% It is assumed that X  is (nearly)
% self-adjoint, and with a reasonable 
% gap in the spectrum at 0.

X = (1/2)*(X + X');

n = length(X);
X2 = [real(X), imag(X); -imag(X), real(X)];
[L,D,p]=ldl(X2,'vector');
sig = 0;
j=1;
while j <= 2*n
	if j<2*n && D(j,j+1) ~= 0
		sig = sig + sum( eig( [D(j,j),D(j,j+1);D(j+1,j),D(j+1,j+1)] ) > 0 );
		j = j+2;
	else
		sig = sig + (D(j,j) > 0);
		j = j+1;
	end
end
sig = sig - n;


%sig = sum(eigs(D,n+4,'LA')>0) - n;

end
