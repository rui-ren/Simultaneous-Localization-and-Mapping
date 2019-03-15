function I = simpsons(f,a,b,n)

%  Written by Juan Camilo Medina  - The University of Notre Dame
%  09/2010 (copyright Dr. Simpson)
% 


if numel(f)>1 % If the input provided is a vector
    n=numel(f)-1; h=(b-a)/(n);
    I= h/3*(f(1)+2*sum(f(3:2:end-2))+4*sum(f(2:2:end))+f(end));
else % If the input provided is an anonymous function
    h=(b-a)/(n); xi=a:h:b;
    I= h/3*(f(xi(1))+2*sum(f(xi(3:2:end-2)))+4*sum(f(xi(2:2:end)))+f(xi(end)));
end

