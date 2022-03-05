%% This code show the matrix (for example A in AX=B)

function f = matrixillustrator(l,d,u)
   
    n = length(d); % n is the number of rows
    f = zeros(n,n); % Specify the size of the matrix
    f(1,1) = d(1); % Specify some matrix elements
    f(1,2) = u(1); % Specify some matrix elements
    f(n,n-1) = l(n); % Specify some matrix elements
    f(n,n) = d(n); % Specify some matrix elements
    
    for i = 2:n-1
        f(i,i-1) = l(i);
        f(i,i) = d(i);
        f(i,i+1) = u(i);
    end
    