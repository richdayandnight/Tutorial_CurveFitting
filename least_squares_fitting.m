# Curve Fitting
# Rich Yap

% Least squares fitting


# Given year (time in years) and population (million)
clear all;
year = [0,10,20,30,40];
population = [35.80, 47.40, 61.95, 77.99, 93.73];

# Initialiyeare matrix A
A = zeros(4,4);

# Initialiyeare year vectors
year_squared = [];
year_cube = [];
year_fourth = [];
year_fifth = [];
year_sixth = [];

# Fill year vectors with the square, cube, ..., sixth power of year
for i = 1:length(year) 
  year_squared(i) = year(i)**2;
  year_cube(i) = year(i)**3;
  year_fourth(i) = year_squared(i)**2;
  year_fifth(i) = year_squared(i) * year_cube(i);
  year_sixth(i) = year_cube(i)**2;
endfor

# Get the sum of the filled vectors above and put it in the year_vector
year_vector = [];
year_vector(end+1) = sum(year);
year_vector(end+1) = sum(year_squared);
year_vector(end+1) = sum(year_cube); 
year_vector(end+1) = sum(year_fourth);
year_vector(end+1) = sum(year_fifth);
year_vector(end+1) = sum(year_sixth);


# Build the A matrix using the year_vector (Note: Ax = b)
A(1,1) = length(year);
A(1,2:end) = year_vector(1:3);
A(2,1:end) = year_vector(1:4);
A(3,1:end) = year_vector(2:5);
A(4,1:end) = year_vector(3:6);

# Build the b matrix 
b = [];
b(end+1) = sum(population);
b(end+1) = sum(year(1:end) .* population(1:end));
b(end+1) = sum(year_squared(1:end) .* population(1:end));
b(end+1) = sum(year_cube(1:end) .* population(1:end));
b = b.';


# Use Cholesky Decomposition to get a0, a1, a2, a3
# Get L and D in the LDL'   
# reference: mathworks.com/matlabcentral/fileexchange/47-ldlt?focused=5033716&tab=function
n = size(A,1);
L = zeros(n,n);
for j=1:n,
  if (j > 1),
    v(1:j-1) = L(j,1:j-1).*d(1:j-1);
    v(j) = A(j,j)-L(j,1:j-1)*v(1:j-1)';
    d(j) = v(j);
    if (j < n),
      L(j+1:n,j) = (A(j+1:n,j)-L(j+1:n,1:j-1)*v(1:j-1)')/v(j);
    endif
  else
    v(1) = A(1,1);
    d(1) = v(1);
    L(2:n,1) = A(2:n,1) / v(1);    
  endif
endfor

D=diag(d);
L=L+eye(n);

# Ly = b
# Solve for y
y = L\(b);

# DL.'x = y
# Solve for x
x = (D*L.')\(y);

printf("\nCurve Fitting\n1.) Matrices formed in least-squares polynomial curve fitting\nA matrix:\n");
A
printf("\nb matrix:\n");
b

# Flip the x polynomial produced to fit the polynomial format
poly = flip(x);
printf("\n======================================================================\n");
printf("\n2.)The following cubic polynomial is interpolates the data points: \n")
polyout(poly, 'x');
printf("where:\n    a_0 = %f\n     a_1 = %f\n     a_2 = %f\n     a_3 = %f\n", poly(1), poly(2), poly(3), poly(4));
printf("SOLUTION: Cholesky LDLt\n");
L
D
L_transpose = L.'
printf("Ly = b // below is the solved y\n");
y
printf("D(Lt)x = y  // below is the solved x\n");
x


# Get first derivative of polynomial generated
first_derivative = polyder(poly.');

# Get second derivative of the polynomial generated
second_derivative = polyder(first_derivative);

printf("\n======================================================================\n");
printf("\n2.) The following first derivative of the polynomial generated:\n");
polyout(first_derivative, 'x');

printf("\n======================================================================\n");
printf("\n3.) See the second subplot generated in figure 2 for the gradient of the first graph\n");
subplot(2,1,1)
hold on;

# Create list of x coordinates
linx = linspace(min(year), max(year));

# Evauate the polynomial with the created list of x coordinates
liny = evalpoly(poly, linx(:));

# Plot the Year vs Population graph
plot(linx(:), liny(:), "b-");

year = [0,5,10,15,20,25,30,35,40];
population = [35.80, 41.29, 47.40, 54.32, 61.95, 69.83, 77.99, 86.27, 93.73];

# Plot the given data points
plot(year(:), population(:), "or");
title("Year vs Population graph");
ylabel("Year");
xlabel("Population");
legend("population(year)", "Data point");
hold off;

subplot(2,1,2)
hold on;

# Evauate the first_derivative with the created list of x coordinates
liny_first_d = evalpoly(first_derivative, linx(:));

# Plot the Year vs Population graph
plot(linx(:), liny_first_d(:), "b-");
legend("Gradient = population(year)'");
title("Gradient");
hold off;


