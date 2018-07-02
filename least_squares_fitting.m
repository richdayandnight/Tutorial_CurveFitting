# Curve Fitting
# Rich Yap

% Least squares fitting


# Given z (time in years) and T (population)
clear all;
z = [0,10,20,30,40];
T = [35.80, 47.40, 61.95, 77.99, 93.73];

# Initialize matrix A
A = zeros(4,4);

# Initialize z vectors
z_squared = [];
z_cube = [];
z_fourth = [];
z_fifth = [];
z_sixth = [];

# Fill z vectors with the square, cube, ..., sixth power of z
for i = 1:length(z) 
  z_squared(i) = z(i)**2;
  z_cube(i) = z(i)**3;
  z_fourth(i) = z_squared(i)**2;
  z_fifth(i) = z_squared(i) * z_cube(i);
  z_sixth(i) = z_cube(i)**2;
endfor

# Get the sum of the filled vectors above and put it in the z_vector
z_vector = [];
z_vector(end+1) = sum(z);
z_vector(end+1) = sum(z_squared);
z_vector(end+1) = sum(z_cube); 
z_vector(end+1) = sum(z_fourth);
z_vector(end+1) = sum(z_fifth);
z_vector(end+1) = sum(z_sixth);


# Build the A matrix using the z_vector (Note: Ax = b)
A(1,1) = length(z);
A(1,2:end) = z_vector(1:3);
A(2,1:end) = z_vector(1:4);
A(3,1:end) = z_vector(2:5);
A(4,1:end) = z_vector(3:6);

# Build the b matrix 
b = [];
b(end+1) = sum(T);
b(end+1) = sum(z(1:end) .* T(1:end));
b(end+1) = sum(z_squared(1:end) .* T(1:end));
b(end+1) = sum(z_cube(1:end) .* T(1:end));
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
linx = linspace(min(z), max(z));

# Evauate the polynomial with the created list of x coordinates
liny = evalpoly(poly, linx(:));

# Plot the Temperature vs Water Depth graph
plot(linx(:), liny(:), "b-");

z = [0,5,10,15,20,25,30,35,40];
T = [35.80, 41.29, 47.40, 54.32, 61.95, 69.83, 77.99, 86.27, 93.73];

# Plot the given data points
plot(z(:), T(:), "or");
title("Temperature vs Water Depth graph");
ylabel("Temperature");
xlabel("Water Depth");
legend("T(z)", "Data point");
hold off;

subplot(2,1,2)
hold on;

# Evauate the first_derivative with the created list of x coordinates
liny_first_d = evalpoly(first_derivative, linx(:));

# Plot the Temperature vs Water Depth graph
plot(linx(:), liny_first_d(:), "b-");
legend("Gradient = T(z)'");
title("Gradient");
hold off;


