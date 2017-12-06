%%							cgist_demo.m
%   This script demonstrates how to use the function cgist.m to solve L1
%   regularized problems.  The script will generate a sparse signal with
%   the specified dynamic range and number of non-zero entries.  The script 
%   then performs signal recovery using the following measurment operators
%   and regularizers:
%		Experiment 1:  The signal is measured using a random Gaussian
%		matrix.  Recovery is done using L1 regularized least squares.  The
%		measurement operator is handed to cgist.m as a dense matlab matrix.
%
%		Experiment 2:  Measurements are taken in a discrete cosine
%		transform (DCT) basis.  Recovery is done using a nonnegativity
%		constraint (i.e. nonnegative least squares regression).  The
%		measurement operator is handed to cgist.m as a function handle.
%
%		Experiment 3:  Measurements are taken in a Fourier (DFT/FFT) basis.
%		Recovery is done using L1 regularized least squares.  Note that
%		this experiment involves complex valued data.
%
%		Experiment 4:  Measurements are taken using a Fourier basis, but
%		this time recovery is performed by minimzing L2 subject to
%		a convex constraint of the form |x| < lambda, where |-|
%		denotes the 1-norm. This experiment also involves complex valued
%		data.
%
%		Experiment 5:  Measurements are taken using a Fourier basis, but
%		this time recovery is performed by minimzing L1 subject to
%		a convex constraint of the form ||Ax -f || < epsilon, where ||-||
%		denotes the 2-norm. This experiment also involves complex valued
%		data.
%
%	Required Files:  cgist.m
%
% Author: 
%		Thomas Goldstein
%		Department of Electrical Engineering
%		Stanford University
%		http://www.stanford.edu/~tagoldst	



%%  Parameters for the tests
N = 500;            % Length of sparse vector
measurements = 100;  % Number of measurements taken via inner pruducts
sparsity = 10;       % Sparsity of the test vector we want to recover 
sigma = 0.01;       % Standard deviation of noise to be added to measurements
decibels = 20;		% The dynamic range of signal, measured in decibels	

disp('Initializing tests...');
%%  construct a sparse signal - the thing we want to recover
signal = zeros(N,1);
perm = randperm(N); %  randomly choose the support of the signal
if decibels>0		% Define the exponents needed to build a signal with the desired dynamic range
	log_intensity = 0:decibels/10/(sparsity-1):(decibels/10);
else
	log_intensity = zeros(sparsity,1);
end
signal(perm(1:sparsity)) = 10.^log_intensity;	% The signal height ranges from 1-10^(decibels/10)



%%   L1 compressed sensing using Gaussian matrix
disp('Begin l1 Gaussian matrix test...')

A = randn(measurements ,N);  % the measurement matrix is random iid Gaussian
f = A*signal;               % form measurements by multiplying by A
f = f + sigma*randn(size(f)); % add noise
    %  Perform signal recovery from noisy data, using non-negative L1
	%  Note: Since "A" is a dense matrix, we don't need to hand in its
	%  transpose
	
mu = 100;	%  Fidelity constant
	
[sol, count, subgrad, out] = cgist(A,[],f,mu,'l1');

disp('    Test Completed:');
disp(['    Matrix Multiplications = ',num2str(count)]);
disp(['    Subgradient norm = ',num2str(subgrad)]);
disp(['    Atoms detected = ',num2str(sum(sol~=0))]);

%close all; plot(rec);


%%  L1 Compressed sensing using DCT Matrix
disp('Begin nonnegative DCT matrix test...')

%  build the row selector matrix
R = zeros(N,1); 
perm = randperm(N);
R(perm(1:measurements)) = 1;
R(1) = 1;  % make sure we measure the lowest order mode.  Otherwise CS doesn't work as well.

% build the function handles for matrix multiplication  
A  = @(x) sqrt(N)*R.*dct(x);
At  = @(x) sqrt(N)*idct(R.*x);

f = A(signal);      %  compute measurements, and add noise
f = R.*(f + sigma*randn(size(f)));

mu = 100;	%  choose mu>0 to use both L1 and nonnegativity

%  Perform signal recovery from noisy data
%  Note:  We must now hand in a function handle to A, as well as to At.
[sol, count, subgrad, out] = cgist(A,At,f,mu,'nonnegative');

disp('    Test Completed:');
disp(['    Matrix Multiplications = ',num2str(count)]);
disp(['    Subgradient norm = ',num2str(subgrad)]);
disp(['    Atoms detected = ',num2str(sum(sol~=0))]);

% close all; plot(rec);


%%  L1 Compressed sensing using complex-valued Fourier Matrix
disp('Begin l1 Complex Fourier matrix test...')

%  build the row selector matrix
R = zeros(N,1); 
perm = randperm(N);
R(perm(1:measurements)) = 1;
R(1) = 1;  % make sure we measure the lowest order mode.  Otherwise CS doesn't work as well.

% build the function handles for matrix multiplication  
A  = @(x) R.*fft(x)/sqrt(N);
At  = @(x) sqrt(N)*ifft(R.*x);

f = A(signal);      %  compute measurements, and add noise
f = R.*(f + sigma*randn(size(f)));

mu = .1;

%  Perform signal recovery from noisy data
%  Note:  We must now hand in a function handle to A, as well as to At.
[sol, count, subgrad, out] = cgist(A,At,f,mu,'l1');

disp('    Test Completed:');
disp(['    Matrix Multiplications = ',num2str(count)]);
disp(['    Subgradient norm = ',num2str(subgrad)]);
disp(['    Atoms detected = ',num2str(sum(sol~=0))]);

% close all; plot(rec);



%%  Sparse Signal Recovery Using A Constraint of the form |x|<lambda
disp('Begin L1 Constrained FFT Test...')

%  We'll re-use the data from the last experiment, non need to rebuild it.

%  Define epsilon from the data noise level
lambda = 0.9*sum(abs(signal));

%  Perform signal recovery from noisy data
%  Note:  This time, we minimizing L2 subject to a constraint
[sol, count, subgrad] = cgist(A,At,f,lambda,'l1_constrained');

disp('    Test Completed:');
disp(['    Matrix Multiplications = ',num2str(count)]);
disp(['    Subgradient norm = ',num2str(subgrad)]);
disp(['    |x|  = ',num2str(sum(abs(sol)))]);
disp(['    lambda  = ',num2str(lambda)]);
disp(['    Atoms detected = ',num2str(sum(sol~=0))]);


%%  Sparse Signal Recovery Using A Constraint of the form ||Ax-f||<epsilon
disp('Begin L2 Constrained FFT Test...')

%  We'll re-use the data from the last experiment, non need to rebuild it.

%  Define epsilon from the data noise level
epsilon = 10*sigma*sqrt(measurements);

% This time, well use some fancy options
opts = [];
opts.record_objective = true;	% record the value of the objective function

%  Perform signal recovery from noisy data
%  Note:  This time, we minimizing L1 subject to a constraint
[sol, count, subgrad, out] = cgist(A,At,f,epsilon,'l2_constrained',opts);

disp('    Test Completed:');
disp(['    Matrix Multiplications = ',num2str(count)]);
disp(['    Subgradient norm = ',num2str(subgrad)]);
disp(['    || A*x-f  ||  = ',num2str(norm(A(sol)-f))]);
disp(['    epsilon  = ',num2str(epsilon)]);
disp(['    Atoms detected = ',num2str(sum(sol~=0))]);



%%  Make a plot
close all; 
figure;
% First, plot signal
subplot(1,2,1);
plot(abs(sol));
title('Recovered Signal');

% Next, plot objective function on last inner solve step
subplot(1,2,2);
%  Find out how many iterations the last step took
iter = out.iterates_per_continuation_step(end);	
%  plot the objective for the iterations on the last step
plot(1:iter, out.objectives( (end-iter+1):end) );
title('Objective Function');

