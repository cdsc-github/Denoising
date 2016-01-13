% Get a clean input image
uexact = imread('knee-mri.png');
uexact = mean(double(uexact)/255,3);
% Simulate Rician noise
sigma = 0.05;
f = ricianrnd(uexact,sigma);

% Print f
fid = fopen('noisy_array', 'w');
dim = ndims(f);
m = size(f, 1);
n = size(f, 2);
fprintf('dim: %d   m: %d   n: %d\n', dim, m, n);
fprintf(fid, '%d %d\n', m, n);

for i = 1:m
	for j = 1:n
		fprintf( fid, '%f ', f(i,j) );
	end
	fprintf( fid, '\n');
end

