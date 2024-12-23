
% Example usage:
N = 8; % Codeword length
K = 4; % Number of information bits

info_bits = [1 0 1 0]; % Example input bits
[u, x] = polar_encode(N, K, info_bits);
disp('Encoded vector u:'), disp(u);
disp('Encoded codeword x:'), disp(x);