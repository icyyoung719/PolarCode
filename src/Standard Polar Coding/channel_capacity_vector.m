function C = channel_capacity_vector(N, alpha)
    % channel_capacity_vector calculates the channel capacity vector for a BEC.
    % N: Code length, must be a power of 2.
    % alpha: Erasure probability of the Binary Erasure Channel (BEC).
    % Returns:
    % C: Channel capacity vector of length N, which implies the proability
    % of error

    % Check if N is a power of 2
    if mod(log2(N), 1) ~= 0
        error('N must be a power of 2.');
    end

    % Initialize the channel capacity matrix
    max_level = log2(N);
    B = zeros(N, N);
    B(1, 1) = 1 - alpha; % Initialize the first channel capacity

    % Recursive calculation of channel capacities
    for level = 1:max_level
        N0 = 2^level; % Current code length
        for j = 1:(N0 / 2)
            % Update channel capacities for the current level
            B(N0, 2*j-1) = 2 * B(N0/2, j) - B(N0/2, j)^2; % Lower branch
            B(N0, 2*j) = B(N0/2, j)^2; % Upper branch
        end
    end

    % Extract the final channel capacity vector
    C = B(N, 1:N);
end