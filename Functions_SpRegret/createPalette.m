function palette = createPalette(colorA, colorB, N)
    % colorA and colorB are [R, G, B] vectors for the start and end colors
    % N is the number of colors to generate

    % Ensure inputs are valid
    if N < 2
        error('N must be at least 2 to create a palette.');
    end
    
    % Interpolate between colorA and colorB
    palette = zeros(N, 3); % Preallocate for speed
    for i = 1:3 % Interpolate each RGB channel
        palette(:, i) = linspace(colorA(i), colorB(i), N);
    end
end
