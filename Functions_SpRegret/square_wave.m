function y = square_wave(t, frequency, dutyCycle, phaseShift)
% SQUARE_WAVE Generates a square wave signal with an optional phase shift
%
% y = SQUARE_WAVE(t, frequency, dutyCycle, phaseShift) returns a square wave signal y
% with the given frequency, duty cycle, and phase shift over the time array t.
%
% Inputs:
%   t - Time array
%   frequency - Frequency of the square wave (in Hz)
%   dutyCycle - Duty cycle of the square wave (percentage)
%   phaseShift - Phase shift of the square wave (in radians)
%
% Output:
%   y - Square wave signal

% Default duty cycle to 50% if not provided
if nargin < 3
    dutyCycle = 50;
end

% Default phase shift to 0 if not provided
if nargin < 4
    phaseShift = 0;
end

% Calculate the period of the wave
T = 1 / frequency;

% Calculate the high time based on the duty cycle
highTime = dutyCycle / 100 * T;

% Adjust the time array for phase shift
t_shifted = t + phaseShift / (2 * pi * frequency);

% Initialize the output array
y = zeros(size(t));

% Generate the square wave
for i = 1:length(t)
    % Compute the time within the period, considering phase shift
    timeInPeriod = mod(t_shifted(i), T);

    % Determine if we are in the "high" or "low" part of the cycle
    if timeInPeriod < highTime
        y(i) = 1;
    else
        y(i) = -1;
    end
end
end
