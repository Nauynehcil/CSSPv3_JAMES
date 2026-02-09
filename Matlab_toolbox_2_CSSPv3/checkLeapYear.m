function [isLeapYear, daysInMonth] = checkLeapYear(year)
% CHECKLEAPYEAR Determines leap years and returns days in each month
%   [isLeapYear, daysInMonth] = checkLeapYear(year)
%
%   Input:
%       year        - Integer year, can be a scalar or a vector (e.g., 2020 or [2000, 2024])
%
%   Output:
%       isLeapYear  - Logical array of the same size as "year", true indicates a leap year
%       daysInMonth - Matrix [length(year) x 12] containing the number of days per month

% Input validation
if ~isnumeric(year) || any(mod(year, 1) ~= 0)
    error('checkLeapYear:InputError', 'Input must be an integer year or an array of years.');
end

% Reshape to column vector for consistent matrix operations
year = year(:);
nYear = length(year);

% Determine leap years
% Rule: (Divisible by 4 AND not divisible by 100) OR (Divisible by 400)
isLeapYear = (mod(year, 4) == 0 & mod(year, 100) ~= 0) | mod(year, 400) == 0;

% Initialize daysInMonth matrix (each row represents a year, each column a month)
% Default February to 28 days
daysInMonth = repmat([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], nYear, 1);

% Update February for identified leap years to 29 days
daysInMonth(isLeapYear, 2) = 29;

end