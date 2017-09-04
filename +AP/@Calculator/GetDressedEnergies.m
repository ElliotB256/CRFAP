function [ eigE, eigV ] = GetDressedEnergies( context, omega0, theta, gamma )
%GETDRESSEDENERGIES Calculates the given dressed eigenenergies.
%   Calculates the dressed eigenenergies for the given RF dressing fields
%   and atom system described by context.
%
%   The energies are calculated at coordinates omega0, theta, gamma. omega0
%   corresponds to the Zeeman splitting (in MHz), while theta and gamma are
%   optional coordinates dictating the angle of the quantisation axis.
%
%   Syntax:
%    [ eigE, eigV ] = GetDressedEnergies( context, omega0, theta, gamma )

% Verify inputs. Check that theta and phi are defined.
if nargin < 3
    theta = zeros(size(omega0));
end

if nargin < 4
    gamma = zeros(size(omega0));
end

% If theta or phi are rank one, resize them to same size as omega0.
if length(theta) == 1
    theta = repmat(theta, length(omega0), 1);
end

if length(gamma) == 1
    gamma = repmat(gamma, length(omega0), 1);
end

% Change all vectors to columns.
gamma = gamma(:); theta = theta(:); omega0 = omega0(:);

% Ensure vectors are all the same length.
if length(gamma) ~= length(theta) || length(omega0) ~= length(theta)
    error('Incorrect size for theta and/or phi. These should either be undefined, length=1 or length=length(omega0).');
end

% Input verification complete.


hs = context.HilbertSpaceSize();

% Initialise arrays for eigE and eigV
eigE = zeros(hs, length(omega0));
eigV = zeros(hs, hs, length(omega0));

if context.Parallel
    parfor i=1:length(omega0)
        [ eigE(:,i), eigV(:,:,i) ] = compute(context, omega0(i), theta(i), gamma(i));
    end
else
    for i=1:length(omega0)
        [ eigE(:,i), eigV(:,:,i) ] = compute(context, omega0(i), theta(i), gamma(i));
    end
end

end

