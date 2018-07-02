function [ q ] = Q_rand(seed)
%Q_RAND
%
% Q_rand() Computes a random unit quaternion
%
% Q_rand(seed) Computes a random unit quaternion with seed 
% T. Reynolds -- RAIN Lab
if( nargin < 1 )
    rng(randi(255,1))
else
    rng(seed)
end

q   = randn(4,1);
q   = q./norm(q);

end

