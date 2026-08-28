% Define EQN
A=[3 1; 2 2];
B = [5; 1];

%solve non-symbolic
solution= A\B;

%defines variable 
x = solution(1);
y = solution(2);

%shows results
fprintf('x = %g\n', x);
fprintf('y = %g\n', y);