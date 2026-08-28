num1 = input('first number: ');
num2 = input('second number: ');

%display result 
result= numOp(num1,num2);
fprintf('Result = %g\n', result);

% creates the function 
function result = numOp(num1, num2)

    % If both nums odd, sum numbers 
    % odd ==> ~=0
    if mod(num1, 2) ~= 0 && mod(num2, 2) ~= 0
        result = num1 + num2;

    % If both nums even, LargerNumber - SmallerNumber
    % even ==> ==0 
    % takes the smaller number and subtracts it from the larger value
    elseif mod(num1, 2) == 0 && mod(num2, 2) == 0
        result = max(num1, num2) - min(num1, num2);

    % If not both odd && not both even
    % one is even one is odd 
    %multiply the numbers together 
    else
        result = num1 * num2;
    end
end