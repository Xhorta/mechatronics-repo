% Create an array to store the 10 numbers
random_numbers = zeros(1, 10);
% Generate random numbers until it has generated 10 random values
for i = 1:10
    random_numbers(i) = rand();
end

% Creates plot
plot(1:10, random_numbers, 'green--o', LineWidth=2);

% axis's labels
xlabel('Number');
ylabel('Random Number');

% creates title
title('Xaviers 10 Random Numbers');
% Turn on the grid
grid on;