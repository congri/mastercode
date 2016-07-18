tic
for n= 1:10000
    for i = 1:10000
        
        p(i) = rand;
        
    end
end
t1 = toc

A = 10*rand(1000) -5;

tic
for n = 1:50
    B = inv(A);
end
t2 = toc