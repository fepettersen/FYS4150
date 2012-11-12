infile = textread('results.txt');
E = infile(1,:);
M = infile(2,:);
N = infile(3,:);

plot(E,N)
figure()
plot(M,N)