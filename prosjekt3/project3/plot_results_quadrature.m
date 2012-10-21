sum_laguerre =[0.186457,0.189759,0.191081,0.19174,0.192113,0.192343,0.192493,0.192595,0.192667]; 
sum_lagendre = [0.0719797,0.239088,0.156139,0.195817,0.177283,0.189923,0.184417,0.189586,0.18756];
length(sum_laguerre);
n = linspace(10,50,9);
result = ones(length(n))*(5*pi*pi/(16*16));

plot(n,sum_lagendre,'b-o')
hold on
plot(n,sum_laguerre,'r-o')
plot(n,result,'g')
legend('lagendre','laguerre','0.192765')
xlabel('number of integration points')
hold off