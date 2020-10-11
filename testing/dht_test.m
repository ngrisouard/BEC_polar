z = linspace(1,100)
%testing = test_f(z)
%testing_f = 0.7*besselj(0) + 0.5*besselj(1)+ 0.45*besselj(5)+ 0.87*besselj(20)+ 0.1*besselj(32)+ 0.05*besselj(40)



trans = dht(@h, 1)
plot(trans)
