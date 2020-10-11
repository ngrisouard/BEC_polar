function z = h(r)

z = 0.7*besselj(0,r) + 0.5*besselj(1,r)+ 0.45*besselj(5,r)+ 0.87*besselj(20,r)+ 0.1*besselj(32,r)+ 0.05*besselj(40,r)
end