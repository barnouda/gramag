function g=f(x,y)
sigma = 0.5;
somme = x.*x+y.*y;
g=exp(-sigma.*somme);
end