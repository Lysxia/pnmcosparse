d = 15;

w1 = 1;
w2 = 1;

p = w1*(d-1)+w2*d

l = 20;

Z = rand(d,1);

Z

[X,L1,L2] = csa_projection_FUS(Z,w1,w2,l)

plot(1:d,X,1:d,Z);