R1 = 1;
C = 0.25;
R2 = 2;
L = 0.2;
R3 = 10;
a = 100;
R4 = 0.1;
R0 = 1000;
G1 = 1/R1 + 1/R2;
G2 = -a*1/R3*(1/R4 + 1/R0);
G3 = 1/R4 + 1/R0;
G = [1/R1 -1/R1 0 0 0 1 0;
    -1/R1 G1 0 0 0 0 1;
    0 0 1/R3 0 0 0 -1;
    0 0 0 1 0 0 -a;
    0 0 0 -1/R4 G3 0 0;
    0 1 -1 0 0 0 0;
    1 0 0 0 0 0 0];
Cvals = [C -C 0 0 0 0 0;
        -C C 0 0 0 0 0;
        0 0 0 0 0 0 0;
        0 0 0 0 0 0 0;
        0 0 0 0 0 0 0;
        0 0 0 0 0 0 -L;
        0 0 0 0 0 0 0];
Vline = linspace(-10, 10, 100);
count = 1;
V0 = zeros(21,1);
V3 = zeros(21,1);
V1 = zeros(21,1);
for i=-10:10
    F = [0; 0; 0; 0; 0; 0; i];
    V = G\F;
    j = i +11;
    V0(j) = V(5);
    V3(j) = V(3);
    
end
figure(1)
plot([-10:10], real(V0))
plot([-10:10], real(V3))


for i=1:100
    F = [0; 0; 0; 0; 0; 0; i];
    V = (G+1i*i*C)\F;
    V3(i) = V(5);
end
figure(2)
plot([1:100],real(V3))

figure(3)
semilogx(w,10*log10(real(V3)/1))

for i = 1:10000
    w = pi+randn(1)*0.05;
    F = [0; 0; 0; 0; 0; 0; i];
    V = (G+1i*w*C)\F;
    V3(i) = V(5);
end
A = 10*log10(real(V3)/1);
figure(4)
%hist(A, 20)


%Part 2
dt = 0.001;
A = C/dt + G;
endtime = 1000*dt;
time = 0;
Vold = zeros(1,1000);
count = 1;
V0B = zeros(1000,1);
V3B = zeros(1000,1);
V1B = zeros(1000,1);
for i = dt:dt:endtime
    if count > 30
        F = [0; 0; 0; 0; 0; 0; 0];
        V = G\F;
        V0B(count) = i;
        V3B(count) = V(3);
    else
        F = [0; 0; 0; 0; 0; 0; 1];
        V = A\(F + C*Vold./dt);
        Vold = V;
        V1B(count) = i;
        V0B(count) = V(5);
        V3B(count) = V(3);
    end
    count = count + 1;
end

figure(5)
plot(V1B, V0B)
plot(V1B, V3B)

V0F = fft(V0B);
V3F = fft(V3B);
V0Fshift = fftshift(V0F);
V3Fshift = fftshift(V3F);

figure(6)
semilogy(abs(V3Fshift))
hold on
semilogy(abs(V0Fshift));

count = 1;
for i = dt:dt:endtime
    s = sin(2*pi*(1/0.03)*i);
    F = [0; 0; 0; 0; 0; 0; s];
    V = A\(F + C*Vold./dt);
    Vold = V;
    V1B(count) = i;
    V0B(count) = V(5);
    V3B(count) = s;
    count = count +1;
end

figure(7)
plot(V1B, V0B)
hold on
plot(V1B, V3B)

V0F = fft(V0B);
V3F = fft(V3B);
V0Fshift = fftshift(V0F);
V3Fshift = fftshift(V3F);

figure(8)
semilogy(abs(V0Fshift));
hold on
semilogy(abs(V3Fshift));

count = 1;
for i = dt:dt:endtime
    s = exp(-(i-0.06).^2/(2*dt));
    V = A\(F + C*Vold./dt);
    Vold = V;
    V1B(count) = i;
    V0B(count) = V(5);
    V3B(count) = s;
    count = count +1;
end

figure(9)
plot(V1B, V0B)
hold on
plot(V1B, V3B)

V0F = fft(V0B);
V3F = fft(V3B);
V0Fshift = fftshift(V0F);
V3Fshift = fftshift(V3F);

figure(10)
semilogy(abs(V0Fshift));
hold on
semilogy(abs(V3Fshift));

%part 3
Cn = 0.0001;
  C = [C -C 0 0 0 0 0;
    -C C 0 0 0 0 0;
     0 0 Cn 0 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 -L;
     0 0 0 0 0 0 0];
 count = 1;
 for i = dt:dt:endtime
     s = exp(-(i-0.06).^2/(2*dt));
     ran = randn(1,1)*0.01;
     F = [0 0 ran 0 0 0 s];
     V = A\(F + C*Vold./dt);
     Vold = V;
     V1B(count) = i;
     V0B(count) = V(5);
     V3B(count) = s;
     count = count + 1;
 end
 
figure(11)
plot(V1B,V0B)
hold on
plot(V1B,V3B)

V0F = fft(V0B);
V3F = fft(V3B);
V0Fshift = fftshift(V0F);
V3Fshift = fftshift(V3F);

figure(12)
semilogy(abs(V0Fshift));
hold on
semilogy(abs(V3Fshift));

    




    
    