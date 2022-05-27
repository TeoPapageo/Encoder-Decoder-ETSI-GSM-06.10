function [s0] = RPE_frame_ST_decoder(LARc, CurrFrmSTResd)
A = [20; 20; 20; 20; 13.637; 15; 8.334; 8.824];
B = [0; 0; 4; -5; 0.184; -3.5; -0.666; -2.235];

%calculating LAR from quantized LARc
LARd = zeros(8,1);

for i = 1:8
   LARd(i) = (LARc(i) - B(i)) / A(i);
end


%calculating reflection coefficients
r_d = zeros(8,1);

for i=1:8
    if abs(LARd(i)) < 0.675
        r_d(i) = LARd(i);
    elseif (abs(LARd(i))>= 0.675) && (abs(LARd(i))<1.225)
        r_d(i) = sign(LARd(i)) * (0.5 * abs(LARd(i)) + 0.3375);
    elseif (abs(LARd(i))>= 1.225) && (abs(LARd(i))<=1.625)
        r_d(i) = sign(LARd(i)) * (0.125 * abs(LARd(i)) + 0.796875);
    end
end



%calculating a_k dequantized 

x = rc2poly(r_d);
x = x';
w_d = -x(2:9);



%using the short time synthesis filter to get the decoded version of s(n)

f1 = [1 -w_d(1) -w_d(2) -w_d(3) -w_d(4) -w_d(5) -w_d(6) -w_d(7) -w_d(8)];
f2 = 1;

s_d = filter(f2, f1, CurrFrmSTResd);


%postprocessing
beta = 28180*2^(-15);
s0 = zeros(160,1);
s0(1) = s_d(1);


for i=2:160
    s0(i) = s_d(i) + beta*s0(i-1);
end


end

