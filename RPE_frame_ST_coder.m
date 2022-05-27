function [LARc, CurrFrmResd ] = RPE_frame_ST_coder(s0)

%Offset compensation
alpha = 32735* 2^(-15);
s_of = zeros(160,1);
s_of(1) = s0(1);

for i = 2:160
    s_of(i) = s0(i) - s0(i-1) + alpha * s_of(i-1);

end 



%Premphasis
beta = 28180*2^(-15);
s = zeros(160,1);
s(1) = s_of(1);
for i = 2:160
    
    s(i) = s_of(i) - beta * s_of(i-1);
    
end 

%Autocorrelation


ACF = zeros(9,1);


for k = 1:9
    for i = k:160
        ACF(k) = ACF(k) + s(i) * s(i-(k-1));
    end
end



%calculating a_k (a = 1,2,...,8)

%first, calculating r and R
r = zeros(8,1);
for i = 2:9
    r(i-1) =  ACF(i);
end

R = zeros(8,8);

for i = 1:8
    for k = 1:8
        if i==k
            R(i,k) = ACF(1);
        end
        if i<k
            R(i,k) = ACF(k-i+1);
        end
        if i>k
            R(i,k) = ACF(i-k+1);
        end
    end
end    


%calculating w = a_k using Rw = r
    
w = R^(-1)*r;


x = zeros(9,1);
x(1) = 1;

for i=2:9
    x(i) = -w(i-1);
end

rc = poly2rc(x);


%calculating LAR
LAR = zeros(8,1);

for i=1:8
   if abs(rc(i)) <0.675
       LAR(i) = rc(i);
   elseif (abs(rc(i))>= 0.675) && (abs(rc(i))<0.950)
       LAR(i) = sign(rc(i)) * (2*abs(rc(i)) - 0.675);
   elseif (abs(rc(i))>= 0.950) && (abs(rc(i))<=1)
       LAR(i) = sign(rc(i)) * (8*abs(rc(i)) - 6.375);
   end
end


%LAR quantization
A = [20; 20; 20; 20; 13.637; 15; 8.334; 8.824];
B = [0; 0; 4; -5; 0.184; -3.5; -0.666; -2.235];

LARc = zeros(8,1);

for i = 1:8
    LARc(i) = A(i) * LAR(i) + B(i);
    LARc(i) = round( LARc(i) + sign(LARc(i)) * 0.5 );
end


%I've computed quantized LAR. 
%Now I have to "decode while encode" to produce d'(n)

%next there is nothing less than my code in the decoder file copied and 
%pasted
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
x1 = rc2poly(r_d);
x1 = x1';

w_d = -x1(2:9);

%now I have the so much desired ak' and I can move on ^^

s_p = zeros(160,1);
s_p(1) = s(1);

for i=2:160
    for k=1:8
        if i > k
            s_p(i) = s_p(i) + w_d(k) * s(i-k);
        end
    end    
end

%residual
CurrFrmResd = zeros(160,1);
for i=1:160
    CurrFrmResd(i) = s(i) - s_p(i);
end


end

