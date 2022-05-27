function [LARc,Nc,bc,CurrFrmExFull,CurrFrmSTResd] = RPE_frame_SLT_coder(s0, PrevFrmSTResd)

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
        ACF(k) =  ACF(k) + s(i) * s(i-(k-1));
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

%residual d(n)
PreCurrFrmSTResd = zeros(160,1);
for i=1:160
    PreCurrFrmSTResd(i) = s(i) - s_p(i);
end
CurrFrmSTResd = zeros(160,1);
%xorizo ta 160 torina residuals kai ta residuals tou proigoumenou frame se 40ades
resids = buffer(PreCurrFrmSTResd, 40); %ta d xorismena se 40ades // epipedo subframes
prevresids = buffer(flip(PrevFrmSTResd), 40); %ta d' tou proigoumenou frame, anapoda kai xorismena se subframes
flippedprevresids = flip(PrevFrmSTResd); %ta d' tou proigoumenou frame, apla anapoda
flippedcurrResids = buffer(flip(PreCurrFrmSTResd), 40); %ta d' tou parontos subframe, anapoda kai xorismena se 40ades
Nc = zeros(4,1); %ousiastika to N den kvantizetai
b = zeros(4,1);
bc = zeros(4,1); %quantized b
bd = zeros(4,1); %dequantized b
predFrmRstd = zeros(160,1); %d^n
CurrFrmExFull = zeros(160,1); %e(n)

%sti sinexeia epanalmvano tin idia diadikasia alla diaforetika gia kathe
%subframe

            
            a = reshape([prevresids(:,1) prevresids(:,2) prevresids(:,3)],[],1);

            
            [Nc(1),b(1)] = RPE_subframe_LTE(resids(:,1), a);
           
            
            
            
               %kvantismos twn b
               if b(1) <= 0.2
                    bc(1) = 0; 
               end
               if b(1) > 0.2 && b(1) <= 0.5
                    bc(i) = 1; 
               end
               if b(1) > 0.5 && b(1) <= 0.8
                    bc(1) = 2; 
               end
               if b(1) > 0.8
                    bc(1) = 3; 
               end
               
            %d^n
            for k=1:40
                predFrmRstd(k) = bc(1)*flippedprevresids(Nc(1)+1-k); 
            end
            
            %e = d-d^
            for k = 1:40
                CurrFrmExFull(k) = PreCurrFrmSTResd(k) - predFrmRstd(k);
            end
            
            %dequantizing bc
               if bc(1) == 0
                    bd(1) = 0.1; 
               end
               if bc(1) == 1
                    bd(1) = 0.35;  
               end
               if bc(1) == 2
                    bc(1) = 0.65; 
               end
               if bc(1) == 3
                    bd(1) = 1; 
               end
            
            %d'
            for k=1:40
                CurrFrmSTResd(k) = CurrFrmExFull(k) + bd(1)*flippedprevresids(Nc(1)+1-k);
            end
            
    
            a = reshape([flippedcurrResids(:,4); prevresids(:,1); prevresids(:,2)],[],1);
            [Nc(2),b(2)] = RPE_subframe_LTE(resids(:,2),a);
            
            if b(2) <= 0.2
                bc(2) = 0; 
            end
            if b(2) > 0.2 && b(2) <= 0.5
                bc(2) = 1; 
            end
            if b(2) > 0.5 && b(2) <= 0.8
                bc(2) = 2; 
            end
            if b(2) > 0.8
                bc(2) = 3; 
            end
            
            %d^n
            for k=41:80
                if k-Nc(2) <= 0
                    predFrmRstd(k) = bc(2) * flippedprevresids(abs(k-Nc(2))+1); 
                else
                    predFrmRstd(k)= bc(2) * CurrFrmSTResd(k-Nc(2));
                end    
            end
            
            %e = d-d^
            for k = 41:80
                CurrFrmExFull(k) = PreCurrFrmSTResd(k) - predFrmRstd(k);
            end            

            %dequantizing bc
               if bc(2) == 0
                    bd(2) = 0.1; 
               end
               if bc(2) == 1
                    bd(2) = 0.35;  
               end
               if bc(2) == 2
                    bc(2) = 0.65; 
               end
               if bc(2) == 3
                    bd(2) = 1; 
               end
            
            
            %d'
            for k = 41:80
                if k-Nc(2) <=0
                    CurrFrmSTResd(k) = CurrFrmExFull(k) + bd(2) * flippedprevresids(abs(k-Nc(2))+1);
                else
                    CurrFrmSTResd(k) = CurrFrmExFull(k)+ bd(2) * CurrFrmSTResd(k-Nc(2));
                end        
            end
           
             a = reshape([flippedcurrResids(:,3); flippedcurrResids(:,4); prevresids(:,1)],[],1);          
            
            [Nc(3),b(3)] = RPE_subframe_LTE(resids(:,3),a);
            
            if b(3) <= 0.2
                bc(3) = 0; 
            end
            if b(3) > 0.2 && b(3) <= 0.5
                bc(3) = 1; 
            end
            if b(3) > 0.5 && b(3) <= 0.8
                bc(3) = 2; 
            end
            if b(3) > 0.8
                bc(3) = 3; 
            end
            
            %d^n
            for k=81:120
                if k-Nc(3) <= 0 
                    predFrmRstd(k) = bc(3) * flippedprevresids(abs(k-Nc(3))+1); 
                else
                    predFrmRstd(k)= bc(3) * CurrFrmSTResd(k-Nc(3));
                end    
            end
            
            %e = d-d^
            for k = 81:120
                CurrFrmExFull(k) = PreCurrFrmSTResd(k) - predFrmRstd(k);
            end            

            %dequantizing bc
               if bc(3) == 0
                    bd(3) = 0.1; 
               end
               if bc(3) == 1
                    bd(3) = 0.35;  
               end
               if bc(3) == 2
                    bc(3) = 0.65; 
               end
               if bc(3) == 3
                    bd(3) = 1; 
               end
            
            
            %d'
            for k = 81:120
                if k-Nc(3) <=0
                    CurrFrmSTResd(k) = CurrFrmExFull(k) + bd(3) * flippedprevresids(abs(k-Nc(3))+1);
                else
                    CurrFrmSTResd(k) = CurrFrmExFull(k)+ bd(3) * CurrFrmSTResd(k-Nc(3));
                end        
            end
            
            a = reshape([flippedcurrResids(:,2); flippedcurrResids(:,3); flippedcurrResids(:,4)],[],1);
            
            [Nc(4),b(4)] = RPE_subframe_LTE(resids(:,4),a);
               if b(4) <= 0.2
                    bc(4) = 0; 
               end
               if b(4) > 0.2 && b(4) <= 0.5
                    bc(4) = 1; 
               end
               if b(4) > 0.5 && b(4) <= 0.8
                    bc(4) = 2; 
               end
               if b(4) > 0.8
                    bc(4) = 3; 
               end
               
            %d^n
            for k=121:160
                predFrmRstd(k)= bc(4) * CurrFrmSTResd(k-Nc(4)); 
            end
            
            %e = d-d^
            for k = 121:160
                CurrFrmExFull(k) = PreCurrFrmSTResd(k) - predFrmRstd(k);
            end
            
            %dequantizing bc
               if bc(4) == 0
                    bd(4) = 0.1; 
               end
               if bc(4) == 1
                    bd(4) = 0.35;  
               end
               if bc(4) == 2
                    bd(4) = 0.65; 
               end
               if bc(4) == 3
                    bd(4) = 1; 
               end
            
            %d'
            for k=121:160
                CurrFrmSTResd(k) = CurrFrmExFull(k)+ bd(4) * CurrFrmSTResd(k-Nc(4));
            end
            
       
end



