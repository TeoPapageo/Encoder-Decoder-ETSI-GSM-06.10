function [s0, CurrFrmSTResd] = RPE_frame_SLT_decoder(LARc,Nc,bc,CurrFrmExFull, PrevFrmSTResd)
bd = zeros(4,1); %dequantized b 
CurrFrmSTResd = zeros(160,1); %d'
flippedprevresids = flip(PrevFrmSTResd);

for i = 1:4
%dequantizing bc
     if bc(i) == 0
        bd(i) = 0.1; 
     end
     if bc(i) == 1
        bd(i) = 0.35;  
     end
     if bc(i) == 2
        bc(i) = 0.65; 
     end
     if bc(i) == 3
        bd(i) = 1; 
     end    
end
    
    
            for k=1:40
                CurrFrmSTResd(k) = CurrFrmExFull(k) + bd(1)*flippedprevresids(Nc(1)+1-k);
            end            
            
            for k = 41:80
                if k-Nc(2) <=0
                    CurrFrmSTResd(k) = CurrFrmExFull(k) + bd(2) * flippedprevresids(abs(k-Nc(2))+1);
                else
                    CurrFrmSTResd(k) = CurrFrmExFull(k)+ bd(2) * CurrFrmSTResd(k-Nc(2));
                end        
            end            
            
           
            for k = 81:120
                if k-Nc(3) <=0
                    CurrFrmSTResd(k) = CurrFrmExFull(k) + bd(3) * flippedprevresids(abs(k-Nc(3))+1);
                else
                    CurrFrmSTResd(k) = CurrFrmExFull(k)+ bd(3) * CurrFrmSTResd(k-Nc(3));
                end        
            end            
            
            
            for k=121:160
                CurrFrmSTResd(k) = CurrFrmExFull(k)+ bd(4) * CurrFrmSTResd(k-Nc(4));
            end            


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

