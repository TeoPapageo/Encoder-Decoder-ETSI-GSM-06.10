    function [N,b] = RPE_subframe_LTE(d,Prevd)

    
%vlepo pos tha exo 80 times sto R
R = zeros(80, 1);
%theto lamda = 40
lamda = 40;
%kai fysika arxikopoio to N iso me tin arxiki timi tou lamda
N = 40;
max = 0;
%ypologismos N
for i = 1:80
    for k = 1:40
        R(i) = R(i) + d(k)*Prevd(lamda-k+1);
    end
    if i ~= 1
      if R(i) > max
          N = lamda;
          R(i) = max;
      end
    else
        R(i) = max;
    end    
    lamda = lamda + 1;
  
end



%ypologismos b

c = 0;
for i = 1:40 
   c = c + d(i)*Prevd(N-k+1); 
end
d = 0;
for i = 1:40 
   d = d + (Prevd(N-k+1))^2; 
end

b = c/d;



end

