%here's my little program to demonstrate that my 2 functions are working 
%properly

%first of all, I'm saving a voice .wav file with Fs = 8 kHz that I've
%found on the internet. Samples are on variable y.

filename = 'male.wav';  
[y, Fs] = audioread(filename);

%then of course I don't know if my samples are dividable by 160 so I have 
%to make sure that it is

k = fix(size(y)/160); 

%d is going to be the output signal from my decoder. That's why it has to
%be dividable by 160.
d = zeros(k(1)*160,1);


%a little explanation of my code: z0 and d0 are set in each iteration as
%the next 160 samples of the .wav file, saved on y. each time I add the
%output of my RPE_frame_ST_decoder function to d
z0 = zeros(160, 1);
d0 = zeros(160, 1);
for i = 1:k(1)
    
     z0 = y(1 + (160*(i-1)) : 160*(i-1) + 160  );
     
     [a, b] = RPE_frame_ST_coder(z0);
     d0 = RPE_frame_ST_decoder(a, b);
     
     d(1 + (160*(i-1)) :  160*(i-1) + 160) = d0;
    
end

%now d is my decoded voice signal. enjoy :)

sound(d,8000);


