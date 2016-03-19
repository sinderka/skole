% Make sure FFT is set correctly 
FFT

disp('e: 329.6Hz')
disp('B: 246.9Hz')
disp('G: 196.0Hz')
disp('D: 146.8Hz')
disp('A: 110.0Hz')
disp('E:  82.4Hz')

disp('RealFreq:')
freq = 0;
a = max(abs(real(X(60:end))));
for i = 60:1:length(X)
  if  abs(real(X(i))) > a*0.5
      freq = i
  end
end

disp('ImagFreq:')
freq = 0;
a = max(abs(imag(X(60:end))));
for i = 60:1:length(X)
  if  abs(real(X(i))) > a*0.5
      freq = i
  end
end



