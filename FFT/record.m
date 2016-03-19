 function y = record(sek,Fs)
 aud = audiorecorder(Fs,16,1);
 disp('Recording...');
 recordblocking(aud,sek);
 disp('Recording complete');
 y = getaudiodata(aud);
 end