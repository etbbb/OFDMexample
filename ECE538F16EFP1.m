%%ECE538F22016EFP1 : OFDM example

y = [13-4j 19+j 10+4j 7-j]
h = [2 4 1]

s3 = (1/4)*exp(j*2*pi*(3/4).*(0:3));
s2 = (1/4)*exp(j*2*pi*(2/4).*(0:3));
s1 = (1/4)*exp(j*2*pi*(1/4).*(0:3));
s0 = (1/4)*exp(j*2*pi*(0/4).*(0:3));

H = [2 2 2 2] + 4*exp(-j*2*pi*(1/4).*(0:3)) + exp(-2*j*2*pi*(1/4).*(0:3));
HH = fft(h,4);
delll = H - HH

b0 = 4*(1/H(0+1))*y*s0'
b1 = 4*(1/H(1+1))*y*s1'
b2 = 4*(1/H(2+1))*y*s2'
b3 = 4*(1/H(3+1))*y*s3'

disp(delll)
function newy = dumpCP(y,h)
   L = length(h)
   newy = y(L+1:end-L+1)
end

function H = HHH(h,N)
   %Need just the first bit of the data for it to work
   %The question is how I can adjust the numerical fft in order to
   %get the data working properly. 
   H = fft(h,N)
end

function bk = symbolsOFDM(y,h,N)
   yy = dumpCP(y,h)
   H  = HHH(h,N)
   bktemp = [];
   for u = (0:(N-1))
      bb = N*(1/H(u+1))*yy*(1/N)*exp(j*2*pi*(u/N).*(0:(N-1)))'
      bktemp = [bktemp bb];
   end
   bk = bktemp;
end

ytest = [2+2j, 6+4j, 7-1j, 13-4j, 19+1j, 10+4j, 7-1j, 5-4j, 1-1j];
bkbk = (real(symbolsOFDM(ytest,h,4)))
