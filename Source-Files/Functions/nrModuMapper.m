function symbOut= nrModuMapper(bitsIn,moduType)
% symbOut = nrModuMapper(bitsIn,moduType) maps the bit sequence into 
% complex modulation symbols using modulation scheme specified in TS 38.211 Section 5.1
% The modulation scheme, moduType must be one of  'BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'.

b = bitsIn;
switch lower(moduType)
    case 'bpsk'
        symbOut = zeros(size(b));
        for i = 0:length(bitsIn)-1
           symbOut(i+1) =  0.5*sqrt(2)*((1 - 2*b(i+1)) + 1j*(1-2*b(i+1)));
        end
        
    case 'qpsk'
        nOfBitsIn = length(bitsIn);
        assert(mod(nOfBitsIn,2) == 0);
        symbOut = zeros(nOfBitsIn/2,1);
        for i = 0: nOfBitsIn/2-1
            symbOut(i+1) = 1/sqrt(2)*((1-2*b(2*i+1))+1j*(1-2*b(2*i+2)));
        end
        
    case '16qam'
        
        nOfBitsIn = length(bitsIn);
        assert(mod(nOfBitsIn,4) == 0);
        symbOut = zeros(nOfBitsIn/4,1);
        for i = 0: nOfBitsIn/4-1
            symbOut(i+1) = (1-2*b(4*i+1))*(2-(1-2*b(4*i+3)))+...
                          1j*((1-2*b(4*i+2))*(2-(1-2*b(4*i+4)))); 
        end
        symbOut = symbOut/sqrt(10);
        
        
    case  '64qam'
        nOfBitsIn = length(bitsIn);
        assert(mod(nOfBitsIn,6) == 0);
        symbOut = zeros(nOfBitsIn/6,1);
        for i = 0: nOfBitsIn/6-1
            symbOut(i+1) = (1-2*b(6*i+1))*(4-(1-2*b(6*i+3))*(2-(1-2*b(6*i+5))))+...
                          1j*(1-2*b(6*i+2))*(4-(1-2*b(6*i+4))*(2-(1-2*b(6*i+6))));
        end
         symbOut = symbOut/sqrt(42);
         
         
    case '256qam'
        nOfBitsIn = length(bitsIn);
        assert(mod(nOfBitsIn,8) == 0);
        symbOut = zeros(nOfBitsIn/8,1);
        for i = 0: nOfBitsIn/8-1
            symbOut(i+1) = (1-2*b(8*i+1))*(8-(1-2*b(8*i+3))*(4-(1-2*b(8*i+5))*(2-(1-2*b(8*i+7)))))+...
                          1j*(1-2*b(8*i+2))*(8-(1-2*b(8*i+4))*(4-(1-2*b(8*i+6))*(2-(1-2*b(8*i+8))))); 
        end
        symbOut = symbOut/sqrt(170);
end



end

