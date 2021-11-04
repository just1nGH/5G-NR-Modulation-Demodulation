function softBits = nrSoftModuDemapper(symbsIn,moduType,N0,method)
% nrSoftModuDemapper demodulates symbols to soft bits
%   softBits = nrSoftModuDemapper(symbsIn,moduType,N0,method) demodulates the complex symbols
%   symbsIn using soft decision. The modulation scheme, moduType must be one
%   of  'BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'. symbsIn must be
%   a column vector.
%   'mothod'   -   Specified as 'max-log-map' or 'approx'.
%
%   Demodulation is performed according to the constellations given in
%   TS 38.211 section 5.1 including the power normalization factors
%   specified. The factors: 1/sqrt(2) for BPSK, and QPSK,
%   1/sqrt(10) for 16QAM, 1/sqrt(42) for 64QAM and 1/sqrt(170) for 256QAM.
%
%   Source Paper: Mao, Juquan, et al. "A low complexity 256QAM soft demapper for 5G mobile system." 
%   2016 European Conference on Networks and Communications (EuCNC). IEEE, 2016.
%   Author Dr J Mao
%   Email: juquan.justin.mao@gmail.com
%   2021 Nov
%   Example:
%   %Demonstrate 16-QAM  demodulation.
%   msg = randi([0 1],100,1,'int8');
%   symb = nrModuMapper(msg,'16QAM');
%   N0 = 0.1;
%   rxsymb = symb + sqrt(N0/2)*randn(size(symb)) ;
%   msg_hat = nrSoftModuDemapper(msg_hat,'16QAM',N0,'max-log-map');
%   numErr = sum(msg ~=(msg_hat < 0));



    switch lower(moduType)
        case 'bpsk'
            softBits = 4/N0 * real(exp(-1j*pi/4)*symbsIn);
        case 'qpsk'
            softBits = zeros(2,length(symbsIn));
            A = 1/sqrt(2);
            softBits(1,:) = 4*A/N0 * real(symbsIn);
            softBits(2,:) = 4*A/N0 * imag(symbsIn);
            softBits = softBits(:);
        case '16qam'
            softBits = qam16SoftDemodulation(symbsIn,N0,method);
        case  '64qam'
            softBits = qam64SoftDemodulation(symbsIn,N0,method);
        case '256qam'
            softBits = qam256SoftDemodulation(symbsIn,N0,method);
    end

end

%%------------------------------------------------------------
% 16QAM
function softBits = qam16SoftDemodulation(symbsIn,N0,method)

    softBits = zeros(4,length(symbsIn));
    
    A = 1/sqrt(10);
    reX = real(symbsIn);
    imX =  imag(symbsIn);
    if strcmpi(method,'approx') 
        softBits(1,:) = 4*A/N0 * reX;
        softBits(2,:) = 4*A/N0 * imX;
    else % max-log-MAP
        softBits(1,:) = maxLogMapQam16SoftDemodulation_1_2(reX,N0);
        softBits(2,:) = maxLogMapQam16SoftDemodulation_1_2(imX,N0);
    end
    
    softBits(3,:) = 4*A/N0 * (-abs(reX)+ 2*A);
    softBits(4,:) = 4*A/N0 * (-abs(imX)+ 2*A);
    
    softBits = softBits(:);

end


function  softBits = maxLogMapQam16SoftDemodulation_1_2(realSymbIn, N0)

    A = 1/sqrt(10);
    X = realSymbIn;
    Y = zeros(size(X));
    for i = 1: length(X)
        if X(i) < -2*A
            Y(i) = 8*A/N0 * (X(i) + A);
        elseif X(i) < 2*A
            Y(i) = 4*A/N0 * X(i);
        else
            Y(i) = 8*A/N0 * (X(i) - A);
        end
             
    end
    softBits = Y;

end

%%------------------------------------------------------------
% 16QAM
function softBits = qam64SoftDemodulation(symbsIn,N0,method)

    softBits = zeros(4,length(symbsIn));
    
    A = 1/sqrt(42);
    reX = real(symbsIn);
    imX =  imag(symbsIn);
    
    if strcmpi(method,'approx') 
        softBits(1,:) = 4*A/N0 * reX;
        softBits(2,:) = 4*A/N0 * imX;
        softBits(3,:) = 4*A/N0 * (-abs(reX) + 4*A);
        softBits(4,:) = 4*A/N0 * (-abs(imX) + 4*A);
    else % max-log-MAP
        softBits(1,:) = maxLogMapQam64SoftDemodulation_1_2(reX, N0);
        softBits(2,:) = maxLogMapQam64SoftDemodulation_1_2(imX, N0);
        softBits(3,:) = maxLogMapQam64SoftDemodulation_3_4(reX, N0);
        softBits(4,:) = maxLogMapQam64SoftDemodulation_3_4(imX, N0);
    end
    
    softBits(5,:) = 4*A/N0 * (-abs(-abs(reX)+ 4*A)+2*A);
    softBits(6,:) = 4*A/N0 * (-abs(-abs(imX)+ 4*A)+2*A);
    
    softBits = softBits(:);

end
function  softBits = maxLogMapQam64SoftDemodulation_1_2(realSymbIn, N0)

    A = 1/sqrt(42);
    X = realSymbIn;
    Y = zeros(size(X));
    for i = 1: length(X)
        if X(i) < -6*A
            Y(i) = 16*A/N0 * (X(i) + 3*A);
        elseif X(i) < -4*A
            Y(i) = 12*A/N0 * (X(i) + 2*A);
        elseif X(i) < -2*A
            Y(i) = 8*A/N0 * (X(i) + A);
        elseif X(i) < 2*A
            Y(i) = 4*A/N0 * X(i) ;
        elseif X(i) < 4*A
            Y(i) = 8*A/N0 * (X(i) - A);
        elseif X(i) < 6*A
            Y(i) = 12*A/N0 * (X(i) - 2*A);
        else
            Y(i) = 16*A/N0 * (X(i) - 3*A);
        end
             
    end
    softBits = Y;

end

function  softBits = maxLogMapQam64SoftDemodulation_3_4(realSymbIn, N0)

    A = 1/sqrt(42);
    X = realSymbIn;
    Y = zeros(size(X));
    for i = 1: length(X)
        if X(i) < -6*A
            Y(i) = 8*A/N0 * (X(i) + 5*A);
        elseif X(i) < -2*A
            Y(i) = 4*A/N0 * (X(i) + 4*A);
        elseif X(i) < 0
            Y(i) = 8*A/N0 * (X(i) + 3*A);
        elseif X(i) < 2*A
            Y(i) = 8*A/N0 * (-X(i) + 3*A);
        elseif X(i) < 6*A
            Y(i) = 4*A/N0 * (-X(i) + 4*A);
        else
            Y(i) = 8*A/N0 * (-X(i) + 5*A);
        end     
    end
    softBits = Y;

end

%%-----------------------------------------------------------
% 256QAM
function softBits = qam256SoftDemodulation(symbsIn,N0,method)

    softBits = zeros(8,length(symbsIn));
    
    A = 1/sqrt(170);
    reX = real(symbsIn);
    imX =  imag(symbsIn);
    
    if strcmpi(method,'approx') 
        softBits(1,:) = 4*A/N0 * reX;
        softBits(2,:) = 4*A/N0 * imX;
        softBits(3,:) = 4*A/N0 * (-abs(reX) + 8*A);
        softBits(4,:) = 4*A/N0 * (-abs(imX) + 8*A);
        softBits(5,:) = 4*A/N0 * (-abs(-abs(reX)+8*A) + 4*A);
        softBits(6,:) = 4*A/N0 * (-abs(-abs(imX)+8*A) + 4*A);
    else % max-log-MAP
        softBits(1,:) = maxLogMapQam256SoftDemodulation_1_2(reX, N0);
        softBits(2,:) = maxLogMapQam256SoftDemodulation_1_2(imX, N0);
        softBits(3,:) = maxLogMapQam256SoftDemodulation_3_4(reX, N0);
        softBits(4,:) = maxLogMapQam256SoftDemodulation_3_4(imX, N0);
        softBits(5,:) = maxLogMapQam256SoftDemodulation_5_6(reX, N0);
        softBits(6,:) = maxLogMapQam256SoftDemodulation_5_6(imX, N0);
    end
    
    softBits(7,:) = 4*A/N0 * (-abs(-abs(-abs(reX)+8*A) + 4*A)+2*A);
    softBits(8,:) = 4*A/N0 * (-abs(-abs(-abs(imX)+8*A) + 4*A)+2*A);
    
    softBits = softBits(:);

end
function  softBits = maxLogMapQam256SoftDemodulation_1_2(realSymbIn, N0)
    A = 1/sqrt(170);
    X = realSymbIn;
    Y = zeros(size(X));
    for i = 1: length(X)
        if X(i) < -14*A
            Y(i) = 32*A/N0 * (X(i) + 7*A);
        elseif X(i) < -12*A
            Y(i) = 28*A/N0 * (X(i) + 6*A);
        elseif X(i) < -10*A
            Y(i) = 24*A/N0 * (X(i) + 5*A);
        elseif X(i) < -8*A
            Y(i) = 20*A/N0 * (X(i) + 4*A);
        elseif X(i) < -6*A
            Y(i) = 16*A/N0 * (X(i) + 3*A);
        elseif X(i) < -4*A
            Y(i) = 12*A/N0 * (X(i) + 2*A);
        elseif X(i) < -2*A
            Y(i) = 8*A/N0 * (X(i) + A);
        elseif X(i) < 2*A
            Y(i) = 4*A/N0 * X(i);
        elseif X(i) < 4*A
            Y(i) = 8*A/N0 * (X(i) - A);
        elseif X(i) < 6*A
            Y(i) = 12*A/N0 * (X(i) - 2*A);
        elseif X(i) < 8*A
            Y(i) = 16*A/N0 * (X(i) - 3*A);
        elseif X(i) < 10*A
            Y(i) = 20*A/N0 * (X(i) - 4*A);
        elseif X(i) < 12*A
            Y(i) = 24*A/N0 * (X(i) - 5*A);
        elseif X(i) < 14*A
            Y(i) = 28*A/N0 * (X(i) - 6*A);
        else
            Y(i) = 32*A/N0 * (X(i) - 7*A);
        end
             
    end
    softBits = Y;

end

function  softBits = maxLogMapQam256SoftDemodulation_3_4(realSymbIn, N0)

    A = 1/sqrt(170);
    X = realSymbIn;
    Y = zeros(size(X));
    for i = 1: length(X)
        if X(i) < -14*A
            Y(i) = 16*A/N0 * (X(i) + 11*A);
        elseif X(i) < -12*A
            Y(i) = 12*A/N0 * (X(i) + 10*A);
        elseif X(i) < -10*A
            Y(i) = 8*A/N0 * (X(i) + 9*A);
        elseif X(i) < -6*A
            Y(i) = 4*A/N0 * (X(i) + 8*A);
        elseif X(i) < -4*A
            Y(i) = 8*A/N0 * (X(i) + 7*A);
        elseif X(i) < -2*A
            Y(i) = 12*A/N0 * (X(i) + 6*A);
        elseif X(i) < 0
            Y(i) = 16*A/N0 * (X(i) + 5*A);
        elseif X(i) < 2*A
            Y(i) = 16*A/N0 * (-X(i) + 5*A);
        elseif X(i) < 4*A
            Y(i) = 12*A/N0 * (-X(i) + 6*A);
        elseif X(i) < 6*A
            Y(i) = 8*A/N0 * (-X(i) + 7*A);
        elseif X(i) < 10*A
            Y(i) = 4*A/N0 * (-X(i) + 8*A);
        elseif X(i) < 12*A
            Y(i) = 8*A/N0 * (-X(i) + 9*A);
        elseif X(i) < 14*A
            Y(i) = 12*A/N0 * (-X(i) + 10*A);
        else
            Y(i) = 16*A/N0 * (-X(i) + 11*A);
        end
             
    end
    softBits = Y;


end

function  softBits = maxLogMapQam256SoftDemodulation_5_6(realSymbIn, N0)

    A = 1/sqrt(170);
    X = realSymbIn;
    Y = zeros(size(X));
    for i = 1: length(X)
        if X(i) < -14*A
            Y(i) = 8*A/N0 * (X(i) + 13*A);
        elseif X(i) < -10*A
            Y(i) = 4*A/N0 * (X(i) + 12*A);
        elseif X(i) < -8*A
            Y(i) = 8*A/N0 * (X(i) + 11*A);
        elseif X(i) < -6*A
            Y(i) = 8*A/N0 * (-X(i) - 5*A);
        elseif X(i) < -2*A
            Y(i) = 4*A/N0 * (-X(i) - 4*A);
        elseif X(i) < 0
            Y(i) = 8*A/N0 * (-X(i) - 3*A);
        elseif X(i) < 2*A
            Y(i) = 8*A/N0 * (X(i) - 3*A);
        elseif X(i) < 6*A
            Y(i) = 4*A/N0 * (X(i) - 4*A);
        elseif X(i) < 8*A
            Y(i) = 8*A/N0 * (X(i) - 5*A);
        elseif X(i) < 10*A
            Y(i) = 8*A/N0 * (-X(i) + 11*A);
        elseif X(i) < 14*A
            Y(i) = 4*A/N0 * (-X(i) + 12*A);
        else
            Y(i) = 8*A/N0 * (-X(i) + 13*A);
        end
             
    end
    softBits = Y;


end





