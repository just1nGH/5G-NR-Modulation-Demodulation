function [] = plotNrConstellation(moduType)


close all;
switch lower(moduType)
    case 'bpsk'
        A =1;
        M =2;
        xMax = 1.2;
        figHeight = 250;
    case 'qpsk'
        A = 1/sqrt(2);
        M = 4;
        xMax = 1.2;
        figHeight = 250;

    case '16qam'
        A = 1/sqrt(10);
        M = 16;
        xMax = 3.2;
        figHeight = 600;
    case  '64qam'
        A = 1/sqrt(42);
        M = 64;
        xMax = 7.5;
        figHeight = 800;
    case '256qam'
        A = 1/sqrt(170);
        M = 256;
        xMax = 15.9;
        figHeight = 1200;

end

% 256QAM
symbBits = de2bi(0:M-1,'left-msb');
symbBitsIn = symbBits.';
symbs = nrModuMapper(symbBitsIn(:),lower(moduType));
% close all;
close all;
plot(symbs/A,'or','MarkerFaceColor','r');
hold on;

for  i = 1: M
    str = num2str(symbBits(i,:));
    str = str(find(~isspace(str)));
    text(real(symbs(i))/A, imag(symbs(i))/A -0.3,str ,'HorizontalAlignment' ,'center');
    %text(real(symbs(i)), imag(symbs(i))-0.05,num2str(i-1),'color','r' );
end

xlim([-xMax, xMax]); ylim([-xMax-0.3,xMax])
xlabel({'Inphase','(A)'});
ylabel({'Qadrature','(A)'});
xticks([-floor(xMax):2:-1 0 1:2:floor(xMax)]); yticks([-floor(xMax):2:-1 0 1:2:floor(xMax)]);
grid on;
set(gcf,'Position',[100 100 figHeight figHeight])
title([moduType ' Constellation'])
end

