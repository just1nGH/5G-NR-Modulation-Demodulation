function [] = plotNrLLR(moduType)

close all;
switch lower(moduType)
    case 'bpsk'
        A =1;
        K =1;
        maxX = 2;
        maxXtick = 1;
        figwidth = 250;
        figHeight = 250;
    case 'qpsk'
        A = 1/sqrt(2);
        K = 2;
        maxX = 2;
        maxXtick = 1;
        figwidth = 500;
        figHeight = 250;

    case '16qam'
        A = 1/sqrt(10);
        K = 4;
        maxX = 6;
        maxXtick = 4;
        figwidth = 500;
        figHeight = 500;
    case  '64qam'
        A = 1/sqrt(42);
        K = 6;
        maxX = 9;
        maxXtick = 7;
        figwidth = 500;
        figHeight = 800;
    case '256qam'
        A = 1/sqrt(170);
        K = 8;
        maxX = 16;
        maxXtick = 15;
        figwidth = 800;
        figHeight = 1100;

end


x = (-maxX:maxX)*A;

% max-log-map
y1 = nrSoftModuDemapper(x+1i*x,lower(moduType),1,'max-log-map');
y1 = reshape(y1,K,[]);

% linear approximation
y2 = nrSoftModuDemapper(x+1i*x,lower(moduType),1,'approx');
y2 = reshape(y2,K,[]);
out = nrSymbolDemodulate((x+1i*x).',moduType,1);
out = reshape(out,K,[]);
figure;

for i = 1:K
    if(K >1)
        subplot(K/2,2,i);
    end
    plot(x/A, y1(i,:),'r-');
    
    grid on;
    hold on;
    plot(x/A, out(i,:),'ro');
    plot(x/A, y2(i,:),'b--');
    %plot(x(2:end-1)/A, y1(i,2:end-1),'ro');

    if mod(i,2)
        xlabel({'Re(x)','(A)'});
    else
        xlabel({'Im(x)','(A)'});
    end
    if (K >1)
        ylabel({['LLR_' int2str(i)], '(1/N0)'});
    else
        ylabel({'LLR', '(1/N_0)'});
    end

    xticks(-maxXtick:2:maxXtick);
    xlim([-maxXtick-0.9,maxXtick+0.9]);
    if i == 1
        legend('max-log-map','matlab nr tool box','linear approx.','Location','north');
    end
end



set(gcf,'Position',[100 100 figwidth figHeight])
sgtitle([moduType ' LLR'])

end

