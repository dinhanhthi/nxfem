tx=zeros(1,nStep); ty=zeros(4,nStep);
    for i=1:nStep
       tx(i)=msh(i).hTmax;
       ty(1,i)=err(i).L2;
%        ty(2,i)=err(i).H1;
%        ty(3,i)=err(i).L2G;
       ty(4,i)=err(i).ENorm;
    end
    tmp = polyfit(log(tx),log(ty(1,:)),1);
    order.L2 = tmp(1);
%     tmp = polyfit(log(tx),log(ty(2,:)),1);
%     order.H1 = tmp(1);
%     tmp = polyfit(log(tx),log(ty(3,:)),1);
%     order.L2G = tmp(1);
    tmp = polyfit(log(tx),log(ty(4,:)),1);
    order.ENorm = tmp(1);
    
    
    plot(log(tx),log(ty(1,:)),'-r.',...
            log(tx),log(ty(4,:)),'-b.');
        legend('L2-norm','H-norm');
        xlabel('log(h)'); 
        ylabel('log(error)');