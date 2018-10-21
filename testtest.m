% % system('/users/home/blouza/MshDist/build/mshdist mshdist/phi')
% system('mshdist mshdist/phi');
% 
% NT = 5;
% idx=0;
% CP = zeros(1,NT);
% CP(1:NT) = CP(1:NT) + 1;
% CP

    
% CP = zeros(1,nCT);
phiTri = zeros(3,nCT);
phiTri = phi(CTs(1,1:nCT)) .* phi(CTs(2,1:nCT)) <0... 
            & abs(phi(CTs(1,1:nCT)))>tole...
            & abs(phi(CTs(2,1:nCT)))>tole



% eee = phiTri(1,1:nCT) < 0


function check = equal0(num1,num2,tole)
    if (abs(num1)<tole)||(abs(num2)<tole)
        check=1;
    else
        check=0;
    end
end
function check = negative(num1,num2,tole)
    if (num1*num2<0)&&(abs(num1)>tole)&&(abs(num2)>tole)
        check=1; 
    else
        check=0;
    end
end