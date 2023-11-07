function [ux,uy]=gradu(u,h,ind)
% compute gradient
% ind: 1-forward; 2-backward

uex=expandf(u);

switch ind
    case 1
        ux=(uex(3:end,2:end-1)-uex(2:end-1,2:end-1))/h;
        uy=(uex(2:end-1,3:end)-uex(2:end-1,2:end-1))/h;
        
    case 2
        ux=(uex(2:end-1,2:end-1)-uex(1:end-2,2:end-1))/h;
        uy=(uex(2:end-1,2:end-1)-uex(2:end-1,1:end-2))/h;
    case 0
        ux=(uex(3:end,2:end-1)-uex(1:end-2,2:end-1))/(2*h);
        uy=(uex(2:end-1,3:end)-uex(2:end-1,1:end-2))/(2*h);
end