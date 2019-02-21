function [Jacobian]=fekine2d(Xcoord,Ycoord,ShapFun)
syms 'x';syms 'y';
CoordMatr=zeros(4,2);
for i=1:4
CoordMatr(i,1)=Xcoord(i);
CoordMatr(i,2)=Ycoord(i);
end
[Jacobian]=[diff(ShapFun(1),x),diff(ShapFun(2),x),diff(ShapFun(3),x),diff(ShapFun(4),x);
    diff(ShapFun(1),y),diff(ShapFun(2),y),diff(ShapFun(3),y),diff(ShapFun(4),y)]*CoordMatr;
end
