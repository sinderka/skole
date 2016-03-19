function [Ax,Ay,Bx,By] = setup_lhs(nx,ny,hx,hy)
Ax = kron(speye(ny)*hy,K1(nx-1,hx,2))+kron(K1(ny,hy,3),speye(nx-1)*hx);
Ay = kron(speye(ny-1)*hy,K1(nx,hx,3))+kron(K1(ny-1,hy,2),speye(nx)*hx);
Bx = kron(speye(ny)*hy,K2(nx));
By = kron(K2(ny),speye(nx)*hx);