function A = A_Matrix(A,CD,H,J2,Re,m,mu,r0,rho0,theta_dot,x,x_dot,y,y_dot,z,z_dot)
%A_Matrix
%    A = A_Matrix(A,CD,H,J2,Re,M,MU,R0,RHO0,THETA_DOT,X,X_DOT,Y,Y_DOT,Z,Z_DOT)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    28-Oct-2023 18:26:36

t2 = abs(x);
t3 = abs(y);
t4 = abs(z);
t5 = abs(z_dot);
t6 = sign(x);
t7 = sign(y);
t8 = sign(z);
t9 = sign(z_dot);
t10 = theta_dot.*x;
t11 = theta_dot.*y;
t12 = Re.^2;
t13 = z.^2;
t14 = z.^3;
t19 = 1.0./H;
t20 = 1.0./m;
t15 = t2.^2;
t16 = t3.^2;
t17 = t4.^2;
t18 = t5.^2;
t21 = t11+x_dot;
t22 = -t10;
t23 = abs(t21);
t24 = sign(t21);
t25 = t22+y_dot;
t30 = t15+t16+t17;
t26 = abs(t25);
t27 = sign(t25);
t28 = t23.^2;
t31 = 1.0./t30;
t33 = sqrt(t30);
t29 = t26.^2;
t32 = t31.^2;
t34 = 1.0./t33;
t40 = t31.*z.*1.0e+1;
t41 = -t33;
t42 = t13.*t31.*5.0;
t35 = t34.^3;
t36 = t34.^5;
t37 = t34.^7;
t44 = r0+t41;
t45 = t18+t28+t29;
t46 = t42-1.0;
t48 = t4.*t8.*t13.*t32.*1.0e+1;
t38 = t35.^3;
t39 = mu.*t35;
t47 = t46-2.0;
t49 = -t48;
t50 = t19.*t44;
t51 = sqrt(t45);
t55 = J2.*mu.*t12.*t36.*t46.*(3.0./2.0);
t43 = -t39;
t52 = exp(t50);
t53 = 1.0./t51;
t54 = t40+t49;
t56 = (A.*CD.*rho0.*t20.*t51.*t52)./2.0;
t57 = t56.*theta_dot;
t58 = -t56;
mt1 = [0.0,0.0,0.0,t43+t55+mu.*t2.*t6.*t36.*x.*3.0+t2.*t6.*t19.*t21.*t34.*t56-J2.*mu.*t2.*t6.*t12.*t13.*t38.*x.*1.5e+1-J2.*mu.*t2.*t6.*t12.*t37.*t46.*x.*(1.5e+1./2.0)+(A.*CD.*rho0.*t20.*t21.*t26.*t27.*t52.*t53.*theta_dot)./2.0,t57+mu.*t2.*t6.*t36.*y.*3.0-J2.*mu.*t2.*t6.*t12.*t13.*t38.*y.*1.5e+1-J2.*mu.*t2.*t6.*t12.*t37.*t46.*y.*(1.5e+1./2.0)-(A.*CD.*rho0.*t20.*t26.*t27.*t52.*t53.*theta_dot.*(t10-y_dot))./2.0-(A.*CD.*rho0.*t2.*t6.*t19.*t20.*t34.*t51.*t52.*(t10-y_dot))./2.0,mu.*t2.*t6.*t36.*z.*3.0+t2.*t6.*t19.*t34.*t56.*z_dot-J2.*mu.*t2.*t6.*t12.*t14.*t38.*1.5e+1-J2.*mu.*t2.*t6.*t12.*t37.*t47.*z.*(1.5e+1./2.0)+(A.*CD.*rho0.*t20.*t26.*t27.*t52.*t53.*theta_dot.*z_dot)./2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
mt2 = [0.0,0.0,0.0,0.0,0.0,mu.*t3.*t7.*t36.*x.*3.0+t3.*t7.*t19.*t21.*t34.*t56-(A.*CD.*rho0.*t20.*t51.*t52.*theta_dot)./2.0-J2.*mu.*t3.*t7.*t12.*t13.*t38.*x.*1.5e+1-J2.*mu.*t3.*t7.*t12.*t37.*t46.*x.*(1.5e+1./2.0)-(A.*CD.*rho0.*t20.*t21.*t23.*t24.*t52.*t53.*theta_dot)./2.0,t43+t55+mu.*t3.*t7.*t36.*y.*3.0-J2.*mu.*t3.*t7.*t12.*t13.*t38.*y.*1.5e+1-J2.*mu.*t3.*t7.*t12.*t37.*t46.*y.*(1.5e+1./2.0)+(A.*CD.*rho0.*t20.*t23.*t24.*t52.*t53.*theta_dot.*(t10-y_dot))./2.0-(A.*CD.*rho0.*t3.*t7.*t19.*t20.*t34.*t51.*t52.*(t10-y_dot))./2.0];
mt3 = [mu.*t3.*t7.*t36.*z.*3.0+t3.*t7.*t19.*t34.*t56.*z_dot-J2.*mu.*t3.*t7.*t12.*t14.*t38.*1.5e+1-J2.*mu.*t3.*t7.*t12.*t37.*t47.*z.*(1.5e+1./2.0)-(A.*CD.*rho0.*t20.*t23.*t24.*t52.*t53.*theta_dot.*z_dot)./2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,mu.*t4.*t8.*t36.*x.*3.0+t4.*t8.*t19.*t21.*t34.*t56+J2.*mu.*t12.*t36.*t54.*x.*(3.0./2.0)-J2.*mu.*t4.*t8.*t12.*t37.*t46.*x.*(1.5e+1./2.0),mu.*t4.*t8.*t36.*y.*3.0+J2.*mu.*t12.*t36.*t54.*y.*(3.0./2.0)-J2.*mu.*t4.*t8.*t12.*t37.*t46.*y.*(1.5e+1./2.0)-(A.*CD.*rho0.*t4.*t8.*t19.*t20.*t34.*t51.*t52.*(t10-y_dot))./2.0];
mt4 = [t43+J2.*mu.*t12.*t36.*t47.*(3.0./2.0)+mu.*t4.*t8.*t36.*z.*3.0+t4.*t8.*t19.*t34.*t56.*z_dot+J2.*mu.*t12.*t36.*t54.*z.*(3.0./2.0)-J2.*mu.*t4.*t8.*t12.*t37.*t47.*z.*(1.5e+1./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,t58-(A.*CD.*rho0.*t20.*t21.*t23.*t24.*t52.*t53)./2.0,(A.*CD.*rho0.*t20.*t23.*t24.*t52.*t53.*(t10-y_dot))./2.0,A.*CD.*rho0.*t20.*t23.*t24.*t52.*t53.*z_dot.*(-1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,A.*CD.*rho0.*t20.*t21.*t26.*t27.*t52.*t53.*(-1.0./2.0),t58+(A.*CD.*rho0.*t20.*t26.*t27.*t52.*t53.*(t10-y_dot))./2.0,A.*CD.*rho0.*t20.*t26.*t27.*t52.*t53.*z_dot.*(-1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0];
mt5 = [A.*CD.*rho0.*t5.*t9.*t20.*t21.*t52.*t53.*(-1.0./2.0),(A.*CD.*rho0.*t5.*t9.*t20.*t52.*t53.*(t10-y_dot))./2.0,t58-(A.*CD.*rho0.*t5.*t9.*t20.*t52.*t53.*z_dot)./2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t35.*x+J2.*t12.*t36.*t46.*x.*(3.0./2.0),-t35.*y+J2.*t12.*t36.*t46.*y.*(3.0./2.0),-t35.*z+J2.*t12.*t36.*t47.*z.*(3.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,mu.*t12.*t36.*t46.*x.*(3.0./2.0),mu.*t12.*t36.*t46.*y.*(3.0./2.0),mu.*t12.*t36.*t47.*z.*(3.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,A.*rho0.*t20.*t21.*t51.*t52.*(-1.0./2.0),(A.*rho0.*t20.*t51.*t52.*(t10-y_dot))./2.0];
mt6 = [A.*rho0.*t20.*t51.*t52.*z_dot.*(-1.0./2.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
A = reshape([mt1,mt2,mt3,mt4,mt5,mt6],18,18);