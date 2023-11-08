function H_tilde = H_tilde_xs1(theta,theta_dot,x,x_dot,xs1,y,y_dot,ys1,z,z_dot,zs1)
%H_tilde_xs1
%    H_tilde = H_tilde_xs1(THETA,THETA_DOT,X,X_DOT,XS1,Y,Y_DOT,YS1,Z,Z_DOT,ZS1)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    31-Oct-2023 11:07:41

t2 = cos(theta);
t3 = sin(theta);
t4 = x.*x_dot;
t5 = x.*xs1;
t6 = x_dot.*xs1;
t7 = x.*ys1;
t8 = xs1.*y;
t9 = x_dot.*ys1;
t10 = xs1.*y_dot;
t11 = y.*y_dot;
t12 = y.*ys1;
t13 = y_dot.*ys1;
t14 = z.*z_dot;
t15 = z_dot.*zs1;
t16 = x.*2.0;
t17 = x.^2;
t18 = xs1.*2.0;
t19 = xs1.^2;
t20 = y.*2.0;
t21 = y.^2;
t22 = ys1.*2.0;
t23 = ys1.^2;
t24 = z.*2.0;
t25 = z.^2;
t26 = zs1.*2.0;
t27 = zs1.^2;
t39 = z.*zs1.*-2.0;
t28 = t5.*2.0;
t29 = t7.*2.0;
t30 = t8.*2.0;
t31 = t12.*2.0;
t32 = t24.*zs1;
t33 = -t18;
t34 = -t20;
t35 = -t26;
t36 = -t8;
t38 = -t10;
t40 = -t15;
t41 = t2.*t16;
t42 = t2.*t18;
t43 = t2.*t20;
t44 = t2.*t22;
t45 = t3.*t16;
t46 = t3.*t18;
t47 = t3.*t20;
t48 = t3.*t22;
t49 = t2.*xs1.*-2.0;
t50 = t2.*y.*-2.0;
t51 = t5+t12;
t52 = t6+t13;
t37 = -t30;
t53 = t24+t35;
t54 = t7+t36;
t55 = t28+t31;
t56 = t9+t38;
t58 = t2.*t52;
t59 = t3.*t51.*theta_dot;
t66 = t33+t41+t47;
t67 = t16+t48+t49;
t68 = t22+t45+t50;
t69 = t34+t44+t46;
t57 = t29+t37;
t60 = t2.*t55;
t61 = -t58;
t62 = t3.*t56;
t64 = t2.*t54.*theta_dot;
t63 = t3.*t57;
t65 = -t60;
t74 = t4+t11+t14+t40+t59+t61+t62+t64;
t70 = t17+t19+t21+t23+t25+t27+t39+t63+t65;
t71 = 1.0./sqrt(t70);
t72 = t71.^3;
t73 = t71.*z_dot;
t75 = (t53.*t71)./2.0;
t76 = (t53.*t72.*t74)./2.0;
H_tilde = reshape([(t67.*t71)./2.0,t71.*(x_dot+t3.*theta_dot.*xs1+t2.*theta_dot.*ys1)-(t67.*t72.*t74)./2.0,t69.*t71.*(-1.0./2.0),t71.*(y_dot-t2.*theta_dot.*xs1+t3.*theta_dot.*ys1)+(t69.*t72.*t74)./2.0,t75,t73-t76,0.0,t71.*(x-t2.*xs1+t3.*ys1),0.0,-t71.*(-y+t3.*xs1+t2.*ys1),0.0,t71.*(z-zs1),0.0,0.0,0.0,0.0,0.0,0.0,t66.*t71.*(-1.0./2.0),-t71.*(t2.*x_dot+t3.*y_dot-t3.*theta_dot.*x+t2.*theta_dot.*y)+(t66.*t72.*t74)./2.0,(t68.*t71)./2.0,t71.*(t3.*x_dot-t2.*y_dot+t2.*theta_dot.*x+t3.*theta_dot.*y)-(t68.*t72.*t74)./2.0,-t75,-t73+t76,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[2,18]);
