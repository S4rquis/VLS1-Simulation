% Nonlinear rigid body VLS-1 dynamic

% x = [ x_L; y_L; z_L; u; v; w; phi; theta; psi; p; q; r]

function dx = vls6dof_rhs(t,x,u_cmd, param)
    
    %% States % Checked
    pos   = x(1:3);
    Vb    = x(4:6);

    phi   = x(7);
    theta = x(8);
    psi   = x(9);

    omega = x(10:12);


    p = omega(1);
    q = omega(2);
    r = omega(3);

    %% Nozzle command % Checked
    beta_y = u_cmd(1);
    beta_z = u_cmd(2);
    
    % Saturation (limited deflection of +-4 degrees)
    beta_max = deg2rad(4);    
    beta_y = clip(beta_y, -beta_max, beta_max);
    beta_z = clip(beta_z, -beta_max, beta_max);

    %% Parameters % Checked
    mdot = param.mdot;
    m0 = param.m0;

    Ixx     = param.Ixx;
    Iyy     = param.Iyy;
    Izz     = param.Izz;

    Ixx_dot = param.Ixx_dot;
    Iyy_dot = param.Iyy_dot;
    Izz_dot = param.Izz_dot;

    Tvac    = param.Tvac;

    I    = diag([Ixx Iyy Izz]);
    Idot = diag([Ixx_dot Iyy_dot Izz_dot]);

    %% Naive mass variation description % Checked
    m = m0 + mdot*t;
    %% C_body_from_launch (C^{B/L}) 
    C_BL = [ cos(psi)*cos(theta),                                        sin(psi),                    -sin(theta)*cos(psi);
             sin(theta)*sin(phi) - cos(theta)*sin(psi)*cos(phi),         cos(phi)*cos(psi),           cos(theta)*sin(phi) + sin(theta)*sin(psi)*cos(phi);
             sin(theta)*cos(phi) + cos(theta)*sin(psi)*sin(phi),         -cos(psi)*sin(phi),          cos(theta)*cos(phi) - sin(theta)*sin(psi)*sin(phi) ];


    %% Atmosphere
    h = max(pos(1),0);
    [T,a,patm,rho] = atmoscoesa(h);
    %% Wind
    Vwind_L = wind_model(t,pos);
    Vwind_B = C_BL*Vwind_L;

    %% Air relative velocity
    Vrel_B = Vb - Vwind_B;

    ur = Vrel_B(1);
    vr = Vrel_B(2);
    wr = Vrel_B(3);

    Vrel = norm(Vrel_B);
    Vrel_safe = max(Vrel,1e-6);

    alpha = atan2(wr, ur);  
    beta  = atan2(vr, ur);

    q_bar = 0.5*rho*Vrel^2;
    %% Thrust
    T = max(Tvac - patm*param.Ae,0);

    FE = T * [cos(beta_y)*cos(beta_z);
                 -sin(beta_y);
                  sin(beta_z)];

    %% Gravitational force
    
    F_gL = [-m*param.g; 0; 0];  % Gravity in the launch referencial
    Fg = C_BL * F_gL;

    %% Aerodynamic Force

    FA = [-param.Cx0*q_bar*param.Ar;
          -param.CnBeta*q_bar*param.Ar*beta;
          -param.CnAlpha*q_bar*param.Ar*alpha];

    F = FE + Fg + FA;

    %% Momentum

    lc = [param.lcx; 0; 0];
    la = [param.lax; 0; 0];
    re = [param.rex; 0; 0];

    ME = cross(lc,FE);
    MA = cross(la,FA);

    Pstar = q_bar*param.Ar*param.dr^2/(2*Vrel_safe);

    MAA = Pstar*[-param.Clp*p;
                 -param.Cmq*q;
                 -param.Cnr*r];

    MAJ = mdot*cross(re,cross(omega,re));

    M = ME + MA + MAA + MAJ;

    %% Translacional dynamic
    Vb_dot = F/m - cross(omega,Vb);

    %% Rotacional Dynamic % Checked
    omega_dot = I \ (M - Idot*omega - cross(omega,I*omega)); 

    %% Translacional Cinematic
    pos_dot = C_BL' * Vb; 

    %% Angular Cinematic
    euler_dot = euler_yzx_rates(phi,psi,omega);

    %% Output
    dx = zeros(12,1);

    dx(1:3)   = pos_dot;
    dx(4:6)   = Vb_dot;
    dx(7:9)   = euler_dot;
    dx(10:12) = omega_dot;
end