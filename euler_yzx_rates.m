function eulerdot = euler_yzx_rates(phi,psi,omega)

    p = omega(1);
    q = omega(2);
    r = omega(3);

    if abs(cos(psi)) < 1e-5
        warning('Proximo da singularidade da sequencia YZX: psi = +-90 deg');
        cps = sign(cos(psi))*1e-5;
    end

    phidot   = p - cos(phi)*tan(psi)*q + sin(phi)*tan(psi)*r;
    thetadot = (cos(phi)/cos(psi))*q - (sin(phi)/cos(psi))*r;
    psidot   = sin(phi)*q + cos(phi)*r;

    eulerdot = [phidot; thetadot; psidot];
end