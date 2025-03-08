function [x_k, y_k, z_k] = compute_satellite_position(ephemeris, time_of_transmission)
    % Extract the parameters from ephemeris
    toe = ephemeris(4, 1);
    a_sqrt = ephemeris(3, 4);
    M_0 = ephemeris(2, 4);
    omega = ephemeris(5, 3);
    i_0 = ephemeris(5, 1);
    i_dot = ephemeris(6, 1);
    e = ephemeris(3, 2);
    delta_n = ephemeris(2, 3);
    OMEGA_0 = ephemeris(4, 3);
    OMEGA_dot = ephemeris(5, 4);
    Crs = ephemeris(2, 2);
    Crc = ephemeris(5, 2);
    Cis = ephemeris(4, 4);
    Cic = ephemeris(4, 2);
    Cus = ephemeris(3, 3);
    Cuc = ephemeris(3, 1);

    % Calculate the satellite position
    Mu = 3.986005E14;
    OMEGA_e_dot = 7.2921151467E-5;
    PI = 3.1415926535898;
    a = a_sqrt^2;
    n_0 = sqrt(Mu/(a^3));
    t_k = time_of_transmission - toe;
    n = n_0 + delta_n;
    M_k = M_0 + n * t_k;
    
    % Normalize M_k
    while M_k >= 2.0 * PI || M_k < 0.0 
        if M_k >= 2.0 * PI
            M_k = M_k - 2.0 * PI;
        end
        if M_k < 0.0
            M_k = M_k + 2.0 * PI;
        end
    end

    % Calculate E_k using for loop
    E_k = M_k;
    tolerance = 1e-8;
    max_iter = 1000;

    for iter = 1:max_iter
        f_E = E_k - e * sin(E_k) - M_k;
        f_prime_E = 1 - e * cos(E_k);
        E_new = E_k - f_E / f_prime_E;

        if abs(E_new - E_k) < tolerance
            break;
        end
        E_k = E_new;
    end

    % Normalize E_k
    while E_k >= 2.0 * PI || E_k < 0.0 
        if E_k >= 2.0 * PI
            E_k = E_k - 2.0 * PI;
        end
        if E_k < 0.0
            E_k = E_k + 2.0 * PI;
        end
    end

    v_k = atan2(sqrt(1-e^2) * sin(E_k), cos(E_k) - e);
    
    % Normalize v_k
    while v_k >= 2.0 * PI || v_k < 0.0 
        if v_k >= 2.0 * PI
            v_k = v_k - 2.0 * PI;
        end
        if v_k < 0.0
            v_k = v_k + 2.0 * PI;
        end
    end

    phi_k = v_k + omega;
    delta_u_k = Cus * sin(2*phi_k) + Cuc * cos(2*phi_k);
    delta_r_k = Crs * sin(2*phi_k) + Crc * cos(2*phi_k);
    delta_i_k = Cis * sin(2*phi_k) + Cic * cos(2*phi_k);

    u_k = phi_k + delta_u_k;
    r_k = a * (1 - e * cos(E_k)) + delta_r_k;
    i_k = i_0 + delta_i_k + i_dot * t_k;

    x_k_0 = r_k * cos(u_k);
    y_k_0 = r_k * sin(u_k);

    OMEGA_k = OMEGA_0 + (OMEGA_dot - OMEGA_e_dot) * t_k - OMEGA_e_dot * toe;

    x_k = x_k_0 * cos(OMEGA_k) - y_k_0 * cos(i_k) * sin(OMEGA_k);
    y_k = x_k_0 * sin(OMEGA_k) + y_k_0 * cos(i_k) * cos(OMEGA_k);
    z_k = y_k_0 * sin(i_k);
end